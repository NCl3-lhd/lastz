//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlchain_test.c
//
//----------
//
// nlchain 的单元测试（docs/nlchain_design.md §9 第 1 步）。
//
// 直接包含 nlchain.c 以访问其内部函数。参考解与正式实现共用公式层和转移层，
// 只替换搜索层（§8.4 第 8 条）：
//	(1) 公式层：手算的小例子，以及 gap_cost 在表各节点处的取值。
//	(2) 暴力搜索层：对每个 i 枚举全部 j，与 K-d 树 + 跳转池的结果逐状态比较
//	    （分数、前驱、段内 HSP 数、前驱状态、边类型完全一致）。
//	(3) 穷举：n ≤ 12 时枚举全部合法路径，按平局规则取最优，与 DP 比较。
//	(4) 次级链：主拷贝 + 较弱旁系拷贝；长主路径内部的短片段重复拷贝。
//
//----------

#include "nlchain.c"

//----------
//
// 测试框架
//
//----------

static int numFailures = 0;
static int numChecks   = 0;

#define check(cond,...)                                                      \
	{                                                                        \
	numChecks++;                                                             \
	if (!(cond))                                                             \
		{                                                                    \
		numFailures++;                                                       \
		fprintf (stderr, "FAIL %s:%d: ", __FILE__, __LINE__);                \
		fprintf (stderr, __VA_ARGS__);                                       \
		fprintf (stderr, "\n");                                              \
		}                                                                    \
	}

static int near (double a, double b) { return fabs (a-b) <= 1e-9 * (1+fabs(a)+fabs(b)); }

static u64 rngState = 88172645463325252ULL;
static int equalScores = false;			// true => 所有 HSP 同分且 λ=0，制造大量平局

static u32 rnd (u32 n)
	{
	rngState ^= rngState << 13;
	rngState ^= rngState >> 7;
	rngState ^= rngState << 17;
	return (u32) (rngState % n);
	}

// make_hsp--
//	按 target 区间 [tb, tb+len) 构造 nlhsp。

static nlhsp make_hsp
   (unspos	qb,
	unspos	len,
	unspos	tb,
	int		sigma,
	u32		p,
	double	s)
	{
	nlhsp	h;

	h.qb    = qb;
	h.qe    = qb + len;
	h.tin   = (sigma > 0)? tb : tb + len;
	h.tout  = (sigma > 0)? tb + len : tb;
	h.sigma = sigma;
	h.p     = p;
	h.s     = s;
	h.rho   = s / len;
	return h;
	}

static nlchainparams default_test_params (void)
	{
	nlchainparams p;

	p.gapCost            = &nlGapCostMedium;
	p.antiPen            = 0;
	p.maxGap             = 50000;
	p.maxDiag            = 10000;
	p.overlapFrac        = 0.5;
	p.jumpPen            = 3000;
	p.spanCoef           = 500;
	p.spanScale          = 1000;
	p.spanCap            = 8;
	p.soloScore          = 6000;
	p.minSourceSupport   = 0;
	p.maxJumpScan        = 1000000;
	p.secondaryScore     = 6000;
	p.secondaryPerRegion = 5;
	p.secondaryRatio     = 0.5;
	return p;
	}

// setup_ctx--
//	把 HSP 按 DP 次序排好放进 ctx（单一 query 分区）。

static void setup_ctx
   (nlctx*					c,
	const nlhsp*			h,
	u32						n,
	const nlchainparams*	params,
	nlchainresult*			out)
	{
	nlsortrec*	recs;
	u32			k;

	ctx_alloc (c, n);
	c->params = params;
	c->out    = out;
	out->len  = 0;

	recs = (nlsortrec*) malloc (((n>0)?n:1) * sizeof(nlsortrec));
	for (k=0 ; k<n ; k++)
		{ recs[k].h = h[k];  recs[k].part2 = 0;  recs[k].table = 0;  recs[k].ix = k; }
	qsort (recs, n, sizeof(nlsortrec), sortrec_order);
	for (k=0 ; k<n ; k++)
		{
		c->hsp[k]      = recs[k].h;
		c->part2[k]    = 0;
		c->srcTable[k] = 0;
		c->srcIx[k]    = recs[k].ix;
		}
	free (recs);
	}

//----------
//
// (1) 公式层手算例子
//
//----------

static void test_formulas (void)
	{
	nlchainparams	p = default_test_params ();
	nlhsp			j, i;
	u32				k;
	double			g;

	// gap_cost 在表各节点处取表值；d=0 为 0；节点间线性插值；表外按最后一段斜率外推

	for (k=0 ; k<11 ; k++)
		{
		check (gap_cost (&nlGapCostMedium, nlGapPositions[k]) == nlGapMediumCosts[k],
		       "medium g(%u)", nlGapPositions[k]);
		check (gap_cost (&nlGapCostLoose, nlGapPositions[k]) == nlGapLooseCosts[k],
		       "loose g(%u)", nlGapPositions[k]);
		}
	check (gap_cost (&nlGapCostMedium, 0) == 0, "g(0)");
	check (near (gap_cost (&nlGapCostMedium, 7), 525), "g(7)=%f", gap_cost (&nlGapCostMedium, 7));
	check (near (gap_cost (&nlGapCostMedium, 61), 750), "g(61)");
	check (near (gap_cost (&nlGapCostMedium, 352111), 417900), "g(352111)");
	check (near (gap_cost (&nlGapCostLoose, 352111), 81600), "loose g(352111)");

	// 正链，两端都有间隔：dq=50，dt=60，d=10，g(10)=450+150·7/8=581.25；R=2 => +100

	p.antiPen = 2;
	j = make_hsp (100, 100, 1000, +1, 0, 5000);
	i = make_hsp (250, 100, 1160, +1, 0, 4000);
	check (query_gap (&j,&i) == 50, "dq");
	check (target_gap (&j,&i) == 60, "dt");
	check (edge_gain (&j,&i) == 4000, "gain");
	check (near (normal_cost (&j,&i,&p), 581.25 + 100), "C_N=%f", normal_cost (&j,&i,&p));
	check (is_normal_edge (&j,&i,&p), "normal");
	check (is_jump_edge (&j,&i,&p), "jump");
	p.antiPen = 0;

	// 两端同时重叠，ov ≥ ov_t：ov=20，dt=-15 => 只在 gain 中撤回 ρ_i·20

	i = make_hsp (180, 100, 1085, +1, 0, 4000);
	g = gap_cost (&nlGapCostMedium, 5);
	check (edge_gain (&j,&i) == 4000 - 40*20, "gain ov");
	check (near (normal_cost (&j,&i,&p), g), "C_N ov=%f", normal_cost (&j,&i,&p));
	check (!is_jump_edge (&j,&i,&p), "jump needs ov=0");

	// 两端同时重叠，ov_t > ov：ov=10，ov_t=30 => 总撤回 ρ_i·30（gain 10 + C_N 20）

	i = make_hsp (190, 100, 1070, +1, 0, 4000);
	g = gap_cost (&nlGapCostMedium, 20);
	check (near (g, 627), "g(20)");
	check (near ((i.s - edge_gain (&j,&i)) + (normal_cost (&j,&i,&p) - g), 40*30),
	       "withdraw once");

	// dq = dt = -k：只撤回一次 ρk

	i = make_hsp (190, 100, 1090, +1, 0, 4000);
	check (near (i.s - edge_gain (&j,&i) + normal_cost (&j,&i,&p), 40*10), "dq=dt=-10");

	// 反链：target 沿 query 正向递减

	j = make_hsp (100, 100, 5000, -1, 0, 5000);
	i = make_hsp (250, 100, 4850, -1, 0, 4000);
	check (target_gap (&j,&i) == 50, "dt minus=%f", target_gap (&j,&i));
	check (is_normal_edge (&j,&i,&p), "normal minus");
	check (jump_span (&j,&i) == 0, "D minus");

	// 跳转代价：J + λ·min(L, log2(1 + D/D0))

	j = make_hsp (100, 100, 1000, +1, 0, 5000);
	i = make_hsp (250, 100, 2150, +1, 0, 4000);		// dt=1050，dq=50，D=1000
	check (near (jump_span (&j,&i), 1000), "D same strand");
	i = make_hsp (300, 100, 2200, +1, 0, 4000);		// D=1000
	check (near (jump_cost (&j,&i,&p), 3500), "P_J D=D0");
	i = make_hsp (300, 100, 4200, +1, 0, 4000);		// D=3000
	check (near (jump_cost (&j,&i,&p), 4000), "P_J D=3D0");
	i = make_hsp (300, 100, 1000000, +1, 0, 4000);
	check (near (jump_cost (&j,&i,&p), 3000 + 500*8), "P_J cap");
	i = make_hsp (300, 100, 1500, +1, 1, 4000);
	check (isinf (jump_span (&j,&i)) && near (jump_cost (&j,&i,&p), 7000), "P_J other p");
	check (jump_kind (&j,&i) == 'T', "kind T");
	i = make_hsp (300, 100, 3000, -1, 0, 4000);		// 换向：|tin_i - tout_j| = 3100-1100
	check (near (jump_span (&j,&i), 2000), "D inversion");
	check (jump_kind (&j,&i) == 'I', "kind I");
	i = make_hsp (300, 100, 500, +1, 0, 4000);
	check (jump_kind (&j,&i) == 'D', "kind D");
	i = make_hsp (300, 100, 90000, +1, 0, 4000);
	check (jump_kind (&j,&i) == 'G', "kind G");

	// Normal(i) 边界：maxGap、maxDiag、ov_max、target 起点推进

	p.maxGap = 100;  p.maxDiag = 20;
	j = make_hsp (100, 100, 1000, +1, 0, 5000);
	i = make_hsp (300, 100, 1200, +1, 0, 4000);		// dq=dt=100
	check (is_normal_edge (&j,&i,&p), "dq=maxGap");
	i = make_hsp (301, 100, 1201, +1, 0, 4000);
	check (!is_normal_edge (&j,&i,&p), "dq>maxGap");
	i = make_hsp (250, 100, 1170, +1, 0, 4000);		// d=20
	check (is_normal_edge (&j,&i,&p), "d=maxDiag");
	i = make_hsp (250, 100, 1171, +1, 0, 4000);
	check (!is_normal_edge (&j,&i,&p), "d>maxDiag");
	i = make_hsp (150, 100, 1050, +1, 0, 4000);		// ov=50=f_ov·100
	check (is_normal_edge (&j,&i,&p), "ov=ov_max");
	i = make_hsp (149, 100, 1049, +1, 0, 4000);
	check (!is_normal_edge (&j,&i,&p), "ov>ov_max");
	i = make_hsp (210, 100, 990, +1, 0, 4000);		// tin_i < tin_j
	check (!is_normal_edge (&j,&i,&p), "target must advance");
	}

//----------
//
// (2) 暴力搜索层：与正式搜索层逐状态比较
//
//----------

static void search_brute
   (nlctx*	c,
	u32		i)
	{
	u32		j;

	for (j=c->rangeLo ; j<i ; j++)
		{
		if (is_normal_edge (&c->hsp[j], &c->hsp[i], c->params))
			{
			relax_normal (c, j, i);
			relax_local  (c, j, i);
			}
		if ((is_jump_edge (&c->hsp[j], &c->hsp[i], c->params))
		 && (is_jump_source (c, j)))
			relax_jump (c, j, i);
		}
	}

static int same_path (const nlpath* a, const nlpath* b)
	{
	return (a->score == b->score) && (a->pred == b->pred) && (a->segHsps == b->segHsps)
	    && (a->predState == b->predState) && (a->edge == b->edge);
	}

static void compare_states
   (const nlctx*	a,
	const nlctx*	b,
	const char*		what)
	{
	u32		i;

	for (i=0 ; i<a->n ; i++)
		{
		check (same_path (&a->st[i].U, &b->st[i].U), "%s U[%u] %f/%u vs %f/%u", what, i,
		       a->st[i].U.score, a->st[i].U.pred, b->st[i].U.score, b->st[i].U.pred);
		check (same_path (&a->st[i].R, &b->st[i].R), "%s R[%u] %f/%u vs %f/%u", what, i,
		       a->st[i].R.score, a->st[i].R.pred, b->st[i].R.score, b->st[i].R.pred);
		check (same_path (&a->st[i].L, &b->st[i].L), "%s L[%u] %f/%u vs %f/%u", what, i,
		       a->st[i].L.score, a->st[i].L.pred, b->st[i].L.score, b->st[i].L.pred);
		}
	}

// random_hsps--
//	随机 HSP：部分沿前一个 HSP 共线延续（产生正常边），部分随机放置；长度
//	可超过 maxGap，并含两端同时重叠的情形。

static void random_hsps
   (nlhsp*	h,
	u32		n,
	unspos	qSpan,
	unspos	tSpan,
	unspos	step,
	unspos	maxLen)
	{
	u32		k;
	unspos	len, qb, tb;
	s64		sq, st;
	int		sigma;
	u32		p;
	double	s;

	for (k=0 ; k<n ; k++)
		{
		len = (equalScores)? 20 : 10 + rnd (maxLen);
		s   = (equalScores)? 1500 : 100.0 * (1 + rnd (40));
		if ((k > 0) && (rnd (10) < 6))
			{
			const nlhsp* prev = &h[k-1];
			sigma = prev->sigma;
			p     = prev->p;
			sq    = ((s64) prev->qe) + rnd (step) - step/4;
			if (equalScores) sq = ((s64) prev->qe) + 10*rnd (3);
			qb    = (sq <= (s64) prev->qb)? prev->qb + 1 : (unspos) sq;
			if (equalScores)
				{
				// 格点化：dq、dt 常相等（d=0）或相差 10，制造等值候选
				s64 off = ((s64) qb - (s64) prev->qe) + ((rnd (2) == 0)? 0 : 10*rnd (3));
				if (sigma > 0) st = ((s64) hsp_te (prev)) + off;
				          else st = ((s64) hsp_tb (prev)) - len - off;
				}
			else if (sigma > 0)
				st = ((s64) hsp_te (prev)) + rnd (step) - step/4;
			else
				st = ((s64) hsp_tb (prev)) - len - rnd (step) + step/4;
			tb = (st >= 0)? (unspos) st : rnd (tSpan);
			}
		else
			{
			sigma = (rnd (2) == 0)? +1 : -1;
			p     = (rnd (3) == 0)? 1 : 0;
			qb    = rnd (qSpan);
			tb    = rnd (tSpan);
			if (equalScores) { qb -= qb % 10;  tb -= tb % 10; }
			}
		h[k] = make_hsp (qb, len, tb + p * 10*tSpan, sigma, p, s);
		}
	}

static nlchainparams random_params (void)
	{
	nlchainparams p = default_test_params ();

	p.maxGap           = 30 + rnd (60);
	p.maxDiag          = 10 + rnd (40);
	p.antiPen          = (rnd (2) == 0)? 0 : 1.5;
	p.gapCost          = (rnd (2) == 0)? &nlGapCostMedium : &nlGapCostLoose;
	p.jumpPen          = 200 + rnd (800);
	p.spanCoef         = rnd (300);
	p.spanScale        = 20 + rnd (100);
	p.spanCap          = 4;
	p.soloScore        = 1000 + rnd (3000);
	p.minSourceSupport = (rnd (2) == 0)? 0 : 1000 + rnd (3000);
	if (equalScores)
		{
		// λ=0、J=s、T_solo<s：跳转候选的上界 R[j]+s-J 常恰好等于当前最优
		p.spanCoef  = 0;
		p.antiPen   = 0;
		p.jumpPen   = (rnd (2) == 0)? 1500 : 3000;
		p.soloScore = (rnd (2) == 0)? 1000 : p.soloScore;
		}
	return p;
	}

static void test_search_vs_brute
   (u32		trials,
	u32		n,
	unspos	qSpan,
	unspos	tSpan,
	unspos	step,
	unspos	maxLen)
	{
	nlhsp*			h = (nlhsp*) malloc (n * sizeof(nlhsp));
	nlchainparams	p;
	nlctx			fast, slow;
	nlchainresult	outFast = {0,0,NULL}, outSlow = {0,0,NULL};
	u32				t;

	for (t=0 ; t<trials ; t++)
		{
		equalScores = (t % 3 == 0);
		p = random_params ();
		if (n > 100) { p.maxGap = 50 * step;  p.maxDiag = 10 * step; }
		random_hsps (h, n, qSpan, tSpan, step, maxLen);

		setup_ctx (&fast, h, n, &p, &outFast);
		setup_ctx (&slow, h, n, &p, &outSlow);
		run_dp (&fast, 0, n, search_candidates);
		run_dp (&slow, 0, n, search_brute);
		check (fast.jumpTruncations == 0, "unexpected H_max truncation");
		compare_states (&fast, &slow, "search-vs-brute");
		ctx_free (&fast);
		ctx_free (&slow);
		}

	nl_free_result (&outFast);
	nl_free_result (&outSlow);
	free (h);
	}

//----------
//
// (3) 穷举全部合法路径
//
// 每条路径在每个 HSP 处给出元组 (分数, 段内 HSP 数, 前驱, 前驱状态)；DP 的
// 平局规则等价于从路径末端起按元组逐个比较（§5.3）。
//
//----------

#define maxExhaustive 12

typedef struct extuple
	{
	double		score;
	u32			segHsps;
	u32			pred;
	u8			predState;
	} extuple;

typedef struct expath
	{
	int			valid;
	u32			len;
	extuple		t[maxExhaustive];	// t[0] 为末端
	} expath;

static int tuple_better (const extuple* a, const extuple* b)
	{
	nlpath pa = path_of (a->score, a->pred, a->segHsps, a->predState, 0);
	nlpath pb = path_of (b->score, b->pred, b->segHsps, b->predState, 0);
	return better (&pa, &pb);
	}

// a 严格优于 b？

static int expath_better (const expath* a, const expath* b)
	{
	u32 k;

	if (!b->valid) return a->valid;
	if (!a->valid) return false;
	for (k=0 ; (k<a->len)&&(k<b->len) ; k++)
		{
		if (tuple_better (&a->t[k], &b->t[k])) return true;
		if (tuple_better (&b->t[k], &a->t[k])) return false;
		}
	return false;
	}

typedef struct exsearch
	{
	const nlctx*	c;
	expath			bestU[maxExhaustive], bestR[maxExhaustive], bestL[maxExhaustive];
	u32				node[maxExhaustive];
	extuple			tup[maxExhaustive];
	u8				state[maxExhaustive];
	u32				depth;
	} exsearch;

static void ex_record (exsearch* x, expath* best)
	{
	expath	cand;
	u32		k;

	cand.valid = true;
	cand.len   = x->depth;
	for (k=0 ; k<x->depth ; k++) cand.t[k] = x->tup[x->depth-1-k];
	if (expath_better (&cand, best)) *best = cand;
	}

// 主路径（U/R）：沿正常边或跳转边延伸

static void ex_extend (exsearch* x)
	{
	const nlctx*	c  = x->c;
	u32				d  = x->depth;
	u32				j  = x->node[d-1];
	u8				sj = x->state[d-1];
	double			prev = x->tup[d-1].score;
	u32				i;
	const nlhsp*	hj = &c->hsp[j], *hi;
	double			gain, cost;

	if (sj == stateU) ex_record (x, &x->bestU[j]);
	             else ex_record (x, &x->bestR[j]);

	for (i=j+1 ; i<c->n ; i++)
		{
		hi = &c->hsp[i];
		if (is_normal_edge (hj, hi, c->params))
			{
			gain = edge_gain (hj, hi);
			cost = normal_cost (hj, hi, c->params);
			x->node[d]  = i;
			x->state[d] = stateR;
			x->tup[d].score     = prev + gain - cost;
			x->tup[d].segHsps   = x->tup[d-1].segHsps + 1;
			x->tup[d].pred      = j;
			x->tup[d].predState = sj;
			x->depth++;  ex_extend (x);  x->depth--;
			}

		// E5：跳出只能从 R，且 L[j] ≥ T_src（L 用穷举值）

		if ((sj == stateR)
		 && (is_jump_edge (hj, hi, c->params))
		 && ((c->params->minSourceSupport == 0)
		  || (x->bestL[j].t[0].score >= c->params->minSourceSupport)))
			{
			x->node[d]  = i;
			x->state[d] = (hi->s >= c->params->soloScore)? stateR : stateU;
			x->tup[d].score     = prev + edge_gain (hj, hi) - jump_cost (hj, hi, c->params);
			x->tup[d].segHsps   = 1;
			x->tup[d].pred      = j;
			x->tup[d].predState = stateR;
			x->depth++;  ex_extend (x);  x->depth--;
			}
		}
	}

// 局部共线链（L）：只走正常边

static void ex_extend_local (exsearch* x)
	{
	const nlctx*	c  = x->c;
	u32				d  = x->depth;
	u32				j  = x->node[d-1];
	double			prev = x->tup[d-1].score;
	u32				i;
	double			gain, cost;

	ex_record (x, &x->bestL[j]);

	for (i=j+1 ; i<c->n ; i++)
		{
		if (!is_normal_edge (&c->hsp[j], &c->hsp[i], c->params)) continue;
		gain = edge_gain (&c->hsp[j], &c->hsp[i]);
		cost = normal_cost (&c->hsp[j], &c->hsp[i], c->params);
		x->node[d] = i;
		x->tup[d].score     = prev + gain - cost;
		x->tup[d].segHsps   = x->tup[d-1].segHsps + 1;
		x->tup[d].pred      = j;
		x->tup[d].predState = stateL;
		x->depth++;  ex_extend_local (x);  x->depth--;
		}
	}

// DP 状态回溯成元组序列

static expath dp_expath (const nlctx* c, u32 i, u8 state)
	{
	expath			e;
	const nlpath*	path;

	e.len = 0;
	path  = (state == stateU)? &c->st[i].U : (state == stateR)? &c->st[i].R : &c->st[i].L;
	e.valid = (path->score != unreached);
	while (e.valid)
		{
		e.t[e.len].score     = path->score;
		e.t[e.len].segHsps   = path->segHsps;
		e.t[e.len].pred      = path->pred;
		e.t[e.len].predState = path->predState;
		e.len++;
		if (path->pred == noPred) break;
		i     = path->pred;
		state = path->predState;
		path  = (state == stateU)? &c->st[i].U : (state == stateR)? &c->st[i].R : &c->st[i].L;
		}
	return e;
	}

static int same_expath (const expath* a, const expath* b)
	{
	u32 k;

	if (a->valid != b->valid) return false;
	if (!a->valid) return true;
	if (a->len != b->len) return false;
	for (k=0 ; k<a->len ; k++)
		{
		if (a->t[k].score     != b->t[k].score)     return false;
		if (a->t[k].segHsps   != b->t[k].segHsps)   return false;
		if (a->t[k].pred      != b->t[k].pred)      return false;
		if (a->t[k].predState != b->t[k].predState) return false;
		}
	return true;
	}

static u64 coverNormal = 0, coverJumpU = 0, coverJumpR = 0, coverLocal = 0, coverNoSource = 0;

static void test_exhaustive (u32 trials)
	{
	nlhsp			h[maxExhaustive];
	nlchainparams	p;
	nlctx			c;
	nlchainresult	out = {0,0,NULL};
	exsearch*		x = (exsearch*) malloc (sizeof(exsearch));
	expath			dp;
	u32				t, n, i;

	for (t=0 ; t<trials ; t++)
		{
		n = 2 + rnd (maxExhaustive-1);
		equalScores = (t % 3 == 0);
		p = random_params ();
		random_hsps (h, n, 300, 400, 40, 70);
		setup_ctx (&c, h, n, &p, &out);
		run_dp (&c, 0, n, search_candidates);

		memset (x, 0, sizeof(exsearch));
		x->c = &c;

		// 先穷举 L（E5 要用），再穷举 U/R

		for (i=0 ; i<n ; i++)
			{
			x->node[0] = i;
			x->tup[0].score     = c.hsp[i].s;
			x->tup[0].segHsps   = 1;
			x->tup[0].pred      = noPred;
			x->tup[0].predState = stateStart;
			x->depth = 1;
			ex_extend_local (x);
			}
		for (i=0 ; i<n ; i++)
			{
			x->node[0]  = i;
			x->state[0] = (c.hsp[i].s >= p.soloScore)? stateR : stateU;
			x->tup[0].score     = c.hsp[i].s;
			x->tup[0].segHsps   = 1;
			x->tup[0].pred      = noPred;
			x->tup[0].predState = stateStart;
			x->depth = 1;
			ex_extend (x);
			}

		for (i=0 ; i<n ; i++)
			{
			if (c.st[i].R.edge == edgeNormal) coverNormal++;
			if (c.st[i].U.edge == edgeJump)   coverJumpU++;
			if (c.st[i].R.edge == edgeJump)   coverJumpR++;
			if (c.st[i].L.pred != noPred)     coverLocal++;
			if ((c.st[i].R.score != unreached) && (!is_jump_source (&c, i))) coverNoSource++;
			dp = dp_expath (&c, i, stateU);
			check (same_expath (&dp, &x->bestU[i]), "exhaustive U[%u] trial %u: %f vs %f",
			       i, t, (dp.valid)? dp.t[0].score : -1.0,
			       (x->bestU[i].valid)? x->bestU[i].t[0].score : -1.0);
			dp = dp_expath (&c, i, stateR);
			check (same_expath (&dp, &x->bestR[i]), "exhaustive R[%u] trial %u: %f vs %f",
			       i, t, (dp.valid)? dp.t[0].score : -1.0,
			       (x->bestR[i].valid)? x->bestR[i].t[0].score : -1.0);
			dp = dp_expath (&c, i, stateL);
			check (same_expath (&dp, &x->bestL[i]), "exhaustive L[%u] trial %u", i, t);
			}

		ctx_free (&c);
		}

	nl_free_result (&out);
	free (x);
	}

//----------
//
// H_max 截断计数
//
//----------

static void test_truncation (void)
	{
	nlchainparams	p = default_test_params ();
	nlhsp			h[40];
	nlctx			c;
	nlchainresult	out = {0,0,NULL};
	u32				k;

	// 一段 query 在 target 上命中很多位点：前 39 个 HSP 都能跳入第 40 个

	for (k=0 ; k<39 ; k++)
		h[k] = make_hsp (k*200, 100, 1000000*k, +1, 0, 8000);
	h[39] = make_hsp (39*200 + 1000, 100, 50, +1, 0, 8000);
	p.maxJumpScan = 3;
	setup_ctx (&c, h, 40, &p, &out);
	run_dp (&c, 0, 40, search_candidates);
	check (c.jumpTruncations > 0, "H_max truncation counted");
	ctx_free (&c);
	nl_free_result (&out);
	}

//----------
//
// (4) 回溯与次级链
//
//----------

static void run_all
   (nlctx*					c,
	const nlhsp*			h,
	u32						n,
	const nlchainparams*	p,
	nlchainresult*			out)
	{
	setup_ctx (c, h, n, p, out);
	run_dp (c, 0, n, search_candidates);
	trace_primary (c, 0, n, 0);
	trace_secondary (c, 0, n, 0);
	}

static void test_paralog (void)
	{
	nlchainparams	p = default_test_params ();
	nlhsp			h[12];
	nlctx			c;
	nlchainresult	out = {0,0,NULL};
	u32				k;

	// 主拷贝：target 1000 起，6 个共线 HSP；旁系拷贝：target 500000 起，同一段
	// query，分数较低

	for (k=0 ; k<6 ; k++)
		{
		h[k]   = make_hsp (k*300, 200, 1000   + k*300, +1, 0, 9000);
		h[k+6] = make_hsp (k*300, 200, 500000 + k*300, +1, 0, 7000);
		}

	run_all (&c, h, 12, &p, &out);
	check (out.len == 2, "paralog: %u segments", out.len);
	if (out.len == 2)
		{
		check ((!out.seg[0].isSecondary) && (out.seg[0].count == 6)
		    && (out.seg[0].tBeg == 1000) && (out.seg[0].jumpIn == 'S'),
		       "paralog primary");
		check ((out.seg[1].isSecondary) && (out.seg[1].count == 6)
		    && (out.seg[1].tBeg == 500000) && (out.seg[1].chainId == 1),
		       "paralog secondary is one full chain (count=%u)", out.seg[1].count);
		}
	ctx_free (&c);

	// 没有次级链时只有主路径

	p.secondaryScore = 0;
	run_all (&c, h, 12, &p, &out);
	check (out.len == 1, "secondary off");
	ctx_free (&c);
	nl_free_result (&out);
	}

static void test_fragment_in_long_primary (void)
	{
	nlchainparams	p = default_test_params ();
	nlhsp			h[43];
	nlctx			c;
	nlchainresult	out = {0,0,NULL};
	u32				k, numSecondary = 0, secCount = 0;

	// 长主路径：40 个 HSP 覆盖 query 0..40k；其中 query 20k..23k 的 3 个 HSP 在
	// target 另一处有较弱的片段重复拷贝

	for (k=0 ; k<40 ; k++)
		h[k] = make_hsp (k*1000, 500, 100000 + k*1000, +1, 0, 5000);
	for (k=0 ; k<3 ; k++)
		h[40+k] = make_hsp ((20+k)*1000, 500, 900000 + k*1000, +1, 0, 4000);

	run_all (&c, h, 43, &p, &out);
	for (k=0 ; k<out.len ; k++)
		{
		if (out.seg[k].isSecondary) { numSecondary++;  secCount = out.seg[k].count; }
		else check ((out.seg[k].count == 40), "long primary in one segment");
		}
	check ((numSecondary == 1) && (secCount == 3),
	       "fragment copy kept by r_sec (%u chains, %u hsps)", numSecondary, secCount);
	ctx_free (&c);
	nl_free_result (&out);
	}

static void test_rearrangements (void)
	{
	nlchainparams	p = default_test_params ();
	nlhsp			h[9];
	nlctx			c;
	nlchainresult	out = {0,0,NULL};
	u32				k;

	// query 上依次是：正链块、倒位块（反链）、另一个 target 分区的块

	for (k=0 ; k<3 ; k++)
		{
		h[k]   = make_hsp (k*300,        200, 10000 + k*300, +1, 0, 8000);
		h[k+3] = make_hsp (2000 + k*300, 200, 30000 - k*300, -1, 0, 8000);
		h[k+6] = make_hsp (4000 + k*300, 200, 7000000 + k*300, +1, 1, 8000);
		}

	p.secondaryScore = 0;
	run_all (&c, h, 9, &p, &out);
	check (out.len == 3, "rearrangement: %u segments", out.len);
	if (out.len == 3)
		{
		check (out.seg[0].jumpIn == 'S' && out.seg[0].strand == +1, "seg0 S");
		check (out.seg[1].jumpIn == 'I' && out.seg[1].strand == -1
		    && out.seg[1].qBeg == 2000 && out.seg[1].tEnd == 30200, "seg1 I");
		check (out.seg[2].jumpIn == 'T' && out.seg[2].part1 == 1, "seg2 T");
		}
	ctx_free (&c);
	nl_free_result (&out);
	}

//----------
//
// main
//
//----------

int main (int argc, char** argv)
	{
	u32 trials = (argc > 1)? (u32) atoi (argv[1]) : 3000;

	test_formulas ();
	test_search_vs_brute (trials, 12,   300,    400,    40,  70);
	test_search_vs_brute (trials/10, 300, 20000, 30000, 60, 150);
	test_search_vs_brute (5,     5000, 2000000, 3000000, 400, 600);
	test_exhaustive (trials);
	equalScores = false;
	test_truncation ();
	test_paralog ();
	test_fragment_in_long_primary ();
	test_rearrangements ();

	printf ("exhaustive coverage: R<-N %" PRIu64 ", U<-J %" PRIu64 ", R<-J %" PRIu64
	        ", L<-N %" PRIu64 ", R blocked by T_src %" PRIu64 "\n",
	        coverNormal, coverJumpU, coverJumpR, coverLocal, coverNoSource);
	check ((coverNormal > 0) && (coverJumpU > 0) && (coverJumpR > 0)
	    && (coverLocal > 0) && (coverNoSource > 0), "random cases cover all transitions");
	printf ("%d checks, %d failures\n", numChecks, numFailures);
	return (numFailures == 0)? EXIT_SUCCESS : EXIT_FAILURE;
	}
