//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlrechain_test.c
//
//----------
//
// nlrechain 的单元测试（docs/nlchain_rechain_plan.md §4 第 1–4 步）。
//
// 直接包含 nlrechain.c 以访问其内部函数。
//	(1) 公式层：gap 表的节点、插值截断、外推、bothGap 按 dq+dt 查表，并与按
//	    kent 描述（小表预计算 + 浮点插值后截断）独立写的参考实现逐点比较；
//	    gapCalcCost 的单调性（搜索层下界依赖它）；交叉点和 connectCost 手算。
//	(2) 拆块：随机 edit script，块的碱基对集合等于比对的碱基对集合，块分之和
//	    减去 gap 罚分等于 score_alignment。
//	(3) DP：2-d 树与暴力枚举（对每个 i 枚举全部 j）逐块比较 total 和 pred。
//	(4) 剥链、链内整理：与按 kent 链表写法独立实现的参考解比较。
//	(5) 输出重建：输出比对的碱基对（按次数）等于保留链的块的碱基对；分数等于
//	    score_alignment；X1、X3 的具体例子。
//
//----------

#include "nlrechain.c"
#include "identity_dist.h"
#include "coverage_dist.h"
#include "continuity_dist.h"
#include "output.h"

//----------
//
// 链接用的桩函数：gapped_extend.o 中 gappily_extend_hsps 引用了它们，本测试
// 不会调到
//
//----------

static void not_in_test (const char* name)
	{ suicidef ("nlrechain_test: %s should not be called", name); }

alignel* filter_aligns_by_identity (arg_dont_complain(seq* seq1), arg_dont_complain(seq* seq2),
    arg_dont_complain(alignel* alignList), arg_dont_complain(float minIdentity),
    arg_dont_complain(float maxIdentity))
	{ not_in_test ("filter_aligns_by_identity");  return NULL; }
alignel* filter_aligns_by_match_count (arg_dont_complain(seq* seq1), arg_dont_complain(seq* seq2),
    arg_dont_complain(alignel* alignList), arg_dont_complain(u32 minMatchCount))
	{ not_in_test ("filter_aligns_by_match_count");  return NULL; }
alignel* filter_aligns_by_mismatch_count (arg_dont_complain(seq* seq1), arg_dont_complain(seq* seq2),
    arg_dont_complain(alignel* alignList), arg_dont_complain(s32 maxMismatchCount))
	{ not_in_test ("filter_aligns_by_mismatch_count");  return NULL; }
alignel* filter_aligns_by_coverage (arg_dont_complain(seq* seq1), arg_dont_complain(seq* seq2),
    arg_dont_complain(alignel* alignList), arg_dont_complain(float minCoverage),
    arg_dont_complain(float maxCoverage))
	{ not_in_test ("filter_aligns_by_coverage");  return NULL; }
alignel* filter_aligns_by_continuity (arg_dont_complain(alignel* alignList),
    arg_dont_complain(float minContinuity), arg_dont_complain(float maxContinuity))
	{ not_in_test ("filter_aligns_by_continuity");  return NULL; }
alignel* filter_aligns_by_num_gaps (arg_dont_complain(alignel* alignList),
    arg_dont_complain(s32 maxSeparateGapsCount))
	{ not_in_test ("filter_aligns_by_num_gaps");  return NULL; }
alignel* filter_aligns_by_num_gap_columns (arg_dont_complain(alignel* alignList),
    arg_dont_complain(s32 maxGapColumnsCount))
	{ not_in_test ("filter_aligns_by_num_gap_columns");  return NULL; }
void print_align_list (arg_dont_complain(alignel* alignList))
	{ not_in_test ("print_align_list"); }
void print_align_list_segments (arg_dont_complain(alignel* alignList))
	{ not_in_test ("print_align_list_segments"); }

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

static u64 rngState = 88172645463325252ULL;

static u32 rnd (u32 n)
	{
	rngState ^= rngState << 13;
	rngState ^= rngState >> 7;
	rngState ^= rngState << 17;
	return (u32) (rngState % n);
	}

// test_scores--
//	lastz 默认的 HOXD70 矩阵（手算用 A/A=91、A/C=-114）。

static scoreset* test_scores (void)
	{
	score t[4][4] = { {   91, -114,  -31, -123 },
	                  { -114,  100, -125,  -31 },
	                  {  -31, -125,  100, -114 },
	                  { -123,  -31, -114,   91 } };
	return new_dna_score_set (t, -100, -100, 400, 30);
	}

// make_seq--
//	不分区的序列；v 后面补一个 0，与 lastz 的约定一致。

static seq* make_seq
   (const char*	s)
	{
	seq*	sq = (seq*) zalloc_or_die ("test seq", sizeof(seq));
	unspos	len = (unspos) strlen (s);

	sq->v = (u8*) malloc_or_die ("test seq v", len+1);
	memcpy (sq->v, s, len+1);
	sq->len      = len;
	sq->trueLen  = len;
	sq->startLoc = 1;
	sq->header   = NULL;
	sq->shortHeader = NULL;
	return sq;
	}

static seq* random_seq
   (unspos	len,
	u32		alphabet)			// 1..4：只用前几个碱基，越小平局越多
	{
	char*	s = (char*) malloc_or_die ("test random seq", len+1);
	seq*	sq;
	unspos	k;

	for (k=0 ; k<len ; k++) s[k] = "ACGT"[rnd(alphabet)];
	s[len] = 0;
	sq = make_seq (s);
	free (s);
	return sq;
	}

static void free_seq (seq* sq) { free (sq->v);  free (sq); }

//----------
//
// (1) 公式层
//
//----------

// kent_gap_ref--
//	按 kent gapCalc 的描述独立实现：1..110 预计算（浮点插值后截断为 int），
//	111 起在表内浮点插值后截断，≥ 表尾按最后一段斜率外推后截断。

static s64 kent_interp
   (const unspos*	pos,
	const s64*		val,
	u32				n,
	u64				x)
	{
	u32		k;
	double	ds, dv;

	for (k=0 ; k<n ; k++)
		{
		if (x == pos[k]) return val[k];
		if (x < pos[k])
			{
			ds = (double) (pos[k] - pos[k-1]);
			dv = (double) (val[k] - val[k-1]);
			return (s64) (val[k-1] + dv * (double) (x - pos[k-1]) / ds);
			}
		}
	ds = (double) (pos[n-1] - pos[n-2]);
	dv = (double) (val[n-1] - val[n-2]);
	return (s64) (val[n-2] + dv * (double) (x - pos[n-2]) / ds);
	}

static s64 kent_row_ref
   (const nlrgaptable*	t,
	const s64*			val,
	u64					d)
	{
	u32		last = t->len - 1;
	double	slope;

	if (d == 0) return 0;
	if (d < 111) return kent_interp (t->position, val, t->len, d);
	if (d >= t->position[last])
		{
		slope = ((double) (val[last] - val[last-1]))
		      / ((double) (t->position[last] - t->position[last-1]));
		return (s64) (val[last] + slope * (double) (d - t->position[last]));
		}
	return kent_interp (t->position, val, t->len, d);
	}

static s64 kent_gap_ref
   (const nlrgaptable*	t,
	s64					dq,
	s64					dt)
	{
	if (dq < 0) dq = 0;
	if (dt < 0) dt = 0;
	if (dt == 0) return kent_row_ref (t, t->qGap, (u64) dq);
	if (dq == 0) return kent_row_ref (t, t->tGap, (u64) dt);
	return kent_row_ref (t, t->bothGap, (u64) (dq + dt));
	}

static void test_gap_table (void)
	{
	const nlrgaptable* tables[2] = { &nlrGapCostMedium, &nlrGapCostLoose };
	const nlrgaptable* t;
	u32		k, w;
	s64		dq, dt, d, prev;

	// 11 个节点：三行都等于表值

	for (w=0 ; w<2 ; w++)
		{
		t = tables[w];
		for (k=0 ; k<t->len ; k++)
			{
			d = (s64) t->position[k];
			check (gap_calc_cost (t, d, 0) == t->qGap[k], "qGap node %u", k);
			check (gap_calc_cost (t, 0, d) == t->tGap[k], "tGap node %u", k);
			if (d >= 2)
				check (gap_calc_cost (t, 1, d-1) == t->bothGap[k], "bothGap node %u", k);
			}
		}

	// 手算：插值截断、外推、bothGap 按 dq+dt、负值归零

	t = &nlrGapCostMedium;
	check (gap_calc_cost (t, 0, 0)   == 0,     "medium (0,0)");
	check (gap_calc_cost (t, 5, 0)   == 487,   "medium q5 = 450+150*2/8");
	check (gap_calc_cost (t, 50, 0)  == 717,   "medium q50 = 600+300*39/100");
	check (gap_calc_cost (t, 112, 0) == 901,   "medium q112 = 900+2000/2000");
	check (gap_calc_cost (t, 0, 7)   == 525,   "medium t7 = 450+150*4/8");
	check (gap_calc_cost (t, 2, 3)   == 887,   "medium both5 = 850+150*2/8");
	check (gap_calc_cost (t, -5, 7)  == 525,   "medium (-5,7) -> t7");
	check (gap_calc_cost (t, -3, -2) == 0,     "medium (-3,-2) -> 0");
	check (gap_calc_cost (t, 300000, 0) == 365789, "medium q300000 外推斜率 1");
	check (gap_calc_cost (t, 150000, 150000) == 318300 + 47889, "medium both300000");
	t = &nlrGapCostLoose;
	check (gap_calc_cost (t, 300000, 0) == 68572,  "loose q300000 = 56600+0.25*47889 截断");
	check (gap_calc_cost (t, 1, 1)      == 660,    "loose both2");
	check (gap_calc_cost (t, 7, 0)      == 425,    "loose q7 = 400+50*4/8");

	// 与 kent 描述的参考实现逐点比较，并检查单调性

	for (w=0 ; w<2 ; w++)
		{
		t = tables[w];
		for (dq=0 ; dq<=400 ; dq++)
			for (dt=0 ; dt<=400 ; dt++)
				{
				check (gap_calc_cost (t, dq, dt) == kent_gap_ref (t, dq, dt),
				       "ref (%lld,%lld)", (long long) dq, (long long) dt);
				if (dq > 0)
					check (gap_calc_cost (t, dq, dt) >= gap_calc_cost (t, dq-1, dt),
					       "monotone dq (%lld,%lld)", (long long) dq, (long long) dt);
				if (dt > 0)
					check (gap_calc_cost (t, dq, dt) >= gap_calc_cost (t, dq, dt-1),
					       "monotone dt (%lld,%lld)", (long long) dq, (long long) dt);
				}
		for (k=0 ; k<20000 ; k++)
			{
			dq = (s64) rnd (2) * (s64) rnd (400000);
			dt = (s64) rnd (2) * (s64) rnd (400000);
			check (gap_calc_cost (t, dq, dt) == kent_gap_ref (t, dq, dt),
			       "ref (%lld,%lld)", (long long) dq, (long long) dt);
			check (gap_calc_cost (t, dq+1, dt) >= gap_calc_cost (t, dq, dt), "monotone dq");
			check (gap_calc_cost (t, dq, dt+1) >= gap_calc_cost (t, dq, dt), "monotone dt");
			}
		prev = 0;
		for (d=1 ; d<400000 ; d+=7)
			{
			check (gap_calc_cost (t, d, 0) >= prev, "row monotone %lld", (long long) d);
			prev = gap_calc_cost (t, d, 0);
			}
		}
	}

static nlrspan make_span (unspos t, unspos q, unspos len)
	{
	nlrspan s;
	s.tStart = t;  s.tEnd = t + len;
	s.qStart = q;  s.qEnd = q + len;
	return s;
	}

static void test_crossover (void)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq = make_seq ("GGGGAAGGGGGG");
	seq*		qSeq;
	nlrctx		c;
	nlrspan		a, b;
	unspos		pos;
	double		adj;

	c.target  = tSeq;
	c.scoring = ss;
	c.gapCost = &nlrGapCostMedium;

	// a = t[0,6) q[0,6)，b = t[4,10) q[10,16)：dq = 4，dt = -2，ov = 2。
	// 重叠列：左 (t4,q4) (t5,q5)，右 (t4,q10) (t5,q11)。
	// 左 A/A=91 A/C=-114，右 A/C=-114 A/A=91：
	//   pos0 = r0+r1 = -23，pos1 = l0+r1 = 182，pos2 = l0+l1 = -23
	//   → pos=1，调整量 = -23 + -23 - 182 = -228（为负）

	qSeq = make_seq ("GGGGACGGGGCAGGGG");
	c.query = qSeq;
	a = make_span (0, 0, 6);
	b = make_span (4, 10, 6);
	adj = find_crossover (&c, &a, &b, 2, &pos);
	check (pos == 1,     "crossover pos %u", (unsigned) pos);
	check (adj == -228,  "crossover adjust %.0f", adj);
	// connectCost = -228 + qGap(dq+ov = 6, dt+ov = 0) = -228 + 506
	check (connect_cost (&c, &a, &b) == -228 + 506, "connect cost %.0f",
	       connect_cost (&c, &a, &b));
	free_seq (qSeq);

	// 平局：左 A/A A/C，右 A/A A/C：三个位置都是 -23，取靠前的 pos=0，
	// 调整量 = -23 + -23 - (-23) = -23

	qSeq = make_seq ("GGGGACGGGGACGGGG");
	c.query = qSeq;
	adj = find_crossover (&c, &a, &b, 2, &pos);
	check (pos == 0,    "tie crossover pos %u", (unsigned) pos);
	check (adj == -23,  "tie crossover adjust %.0f", adj);
	free_seq (qSeq);

	// 包含：ov ≥ 某块长度时加 1e8，gap 用归零后的 dq、dt

	qSeq = make_seq ("GGGGACGGGGACGGGGGGGGGGGG");
	c.query = qSeq;
	a = make_span (0, 0, 6);
	b = make_span (2, 10, 3);				// dt = -4，ov = 4 ≥ |b| = 3
	check (connect_cost (&c, &a, &b) == enclosedCost + 468,		// dq = 4 → qGap(4)
	       "enclosed %.0f", connect_cost (&c, &a, &b));
	free_seq (qSeq);

	// 不重叠：dq=3, dt=5 → bothGap(8)

	qSeq = make_seq ("GGGGACGGGGACGGGGGGGGGGGG");
	c.query = qSeq;
	a = make_span (0, 0, 2);
	b = make_span (7, 5, 2);
	check (connect_cost (&c, &a, &b) == (double) gap_calc_cost (&nlrGapCostMedium, 3, 5),
	       "plain gap");
	check (gap_calc_cost (&nlrGapCostMedium, 3, 5) == 850 + 150*5/8, "bothGap(8)");
	free_seq (qSeq);

	free_seq (tSeq);
	free (ss);
	}

//----------
//
// (2) 拆块
//
//----------

// random_alignment--
//	在 [t0, t0+...) × [q0, q0+...) 上随机生成一条比对：sub 段与 ins/del 交替，
//	允许 ins 紧接 del、sub 被拆成多个连续 sub 操作。首尾都是 sub。

static alignel* random_alignment
   (unspos	t0,
	unspos	q0,
	u32		numRuns,
	u32		maxRun,
	u32		maxGap)
	{
	alignel*	a = (alignel*) zalloc_or_die ("test alignel", sizeof(alignel));
	editscript*	s = edit_script_new ();
	unspos		i = 0, j = 0, run, g;
	u32			k;

	for (k=0 ; k<numRuns ; k++)
		{
		if (k > 0)
			{
			switch (rnd (3))
				{
				case 0:  g = 1 + rnd (maxGap);  edit_script_ins (&s, g);  j += g;  break;
				case 1:  g = 1 + rnd (maxGap);  edit_script_del (&s, g);  i += g;  break;
				default: g = 1 + rnd (maxGap);  edit_script_ins (&s, g);  j += g;
				         g = 1 + rnd (maxGap);  edit_script_del (&s, g);  i += g;  break;
				}
			}
		run = 1 + rnd (maxRun);
		if ((run > 1) && (rnd (4) == 0))
			{
			// 两个相邻的 sub 操作（edit_script_add 会合并，这里绕过它）
			edit_script_sub (&s, run-1);
			s->tailOp = 0;
			edit_script_sub (&s, 1);
			}
		else
			edit_script_sub (&s, run);
		i += run;  j += run;
		}

	a->beg1   = t0 + 1;
	a->beg2   = q0 + 1;
	a->end1   = t0 + i;
	a->end2   = q0 + j;
	a->script = s;
	a->s      = 0;
	return a;
	}

// alignment_pairs--
//	直接按 script 枚举比对的碱基对 (t,q)，追加到 pairs（每对编码为 t*2^32+q）。

static u32 alignment_pairs
   (const alignel*	a,
	u64*			pairs,
	u32				n)
	{
	unspos	i = a->beg1 - 1, j = a->beg2 - 1;
	u32		k, r, rpt;

	for (k=0 ; k<a->script->len ; k++)
		{
		rpt = edit_op_repeat (a->script->op[k]);
		switch (edit_op_operation (a->script->op[k]))
			{
			case editopSub:
				for (r=0 ; r<rpt ; r++) pairs[n++] = (((u64) (i+r)) << 32) | (u64) (j+r);
				i += rpt;  j += rpt;
				break;
			case editopIns: j += rpt;  break;
			case editopDel: i += rpt;  break;
			}
		}
	return n;
	}

static u32 span_pairs
   (const nlrspan*	b,
	u64*			pairs,
	u32				n)
	{
	unspos r;
	for (r=0 ; r<b->tEnd-b->tStart ; r++)
		pairs[n++] = (((u64) (b->tStart+r)) << 32) | (u64) (b->qStart+r);
	return n;
	}

static int cmp_u64 (const void* va, const void* vb)
	{
	u64 a = *(const u64*) va, b = *(const u64*) vb;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

static int same_pairs
   (u64*	x,
	u32		nx,
	u64*	y,
	u32		ny)
	{
	qsort (x, nx, sizeof(u64), cmp_u64);
	qsort (y, ny, sizeof(u64), cmp_u64);
	return (nx == ny) && ((nx == 0) || (memcmp (x, y, nx * sizeof(u64)) == 0));
	}

static void test_split (u32 trials)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq = random_seq (4000, 4);
	seq*		qSeq = random_seq (4000, 4);
	nlrctx		c;
	alignel*	a;
	nlrblock*	blk;
	u32			numBlocks, blkSize, k, trial, nx, ny;
	u64*		x = (u64*) malloc_or_die ("test pairs x", 8000 * sizeof(u64));
	u64*		y = (u64*) malloc_or_die ("test pairs y", 8000 * sizeof(u64));
	double		sumBlocks, gapPen;
	editop		op;

	c.target = tSeq;  c.query = qSeq;  c.scoring = ss;  c.gapCost = &nlrGapCostMedium;

	for (trial=0 ; trial<trials ; trial++)
		{
		a = random_alignment (rnd (500), rnd (500), 1 + rnd (12), 60, 20);
		blk = NULL;  numBlocks = blkSize = 0;
		split_alignment (tSeq, qSeq, a, 7, &blk, &numBlocks, &blkSize);
		score_blocks (&c, blk, numBlocks);

		nx = alignment_pairs (a, x, 0);
		for (k=ny=0 ; k<numBlocks ; k++) ny = span_pairs (&blk[k].span, y, ny);
		check (same_pairs (x, nx, y, ny), "split pairs (trial %u)", trial);

		sumBlocks = 0;
		for (k=0 ; k<numBlocks ; k++)
			{
			check (blk[k].src == 7, "src");
			check (blk[k].ord == k, "ord");
			check ((k == 0) || (blk[k].span.tStart >= blk[k-1].span.tEnd), "order");
			sumBlocks += blk[k].score;
			check (blk[k].score == span_score (&c, &blk[k].span), "block score");
			}

		gapPen = 0;
		for (k=0 ; k<a->script->len ; k++)
			{
			op = a->script->op[k];
			if (edit_op_operation (op) != editopSub)
				gapPen += ss->gapOpen + edit_op_repeat (op) * ss->gapExtend;
			}
		check (sumBlocks - gapPen
		         == score_alignment (ss, tSeq->v, a->beg1-1, qSeq->v, a->beg2-1, a->script),
		       "block scores - gap penalties = score_alignment (trial %u)", trial);

		free (blk);
		free_align_list (a);
		}

	free (x);  free (y);
	free_seq (tSeq);  free_seq (qSeq);
	free (ss);
	}

//----------
//
// (3) DP 对暴力枚举
//
//----------

static void dp_brute
   (const nlrctx*	c,
	nlrblock*		blk,
	u32				n)
	{
	u32		i, j;
	nlrbest	best;

	for (i=0 ; i<n ; i++)
		{
		best.score = 0;
		best.pred  = noIx;
		for (j=0 ; j<i ; j++)
			if (is_predecessor (&blk[j], &blk[i]))
				offer (&best, chain_candidate (c, blk, j, i), j);
		settle (&blk[i], &best);
		}
	}

// test_dp_rules--
//	暴力参考解与 2-d 树共用 offer()/settle()，这里单独手算它们的规则：
//	候选必须严格大于 0 且严格大于块分本身才接上前驱；同分取下标小的 j。

static void test_dp_rules (void)
	{
	nlrblock	b;
	nlrbest		best;

	// 候选 900 > 0 但低于块分 1000：不接
	best.score = 0;  best.pred = noIx;
	offer (&best, 900, 3);
	b.score = 1000;  settle (&b, &best);
	check ((b.pred == noIx) && (b.total == 1000), "cand below own score");

	// 候选 1200 高于块分 1000：接上
	offer (&best, 1200, 5);
	settle (&b, &best);
	check ((b.pred == 5) && (b.total == 1200), "cand above own score");

	// 块分为负：候选 -10 高于块分 -50 但不大于 0，不接
	best.score = 0;  best.pred = noIx;
	offer (&best, -10, 2);
	offer (&best, 0, 1);
	b.score = -50;  settle (&b, &best);
	check ((b.pred == noIx) && (b.total == -50), "cand not above zero");

	// 块分为负：候选 5 > 0，接上
	offer (&best, 5, 4);
	settle (&b, &best);
	check ((b.pred == 4) && (b.total == 5), "positive cand for negative block");

	// 平局：先给 j=7 再给 j=2，取 2；再给 j=9 同分不变
	best.score = 0;  best.pred = noIx;
	offer (&best, 300, 7);
	offer (&best, 300, 2);
	offer (&best, 300, 9);
	check (best.pred == 2, "tie -> smaller index");
	}

// random_blocks--
//	在 span×span 的范围内随机放 n 个块，按 DP 次序排好、合并完全重复的块、
//	算块分。返回块数。

static u32 random_blocks
   (const nlrctx*	c,
	nlrblock*		blk,
	u32				n,
	unspos			span,
	unspos			maxLen,
	int				diagonal)		// true => 大多沿对角线，模拟真实比对
	{
	u32		k;
	unspos	t, q, len;

	for (k=0 ; k<n ; k++)
		{
		len = 1 + rnd (maxLen);
		t   = rnd (span);
		q   = (diagonal && (rnd (4) != 0))? t + rnd (40) : rnd (span);
		memset (&blk[k], 0, sizeof(nlrblock));
		blk[k].span = make_span (t, q, len);
		blk[k].src  = k;
		}
	qsort (blk, n, sizeof(nlrblock), cmp_blocks);
	n = merge_exact_overlaps (blk, n);
	score_blocks (c, blk, n);
	return n;
	}

static void test_dp (u32 trials)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq;
	seq*		qSeq;
	nlrctx		c;
	nlrblock*	x = (nlrblock*) malloc_or_die ("test blk x", 3000 * sizeof(nlrblock));
	nlrblock*	y = (nlrblock*) malloc_or_die ("test blk y", 3000 * sizeof(nlrblock));
	u32			trial, n, k, size, numPred = 0, numNegAdj = 0;
	unspos		span, maxLen;
	u32			alphabet;

	c.scoring = ss;

	for (trial=0 ; trial<trials ; trial++)
		{
		switch (trial % 4)
			{
			case 0:  size = 1 + rnd (12);   span = 60;    maxLen = 12;  break;
			case 1:  size = 1 + rnd (12);   span = 400;   maxLen = 60;  break;
			case 2:  size = 50 + rnd (200); span = 3000;  maxLen = 80;  break;
			default: size = 500 + rnd (2000); span = 60000; maxLen = 300; break;
			}
		alphabet = 1 + rnd (4);					// 字母表小时平局多
		tSeq = random_seq (span + maxLen + 50, alphabet);
		qSeq = random_seq (span + maxLen + 50, alphabet);
		c.target  = tSeq;
		c.query   = qSeq;
		c.gapCost = (rnd (2) == 0)? &nlrGapCostMedium : &nlrGapCostLoose;

		n = random_blocks (&c, x, size, span, maxLen, (rnd (2) == 0));
		memcpy (y, x, n * sizeof(nlrblock));
		chain_dp (&c, x, n);
		dp_brute (&c, y, n);

		for (k=0 ; k<n ; k++)
			{
			check ((x[k].total == y[k].total) && (x[k].pred == y[k].pred),
			       "dp vs brute (trial %u, block %u): total %.0f/%.0f pred %d/%d",
			       trial, k, x[k].total, y[k].total, (int) x[k].pred, (int) y[k].pred);
			if (x[k].pred != noIx)
				{
				numPred++;
				if (connect_cost (&c, &x[x[k].pred].span, &x[k].span)
				      < (double) gap_calc_cost (c.gapCost,
				            ((s64) x[k].span.qStart) - ((s64) x[x[k].pred].span.qEnd),
				            ((s64) x[k].span.tStart) - ((s64) x[x[k].pred].span.tEnd)))
					numNegAdj++;
				}
			}

		free_seq (tSeq);  free_seq (qSeq);
		}

	check (numPred > 0,   "dp: some predecessors");
	check (numNegAdj > 0, "dp: covered negative crossover adjustment");
	fprintf (stderr, "  dp: %u links, %u with negative crossover adjustment\n",
	         numPred, numNegAdj);

	free (x);  free (y);
	free (ss);
	}

//----------
//
// (4) 剥链、链内整理
//
//----------

// peel_ref--
//	每次选 total 最高（平局取下标小）的未用块，沿前驱回溯。

static u32 peel_ref
   (const nlrblock*	blk,
	u32				n,
	u32*			chainBlk,
	u32*			chainStart)
	{
	char*	hit = (char*) zalloc_or_die ("test hit", n+1);
	u32		used = 0, numChains = 0, len = 0, k, best, b, lo, hi, tmp;

	while (used < n)
		{
		best = noIx;
		for (k=0 ; k<n ; k++)
			if ((!hit[k]) && ((best == noIx) || (blk[k].total > blk[best].total))) best = k;
		chainStart[numChains++] = lo = len;
		for (b=best ; ; )
			{
			hit[b] = true;  used++;
			chainBlk[len++] = b;
			if ((blk[b].pred == noIx) || hit[blk[b].pred]) break;
			b = blk[b].pred;
			}
		for (hi=len-1 ; lo<hi ; lo++,hi--)
			{ tmp = chainBlk[lo];  chainBlk[lo] = chainBlk[hi];  chainBlk[hi] = tmp; }
		}
	chainStart[numChains] = len;
	free (hit);
	return numChains;
	}

// overlaps_ref--
//	按 kent chainRemovePartialOverlaps 的链表写法独立实现（用 next[] 模拟链表）。

static u32 overlaps_ref
   (const nlrctx*	c,
	nlrpiece*		p,
	u32				n)
	{
	u32*	next = (u32*) malloc_or_die ("test next", (n+1) * sizeof(u32));
	u32		head, a, b, k, m, prev;
	s64		dq, dt, ov;
	unspos	x;
	int		trimA, trimB;
	nlrpiece* q;

	for (k=0 ; k<n ; k++) next[k] = (k+1 < n)? k+1 : noIx;
	head = 0;
	for (;;)
		{
		trimA = false;
		a = head;  b = next[a];
		while (b != noIx)
			{
			trimB = false;
			dq = ((s64) p[b].span.qStart) - ((s64) p[a].span.qEnd);
			dt = ((s64) p[b].span.tStart) - ((s64) p[a].span.tEnd);
			if ((dq < 0) || (dt < 0))
				{
				ov = (dq < dt)? -dq : -dt;
				if ((ov >= (s64) (p[a].span.qEnd - p[a].span.qStart))
				 || (ov >= (s64) (p[b].span.qEnd - p[b].span.qStart)))
					trimB = true;
				else
					{
					find_crossover (c, &p[a].span, &p[b].span, (unspos) ov, &x);
					p[b].span.qStart += x;  p[b].span.tStart += x;
					p[a].span.qEnd -= (unspos) ov - x;  p[a].span.tEnd -= (unspos) ov - x;
					if (p[b].span.qEnd <= p[b].span.qStart) trimB = true;
					else if (p[a].span.qEnd <= p[a].span.qStart) trimA = true;
					}
				}
			if (trimA) break;
			if (trimB) { next[a] = next[b];  b = next[a]; }
			else       { a = b;  b = next[b]; }
			}
		if (!trimA) break;
		prev = noIx;					// 删掉所有空块
		for (k=head ; k!=noIx ; k=next[k])
			{
			if (p[k].span.qEnd > p[k].span.qStart) { prev = k;  continue; }
			if (prev == noIx) head = next[k];
			             else next[prev] = next[k];
			}
		}

	q = (nlrpiece*) malloc_or_die ("test q", (n+1) * sizeof(nlrpiece));
	for (m=0,k=head ; k!=noIx ; k=next[k]) q[m++] = p[k];
	memcpy (p, q, m * sizeof(nlrpiece));
	free (q);  free (next);
	return m;
	}

static void test_peel_and_clean (u32 trials)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq;
	seq*		qSeq;
	nlrctx		c;
	nlrblock*	blk = (nlrblock*) malloc_or_die ("test blk", 300 * sizeof(nlrblock));
	u32*		cb1 = (u32*) malloc_or_die ("test cb1", 301 * sizeof(u32));
	u32*		cs1 = (u32*) malloc_or_die ("test cs1", 301 * sizeof(u32));
	u32*		cb2 = (u32*) malloc_or_die ("test cb2", 301 * sizeof(u32));
	u32*		cs2 = (u32*) malloc_or_die ("test cs2", 301 * sizeof(u32));
	nlrpiece*	p1 = (nlrpiece*) malloc_or_die ("test p1", 301 * sizeof(nlrpiece));
	nlrpiece*	p2 = (nlrpiece*) malloc_or_die ("test p2", 301 * sizeof(nlrpiece));
	u32			trial, n, nc1, nc2, ch, k, m1, m2, len, numTrimmed = 0, numDropped = 0;
	int			same;

	c.scoring = ss;
	c.gapCost = &nlrGapCostMedium;

	for (trial=0 ; trial<trials ; trial++)
		{
		tSeq = random_seq (2000, 1 + rnd (4));
		qSeq = random_seq (2000, 1 + rnd (4));
		c.target = tSeq;  c.query = qSeq;

		n = random_blocks (&c, blk, 1 + rnd (60), 600, 80, true);
		chain_dp (&c, blk, n);
		nc1 = peel_chains (blk, n, cb1, cs1);
		nc2 = peel_ref    (blk, n, cb2, cs2);
		same = (nc1 == nc2) && (memcmp (cs1, cs2, (nc1+1) * sizeof(u32)) == 0)
		    && (memcmp (cb1, cb2, n * sizeof(u32)) == 0);
		check (same, "peel vs reference (trial %u)", trial);
		check (cs1[nc1] == n, "every block in exactly one chain");

		// 链内整理：对剥出的链，以及一条把全部块按 DP 次序强行串起来的"链"
		// （DP 不会接上被包含的块，后者用来覆盖包含的情形）

		for (ch=0 ; ch<=nc1 ; ch++)
			{
			if (ch < nc1)
				{
				len = cs1[ch+1] - cs1[ch];
				for (k=0 ; k<len ; k++)
					{ p1[k].blk = cb1[cs1[ch]+k];  p1[k].span = blk[p1[k].blk].span; }
				}
			else
				{
				for (k=len=0 ; k<n ; k++)
					if ((len == 0) || is_predecessor (&blk[p1[len-1].blk], &blk[k]))
						{ p1[len].blk = k;  p1[len].span = blk[k].span;  len++; }
				}
			memcpy (p2, p1, len * sizeof(nlrpiece));
			m1 = remove_partial_overlaps (&c, p1, len);
			m2 = overlaps_ref (&c, p2, len);
			check ((m1 == m2) && (memcmp (p1, p2, m1 * sizeof(nlrpiece)) == 0),
			       "remove_partial_overlaps vs reference (trial %u, chain %u)", trial, ch);
			if (m1 < len) numDropped++;
			for (k=0 ; k<m1 ; k++)
				{
				check (p1[k].span.tEnd > p1[k].span.tStart, "piece non-empty");
				check (p1[k].span.tEnd - p1[k].span.tStart == p1[k].span.qEnd - p1[k].span.qStart,
				       "piece is ungapped");
				if ((p1[k].span.tStart != blk[p1[k].blk].span.tStart)
				 || (p1[k].span.tEnd   != blk[p1[k].blk].span.tEnd)) numTrimmed++;
				if (k > 0)
					check ((p1[k].span.qStart >= p1[k-1].span.qEnd)
					    && (p1[k].span.tStart >= p1[k-1].span.tEnd), "no overlap after cleanup");
				}
			}

		free_seq (tSeq);  free_seq (qSeq);
		}

	check (numTrimmed > 0, "cleanup: covered crossover trimming");
	check (numDropped > 0, "cleanup: covered enclosed-block removal");
	fprintf (stderr, "  cleanup: %u trimmed pieces, %u chains with dropped blocks\n",
	         numTrimmed, numDropped);

	free (blk);  free (cb1);  free (cs1);  free (cb2);  free (cs2);  free (p1);  free (p2);
	free (ss);
	}

// test_remove_overlap_cases--
//	手工构造：部分重叠被交叉点裁开；完全包含的 b 被删掉。

static void test_remove_overlap_cases (void)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq = make_seq ("GGGGAAGGGGGGGGGGGGGG");
	seq*		qSeq = make_seq ("GGGGACGGGGCAGGGGGGGG");
	nlrctx		c;
	nlrpiece	p[3];
	u32			n;

	c.target = tSeq;  c.query = qSeq;  c.scoring = ss;  c.gapCost = &nlrGapCostMedium;

	// 同 test_crossover：交叉点 pos=1 → b 起点 +1，a 终点 -(2-1)

	p[0].blk = 0;  p[0].span = make_span (0, 0, 6);
	p[1].blk = 1;  p[1].span = make_span (4, 10, 6);
	n = remove_partial_overlaps (&c, p, 2);
	check (n == 2, "partial: two pieces");
	check ((p[0].span.tEnd == 5) && (p[0].span.qEnd == 5), "partial: a trimmed to 5");
	check ((p[1].span.tStart == 5) && (p[1].span.qStart == 11), "partial: b starts at t5 q11");

	// 包含：b 在 t 上完全落在 a 里，删掉 b；c 照常保留

	p[0].blk = 0;  p[0].span = make_span (0, 0, 10);
	p[1].blk = 1;  p[1].span = make_span (2, 12, 3);
	p[2].blk = 2;  p[2].span = make_span (14, 16, 3);
	n = remove_partial_overlaps (&c, p, 3);
	check (n == 2, "enclosed: b removed");
	check ((p[0].blk == 0) && (p[1].blk == 2), "enclosed: a and c remain");

	free_seq (tSeq);  free_seq (qSeq);
	free (ss);
	}

//----------
//
// (5) 输出重建
//
//----------

// expected_kept_pairs--
//	不经输出重建，直接按第 1–7 步算出保留链的块的碱基对。

static u32 expected_kept_pairs
   (const nlrctx*	c,
	alignel*		list,
	score			minScore,
	u64*			pairs)
	{
	nlrblock*	blk = NULL;
	u32			numBlocks = 0, blkSize = 0, k, src, ch, numChains, len, m, n = 0;
	u32*		cb;
	u32*		cs;
	nlrpiece*	p;
	alignel*	a;

	for (src=0,a=list ; a!=NULL ; src++,a=a->next)
		split_alignment (c->target, c->query, a, src, &blk, &numBlocks, &blkSize);
	qsort (blk, numBlocks, sizeof(nlrblock), cmp_blocks);
	numBlocks = merge_exact_overlaps (blk, numBlocks);
	score_blocks (c, blk, numBlocks);
	cb = (u32*) malloc_or_die ("test cb", (numBlocks+1) * sizeof(u32));
	cs = (u32*) malloc_or_die ("test cs", (numBlocks+1) * sizeof(u32));
	p  = (nlrpiece*) malloc_or_die ("test p", (numBlocks+1) * sizeof(nlrpiece));
	chain_dp (c, blk, numBlocks);
	numChains = peel_chains (blk, numBlocks, cb, cs);
	for (ch=0 ; ch<numChains ; ch++)
		{
		len = cs[ch+1] - cs[ch];
		for (k=0 ; k<len ; k++) { p[k].blk = cb[cs[ch]+k];  p[k].span = blk[p[k].blk].span; }
		m = remove_partial_overlaps (c, p, len);
		if (chain_score (c, p, m) < minScore) continue;
		for (k=0 ; k<m ; k++) n = span_pairs (&p[k].span, pairs, n);
		}
	free (blk);  free (cb);  free (cs);  free (p);
	return n;
	}

static alignel* copy_list (alignel* list)
	{
	alignel*	head = NULL, *tail = NULL, *a, *b;
	for (a=list ; a!=NULL ; a=a->next)
		{
		b = (alignel*) malloc_or_die ("test copy", sizeof(alignel));
		*b = *a;
		b->next = NULL;
		b->script = edit_script_copy (a->script);
		if (head == NULL) head = b; else tail->next = b;
		tail = b;
		}
	return head;
	}

static void test_rebuild (u32 trials)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq;
	seq*		qSeq;
	nlrctx		c;
	nlrechainparams params;
	alignel*	list, *a, *out, *tail;
	u64*		x = (u64*) malloc_or_die ("test x", 200000 * sizeof(u64));
	u64*		y = (u64*) malloc_or_die ("test y", 200000 * sizeof(u64));
	u32			trial, k, numAligns, nx, ny, numOut = 0, numIn = 0;
	unspos		i, j;
	score		minScore;

	c.scoring = ss;
	c.gapCost = &nlrGapCostMedium;
	params.gapCost    = &nlrGapCostMedium;
	params.reportFile = NULL;

	for (trial=0 ; trial<trials ; trial++)
		{
		tSeq = random_seq (6000, 4);
		qSeq = random_seq (6000, 4);
		c.target = tSeq;  c.query = qSeq;

		// 几条大致共线、互相重叠的比对；分数随便，nl_rechain_aligns 不看 a->s

		list = tail = NULL;
		numAligns = 1 + rnd (6);
		for (k=0 ; k<numAligns ; k++)
			{
			i = rnd (2000);
			j = (rnd (3) == 0)? rnd (2000) : i + rnd (300);
			a = random_alignment (i, j, 1 + rnd (15), 120, 30);
			if (list == NULL) list = a; else tail->next = a;
			tail = a;
			}
		minScore = (rnd (3) == 0)? -1000000 : (score) rnd (4000);
		params.minChainScore = minScore;

		nx = expected_kept_pairs (&c, list, minScore, x);
		out = nl_rechain_aligns (tSeq, qSeq, ss, copy_list (list), &params);

		ny = 0;
		for (a=out ; a!=NULL ; a=a->next)
			{
			numOut++;
			ny = alignment_pairs (a, y, ny);
			check (a->s == score_alignment (ss, tSeq->v, a->beg1-1, qSeq->v, a->beg2-1, a->script),
			       "output score");
			check (edit_op_operation (a->script->op[0]) == editopSub, "output starts with sub");
			check (edit_op_operation (a->script->op[a->script->len-1]) == editopSub,
			       "output ends with sub");
			}
		check (same_pairs (x, nx, y, ny), "output pairs = kept chain pairs (trial %u)", trial);
		numIn += numAligns;

		free_align_list (out);
		free_align_list (list);
		free_seq (tSeq);  free_seq (qSeq);
		}
	fprintf (stderr, "  rebuild: %u alignments in, %u out\n", numIn, numOut);

	free (x);  free (y);
	free (ss);
	}

// test_rebuild_cases--
//	具体例子：
//	(a) 单条强比对，下限 0：原样输出（beg/end/script/分数都不变）。
//	(b) X1：两条比对有一个起点相同的块，较短的那块不再输出。
//	(c) X3：输出按来源比对的输入次序排列。
//	(d) 下限高于唯一链分：输出为空。

static void test_rebuild_cases (void)
	{
	scoreset*	ss = test_scores ();
	seq*		tSeq = random_seq (3000, 4);
	seq*		qSeq;
	nlrechainparams params;
	alignel*	a, *b, *out;
	editscript*	s;
	u64*		x = (u64*) malloc_or_die ("test x", 10000 * sizeof(u64));
	u32			nx, k;
	char*		q;

	// query = t[0,200) + t[202,302) + "ACG" + t[302,3000)，再每 97 个碱基改一个；
	// 于是 q[i] = t[i-1]（i ≥ 303），下面的比对都落在正确的对角线上

	q = (char*) malloc_or_die ("test q", 3002);
	memcpy (q,       tSeq->v,       200);
	memcpy (q + 200, tSeq->v + 202, 100);
	memcpy (q + 300, "ACG",         3);
	memcpy (q + 303, tSeq->v + 302, 3000-302+1);
	for (k=0 ; k<3001 ; k+=97) q[k] = (q[k] == 'A')? 'C' : 'A';
	qSeq = make_seq (q);
	free (q);

	params.gapCost    = &nlrGapCostMedium;
	params.reportFile = NULL;
	params.minChainScore = 0;

	// (a) t[100,400) 与 q[100,400) 中间夹一个 2bp 的 del 和一个 3bp 的 ins

	a = (alignel*) zalloc_or_die ("test a", sizeof(alignel));
	s = edit_script_new ();
	edit_script_sub (&s, 100);  edit_script_del (&s, 2);
	edit_script_sub (&s, 100);  edit_script_ins (&s, 3);
	edit_script_sub (&s, 98);
	a->beg1 = 101;  a->end1 = 400;  a->beg2 = 101;  a->end2 = 401;
	a->script = s;  a->hspId = 42;
	a->s = score_alignment (ss, tSeq->v, 100, qSeq->v, 100, s);
	out = nl_rechain_aligns (tSeq, qSeq, ss, copy_list (a), &params);
	check ((out != NULL) && (out->next == NULL), "(a) one alignment");
	if (out != NULL)
		{
		check ((out->beg1 == a->beg1) && (out->end1 == a->end1)
		    && (out->beg2 == a->beg2) && (out->end2 == a->end2), "(a) same extent");
		check ((out->script->len == a->script->len)
		    && (memcmp (out->script->op, a->script->op, a->script->len * sizeof(editop)) == 0),
		       "(a) same script");
		check (out->s == a->s, "(a) same score");
		check (out->hspId == 42, "(a) hspId kept");
		}
	free_align_list (out);

	// (d) 下限高于链分

	params.minChainScore = a->s + 100000;
	out = nl_rechain_aligns (tSeq, qSeq, ss, copy_list (a), &params);
	check (out == NULL, "(d) everything removed");
	params.minChainScore = 0;

	// (b) b 与 a 的第一个块起点相同但更短（t[100,150)），另有一个远处的块
	// t[1149,1449) q[1150,1450)。合并后 b 的第一块并入 a 的第一块，b 只剩远处
	// 那块；碱基对中不再有重复。

	b = (alignel*) zalloc_or_die ("test b", sizeof(alignel));
	s = edit_script_new ();
	edit_script_sub (&s, 50);  edit_script_del (&s, 999);  edit_script_ins (&s, 1000);
	edit_script_sub (&s, 300);
	b->beg1 = 101;  b->end1 = 1449;  b->beg2 = 101;  b->end2 = 1450;
	b->script = s;  b->hspId = 7;
	a->next = b;
	out = nl_rechain_aligns (tSeq, qSeq, ss, copy_list (a), &params);
	nx = 0;
	for (b=out ; b!=NULL ; b=b->next) nx = alignment_pairs (b, x, nx);
	check (nx == 298 + 300, "(b) duplicate block merged away (%u pairs)", nx);

	// (c) 次序：a 的那条在前，b 剩下的那段（hspId 7）在后

	check ((out != NULL) && (out->hspId == 42), "(c) first output from first source");
	check ((out != NULL) && (out->next != NULL) && (out->next->hspId == 7)
	    && (out->next->beg1 == 1150), "(c) second output from second source");
	free_align_list (out);

	free_align_list (a);
	free (x);
	free_seq (tSeq);  free_seq (qSeq);
	free (ss);
	}

//----------
//
// main
//
//----------

int main
   (int		argc,
	char**	argv)
	{
	u32 trials = 3000;

	if (argc > 1) trials = (u32) atoi (argv[1]);

	test_gap_table ();
	test_crossover ();
	test_dp_rules ();
	test_split (trials);
	test_dp (trials);
	test_peel_and_clean (trials);
	test_remove_overlap_cases ();
	test_rebuild (trials);
	test_rebuild_cases ();

	fprintf (stderr, "nlrechain_test: %d checks, %d failures\n", numChecks, numFailures);
	return (numFailures == 0)? EXIT_SUCCESS : EXIT_FAILURE;
	}
