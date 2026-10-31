//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlchain.c
//
//----------
//
// nlchain--
//	非线性串联。在 gapped 延伸之前，从两条链的 HSP 中选出一条 query 单调的主
//	路径（路径中允许倒位、易位、回跳等异常边），再沿局部共线状态 L 回溯次级链
//	以覆盖旁系同源，链外的 HSP 全部丢弃。
//
// 模型与记号见 docs/nlchain_design.md，下文的 §x.y 均指该文档章节。代码分三层
// （§8.4）：
//	公式层：§4–§5.1 的纯算术函数，不知道 K-d 树、池、批次的存在。
//	转移层：E1–E6 每个方程一个函数；候选值经 offer() 更新，平局只在 better()
//	        里判断。
//	搜索层：K-d 树（正常边，§6.1）、跳转池与待释放堆（跳转边，§6.2），只负责
//	        枚举候选 j 并调用转移层；剪枝条件只出现在这一层。
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <math.h>				// standard C math stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff

#define  nlchain_owner			// (make this the owner of its globals)
#include "nlchain.h"			// interface to this module

//----------
//
// 内置 g(d) 表（§5.1）
//
// 取自 kent/src/lib/gapCalc.c 的 qGap 行：medium 为 originalGapCosts（axtChain
// 用于人–鼠），loose 为 defaultGapCosts（用于人–鸡）。
//
//----------

static const unspos nlGapPositions[11]
	= { 1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111 };
static const double nlGapMediumCosts[11]
	= { 350, 425, 450, 600, 900, 2900, 22900, 57900, 117900, 217900, 317900 };
static const double nlGapLooseCosts[11]
	= { 325, 360, 400, 450, 600, 1100, 3600, 7600, 15600, 31600, 56600 };

const nlgaptable nlGapCostMedium = { 11, nlGapPositions, nlGapMediumCosts };
const nlgaptable nlGapCostLoose  = { 11, nlGapPositions, nlGapLooseCosts  };

//----------
//
// private data structures
//
//----------

// 状态类别；数值即平局规则中前驱状态的固定次序 U < R < 起点（§5.3）。L 只与 L
// 比较，不与 U/R 混排。

enum { stateU = 0, stateR, stateL, stateStart };

// 边类型

enum { edgeStart = 0, edgeNormal, edgeJump };

#define noPred    ((u32) -1)
#define unreached (-INFINITY)

// nlhsp--
//	一个 HSP 的统一坐标（§4），进入 DP 前一次性算好。转移函数里不再出现
//	pos2、rcf_rev 这类原始字段。

typedef struct nlhsp
	{
	unspos		qb, qe;			// 正链 query 上的半开区间
	unspos		tin, tout;		// 沿 query 正向行走时进入/离开 target 的位置
	int			sigma;			// +1 正链，-1 反链
	u32			p;				// target 分区下标
	double		s, rho;			// HSP 分数；ρ = s / len
	} nlhsp;

// nlpath--
//	U、R 或 L 的一个状态（§8.1）。L 也用这个结构，segHsps 为 L 链的 HSP 数，
//	使三种状态共用 better() 的平局规则。

typedef struct nlpath
	{
	double		score;			// U[i] / R[i] / L[i]；未达为 -inf
	u32			pred;			// 前驱 HSP；noPred 表示起点
	u32			segHsps;		// 当前段内 HSP 数（L 为整条局部链的 HSP 数）
	u8			predState;		// 前驱状态：stateU / stateR / stateL / stateStart
	u8			edge;			// 进入 i 的边：edgeStart / edgeNormal / edgeJump
	} nlpath;

// nlstate--
//	一个 HSP 的全部 DP 状态。

typedef struct nlstate
	{
	nlpath		U, R;			// E1–E4
	nlpath		L;				// E6；L.pred 即文档中的 predL
	} nlstate;

// nlkdnode--
//	(p,σ) 一棵的 2-d 树节点，坐标为 j 的出口 (qe, σ·tout)（§6.1）。

typedef struct nlkdnode
	{
	int			isBucket;		// true => 叶子
	u32			loIx, hiIx;		// 叶子：perm[] 中的下标区间（闭区间）
								// 内部节点：hiIx 为 cutVal 对应的下标
	sgnpos		cutVal;			// 内部节点的切分值；lo 子树 ≤ cutVal ≤ hi 子树
	double		maxN;			// 子树内 max(U,R) 的最大值；未计算为 -inf
	double		maxL;			// 子树内 L 的最大值；未计算为 -inf
	struct nlkdnode* loSon, *hiSon;
	} nlkdnode;

#define bucketSize 3

// nlheap--
//	u32 下标的二叉堆；次序由调用者给出的 above() 决定。

typedef struct nlheap
	{
	u32			len;
	u32*		v;
	} nlheap;

// nlbox--
//	K-d 树的 box 查询范围，下标 0 为 qe 轴，1 为 σ·tout 轴。

typedef struct nlbox
	{
	double		lo[2], hi[2];
	} nlbox;

// nlctx--
//	一次 nl_reduce_to_chains 调用的全部工作数据。HSP 已按 query 分区和
//	(qb, p, tb, σ, len) 排好序，下标即平局规则中的"前驱下标"。

typedef struct nlctx
	{
	const nlchainparams* params;
	u32			n;
	nlhsp*		hsp;
	nlstate*	st;
	u32*		part2;			// HSP 所在 query 分区
	u8*			srcTable;		// 0 => stMinus，1 => stPlus
	u32*		srcIx;			// 在原 segtable 中的下标

	// 搜索层

	u32			rangeLo;		// 当前 DP 区间 [rangeLo, rangeHi)
	u32			rangeHi;
	u32*		perm;			// K-d 树的点排列
	u32*		treePos;		// HSP 在 perm[] 中的位置
	nlkdnode**	treeOf;			// HSP 所属 (p,σ) 树的根
	nlkdnode*	nodes;			// 树节点池
	u32			nodesLen;
	nlheap		pending;		// 待释放：按 qe 升序
	nlheap		pool;			// 跳转池：按 R 降序
	nlheap		frontier;		// 按 R 降序遍历 pool 时的前沿（存 pool 中的位置）
	u64			jumpTruncations;// H_max 截断次数

	// 回溯

	u8*			visited;
	u32*		pathIx;			// 主路径 HSP，按 query 升序
	u8*			pathState;		// 主路径上各 HSP 的状态
	u32			pathLen;
	u32*		chainBuf;		// 次级链回溯缓冲
	u32*		cands;			// 次级链候选终点
	u32*		members;		// 输出段的 HSP，按段、段内按路径次序
	u32			membersLen;
	nlchainresult* out;
	} nlctx;

typedef void (*nlsearch) (nlctx* c, u32 i);

//----------
//
// prototypes for private functions
//
//----------

static unspos   hsp_len          (const nlhsp* h);
static unspos   hsp_tb           (const nlhsp* h);
static unspos   hsp_te           (const nlhsp* h);
static double   query_gap        (const nlhsp* j, const nlhsp* i);
static double   target_gap       (const nlhsp* j, const nlhsp* i);
static double   query_overlap    (const nlhsp* j, const nlhsp* i);
static double   edge_gain        (const nlhsp* j, const nlhsp* i);
static double   gap_cost         (const nlgaptable* t, unspos d);
static double   normal_cost      (const nlhsp* j, const nlhsp* i,
                                  const nlchainparams* p);
static double   jump_span        (const nlhsp* j, const nlhsp* i);
static double   jump_cost        (const nlhsp* j, const nlhsp* i,
                                  const nlchainparams* p);
static int      is_advancing     (const nlhsp* j, const nlhsp* i,
                                  const nlchainparams* p);
static int      is_normal_edge   (const nlhsp* j, const nlhsp* i,
                                  const nlchainparams* p);
static int      is_jump_edge     (const nlhsp* j, const nlhsp* i,
                                  const nlchainparams* p);
static char     jump_kind        (const nlhsp* j, const nlhsp* i);

static nlpath   path_of          (double score, u32 pred, u32 segHsps,
                                  u8 predState, u8 edge);
static int      better           (const nlpath* a, const nlpath* b);
static void     offer            (nlpath* acc, const nlpath* cand);
static int      is_jump_source   (const nlctx* c, u32 j);
static void     init_state       (nlctx* c, u32 i);
static void     relax_jump       (nlctx* c, u32 j, u32 i);
static void     relax_normal     (nlctx* c, u32 j, u32 i);
static void     relax_local      (nlctx* c, u32 j, u32 i);
static void     settle_state     (nlctx* c, u32 i);

static void     ctx_alloc        (nlctx* c, u32 n);
static void     ctx_free         (nlctx* c);
static void     heap_push        (nlctx* c, nlheap* h, u32 x,
                                  int (*above) (const nlctx*, u32, u32));
static u32      heap_pop         (nlctx* c, nlheap* h,
                                  int (*above) (const nlctx*, u32, u32));
static void     build_trees      (nlctx* c, u32 a, u32 b);
static void     kd_propagate     (nlctx* c, u32 i);
static void     search_normal    (nlctx* c, const nlkdnode* node, int axis,
                                  u32 i, const nlbox* box);
static void     search_jump      (nlctx* c, u32 i);
static void     search_candidates(nlctx* c, u32 i);
static void     release_pending  (nlctx* c, unspos qb);
static void     run_dp           (nlctx* c, u32 a, u32 b, nlsearch search);

static double   trace_primary    (nlctx* c, u32 a, u32 b, u32 part2);
static void     trace_secondary  (nlctx* c, u32 a, u32 b, u32 part2);
static void     emit_segment     (nlctx* c, const u32* hsps, u32 count,
                                  char jumpIn, int isSecondary, u32 chainId,
                                  u32 part2);

//----------
//
// 公式层（§4、§5.1）
//
//----------

static unspos hsp_len (const nlhsp* h) { return h->qe - h->qb; }
static unspos hsp_tb  (const nlhsp* h) { return (h->sigma > 0)? h->tin  : h->tout; }
static unspos hsp_te  (const nlhsp* h) { return (h->sigma > 0)? h->tout : h->tin;  }

// dq = qb_i - qe_j   (§4)，负数表示 query 重叠

static double query_gap
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	return ((double) i->qb) - ((double) j->qe);
	}

// dt = σ·(tin_i - tout_j)   (§4)，仅同 p、同 σ 时有意义

static double target_gap
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	return i->sigma * (((double) i->tin) - ((double) j->tout));
	}

// ov = max(0, -dq)   (§4)

static double query_overlap
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	double dq = query_gap (j, i);
	return (dq < 0)? -dq : 0;
	}

// gain(j,i) = s_i - ρ_i·ov   (§4)

static double edge_gain
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	return i->s - i->rho * query_overlap (j, i);
	}

// g(d)：代价表上分段线性插值，d = 0 时为 0，表外按最后一段斜率外推   (§5.1)

static double gap_cost
   (const nlgaptable*	t,
	unspos				d)
	{
	u32		k, last;
	double	slope;

	if (d == 0) return 0;
	if (d <= t->position[0]) return t->cost[0];	// （表首为 1，d ≥ 1 时只在 d=1 处取到）

	last = t->len - 1;
	if (d >= t->position[last])
		{
		slope = (t->cost[last] - t->cost[last-1])
		      / (((double) t->position[last]) - ((double) t->position[last-1]));
		return t->cost[last] + slope * (((double) d) - ((double) t->position[last]));
		}

	for (k=1 ; t->position[k]<d ; k++)		// position[k-1] < d ≤ position[k]
		;
	return t->cost[k-1]
	     + (t->cost[k] - t->cost[k-1])
	     * (((double) d) - ((double) t->position[k-1]))
	     / (((double) t->position[k]) - ((double) t->position[k-1]));
	}

// C_N = g(d) + R·max(0,min(dq,dt)) + ρ_i·max(0, ov_t - ov)   (§5.1)

static double normal_cost
   (const nlhsp*			j,
	const nlhsp*			i,
	const nlchainparams*	p)
	{
	double	dq   = query_gap  (j, i);
	double	dt   = target_gap (j, i);
	double	ov   = (dq < 0)? -dq : 0;
	double	ov_t = (dt < 0)? -dt : 0;
	double	d    = fabs (dq - dt);
	double	both = (dq < dt)? dq : dt;

	return gap_cost (p->gapCost, (unspos) d)
	     + p->antiPen * ((both > 0)? both : 0)
	     + i->rho * ((ov_t > ov)? ov_t - ov : 0);
	}

// D(j,i) 的三种情形   (§5.1)
//	同 p、同 σ：|dt - dq|；同 p、换向：|tin_i - tout_j|；不同 p：+∞

static double jump_span
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	if (j->p != i->p)
		return INFINITY;
	if (j->sigma == i->sigma)
		return fabs (target_gap (j, i) - query_gap (j, i));
	return fabs (((double) i->tin) - ((double) j->tout));
	}

// P_J = J + λ·min(L, log2(1 + D/D0))   (§5.1)

static double jump_cost
   (const nlhsp*			j,
	const nlhsp*			i,
	const nlchainparams*	p)
	{
	double	D    = jump_span (j, i);
	double	span = (isinf (D))? p->spanCap
	                          : log2 (1 + D / p->spanScale);

	if (span > p->spanCap) span = p->spanCap;
	return p->jumpPen + p->spanCoef * span;
	}

// 所有边共有的推进约束   (§4)
//	qb_j < qb_i，qe_j < qe_i，ov ≤ f_ov·min(len_i, len_j)

static int is_advancing
   (const nlhsp*			j,
	const nlhsp*			i,
	const nlchainparams*	p)
	{
	unspos	minLen;

	if (j->qb >= i->qb) return false;
	if (j->qe >= i->qe) return false;
	minLen = (hsp_len (i) < hsp_len (j))? hsp_len (i) : hsp_len (j);
	return (query_overlap (j, i) <= p->overlapFrac * minLen);
	}

// Normal(i) 的判定   (§5.1)
//	p_j = p_i，σ_j = σ_i，σ·(tin_i - tin_j) > 0，
//	-ov_max_i ≤ dq ≤ maxGap，-ov_max_i ≤ dt ≤ maxGap，|dq - dt| ≤ maxDiag

static int is_normal_edge
   (const nlhsp*			j,
	const nlhsp*			i,
	const nlchainparams*	p)
	{
	double	dq, dt, ov_max_i;

	if (j->p     != i->p)     return false;
	if (j->sigma != i->sigma) return false;
	if (!is_advancing (j, i, p)) return false;
	if (i->sigma * (((double) i->tin) - ((double) j->tin)) <= 0) return false;

	dq       = query_gap  (j, i);
	dt       = target_gap (j, i);
	ov_max_i = p->overlapFrac * hsp_len (i);
	if ((dq < -ov_max_i) || (dq > p->maxGap)) return false;
	if ((dt < -ov_max_i) || (dt > p->maxGap)) return false;
	return (fabs (dq - dt) <= p->maxDiag);
	}

// 跳转边的判定   (§4、§5.1)：满足推进约束且 query 上不重叠（qe_j ≤ qb_i）

static int is_jump_edge
   (const nlhsp*			j,
	const nlhsp*			i,
	const nlchainparams*	p)
	{
	return (is_advancing (j, i, p)) && (j->qe <= i->qb);
	}

// 跳转边的进入方式   (§7 第 3 步)
//	T 换分区；I 换向；D 同向回跳（dt < 0）；G 同向长间隔

static char jump_kind
   (const nlhsp*	j,
	const nlhsp*	i)
	{
	if (j->p     != i->p)        return 'T';
	if (j->sigma != i->sigma)    return 'I';
	if (target_gap (j, i) < 0)   return 'D';
	return 'G';
	}

//----------
//
// 转移层（§5.3，E1–E6）
//
//----------

static nlpath path_of
   (double	score,
	u32		pred,
	u32		segHsps,
	u8		predState,
	u8		edge)
	{
	nlpath	path;

	path.score     = score;
	path.pred      = pred;
	path.segHsps   = segHsps;
	path.predState = predState;
	path.edge      = edge;
	return path;
	}

// better--
//	平局规则（§5.3），所有状态共用：分数高者优先，其次段内 HSP 数多者，再次
//	前驱下标小者，最后前驱状态按 U < R < 起点。a 严格优于 b 时返回 true。

static int better
   (const nlpath*	a,
	const nlpath*	b)
	{
	if (a->score     != b->score)     return (a->score     > b->score);
	if (a->segHsps   != b->segHsps)   return (a->segHsps   > b->segHsps);
	if (a->pred      != b->pred)      return (a->pred      < b->pred);
	return (a->predState < b->predState);
	}

// offer--
//	用候选更新状态；未达的候选（-inf）不参与。

static void offer
   (nlpath*			acc,
	const nlpath*	cand)
	{
	if (cand->score == unreached) return;
	if (better (cand, acc)) *acc = *cand;
	}

// (E5) Jump(i) 只含 R 已达且 L[j] ≥ T_src 的 j（T_src = 0 时不限制）

static int is_jump_source
   (const nlctx*	c,
	u32				j)
	{
	const nlstate* st = &c->st[j];

	if (st->R.score == unreached) return false;
	if (c->params->minSourceSupport == 0) return true;
	return (st->L.score >= c->params->minSourceSupport);
	}

// (E1) U0[i] = s_i；同时给出 (E6) 的起点项 L[i] ≥ s_i。R 在 E3 之前未达。

static void init_state
   (nlctx*	c,
	u32		i)
	{
	nlstate*	st = &c->st[i];
	double		s  = c->hsp[i].s;

	st->U = path_of (s,         noPred, 1, stateStart, edgeStart);
	st->R = path_of (unreached, noPred, 0, stateStart, edgeStart);
	st->L = path_of (s,         noPred, 1, stateStart, edgeStart);
	}

// (E2) UJ[i] = max_j R[j] + gain(j,i) - P_J(j,i)
//	（UJ 与 U0 在 U 中累积，E4 再定归属）

static void relax_jump
   (nlctx*	c,
	u32		j,
	u32		i)
	{
	const nlhsp*	hj = &c->hsp[j];
	const nlhsp*	hi = &c->hsp[i];
	nlpath			cand;

	cand = path_of (c->st[j].R.score + edge_gain (hj, hi) - jump_cost (hj, hi, c->params),
	                j, 1, stateR, edgeJump);
	offer (&c->st[i].U, &cand);
	}

// (E3) RN[i] = max_j max(U[j],R[j]) + gain(j,i) - C_N(j,i)
//	（max(U[j],R[j]) 拆成两个候选，由 better() 按平局规则取舍）

static void relax_normal
   (nlctx*	c,
	u32		j,
	u32		i)
	{
	const nlhsp*	hj   = &c->hsp[j];
	const nlhsp*	hi   = &c->hsp[i];
	const nlstate*	stj  = &c->st[j];
	double			gain = edge_gain (hj, hi);
	double			cost = normal_cost (hj, hi, c->params);
	nlpath			cand;

	cand = path_of (stj->U.score + gain - cost, j, stj->U.segHsps+1, stateU, edgeNormal);
	offer (&c->st[i].R, &cand);

	cand = path_of (stj->R.score + gain - cost, j, stj->R.segHsps+1, stateR, edgeNormal);
	offer (&c->st[i].R, &cand);
	}

// (E6) L[i] = max( s_i, max_j L[j] + gain(j,i) - C_N(j,i) )
//	（s_i 项由 init_state 给出）

static void relax_local
   (nlctx*	c,
	u32		j,
	u32		i)
	{
	const nlhsp*	hj   = &c->hsp[j];
	const nlhsp*	hi   = &c->hsp[i];
	const nlstate*	stj  = &c->st[j];
	double			gain = edge_gain (hj, hi);
	double			cost = normal_cost (hj, hi, c->params);
	nlpath			cand;

	cand = path_of (stj->L.score + gain - cost, j, stj->L.segHsps+1, stateL, edgeNormal);
	offer (&c->st[i].L, &cand);
	}

// (E4) 若 s_i ≥ T_solo：R[i] = max(RN[i], U0[i], UJ[i])，U[i] = -∞
//      否则：           R[i] = RN[i]，              U[i] = max(U0[i], UJ[i])

static void settle_state
   (nlctx*	c,
	u32		i)
	{
	nlstate*	st = &c->st[i];

	if (c->hsp[i].s >= c->params->soloScore)
		{
		offer (&st->R, &st->U);
		st->U = path_of (unreached, noPred, 0, stateStart, edgeStart);
		}
	}

//----------
//
// 搜索层（§6）
//
//----------

static void ctx_alloc
   (nlctx*	c,
	u32		n)
	{
	size_t	m = (n > 0)? n : 1;

	c->n         = n;
	c->hsp       = (nlhsp*)      malloc_or_die ("nlchain hsp",       m * sizeof(nlhsp));
	c->st        = (nlstate*)    malloc_or_die ("nlchain st",        m * sizeof(nlstate));
	c->part2     = (u32*)        malloc_or_die ("nlchain part2",     m * sizeof(u32));
	c->srcTable  = (u8*)         malloc_or_die ("nlchain srcTable",  m * sizeof(u8));
	c->srcIx     = (u32*)        malloc_or_die ("nlchain srcIx",     m * sizeof(u32));
	c->perm      = (u32*)        malloc_or_die ("nlchain perm",      m * sizeof(u32));
	c->treePos   = (u32*)        malloc_or_die ("nlchain treePos",   m * sizeof(u32));
	c->treeOf    = (nlkdnode**)  malloc_or_die ("nlchain treeOf",    m * sizeof(nlkdnode*));
	c->nodes     = (nlkdnode*)   malloc_or_die ("nlchain nodes",     (2*m+1) * sizeof(nlkdnode));
	c->pending.v = (u32*)        malloc_or_die ("nlchain pending",   m * sizeof(u32));
	c->pool.v    = (u32*)        malloc_or_die ("nlchain pool",      m * sizeof(u32));
	c->frontier.v= (u32*)        malloc_or_die ("nlchain frontier",  m * sizeof(u32));
	c->visited   = (u8*)         zalloc_or_die ("nlchain visited",   m * sizeof(u8));
	c->pathIx    = (u32*)        malloc_or_die ("nlchain pathIx",    m * sizeof(u32));
	c->pathState = (u8*)         malloc_or_die ("nlchain pathState", m * sizeof(u8));
	c->chainBuf  = (u32*)        malloc_or_die ("nlchain chainBuf",  m * sizeof(u32));
	c->cands     = (u32*)        malloc_or_die ("nlchain cands",     m * sizeof(u32));
	c->members   = (u32*)        malloc_or_die ("nlchain members",   m * sizeof(u32));
	c->nodesLen        = 0;
	c->pending.len     = 0;
	c->pool.len        = 0;
	c->frontier.len    = 0;
	c->jumpTruncations = 0;
	c->pathLen         = 0;
	c->membersLen      = 0;
	c->rangeLo         = 0;
	c->rangeHi         = 0;
	}

static void ctx_free
   (nlctx*	c)
	{
	free_if_valid ("nlchain hsp",       c->hsp);
	free_if_valid ("nlchain st",        c->st);
	free_if_valid ("nlchain part2",     c->part2);
	free_if_valid ("nlchain srcTable",  c->srcTable);
	free_if_valid ("nlchain srcIx",     c->srcIx);
	free_if_valid ("nlchain perm",      c->perm);
	free_if_valid ("nlchain treePos",   c->treePos);
	free_if_valid ("nlchain treeOf",    c->treeOf);
	free_if_valid ("nlchain nodes",     c->nodes);
	free_if_valid ("nlchain pending",   c->pending.v);
	free_if_valid ("nlchain pool",      c->pool.v);
	free_if_valid ("nlchain frontier",  c->frontier.v);
	free_if_valid ("nlchain visited",   c->visited);
	free_if_valid ("nlchain pathIx",    c->pathIx);
	free_if_valid ("nlchain pathState", c->pathState);
	free_if_valid ("nlchain chainBuf",  c->chainBuf);
	free_if_valid ("nlchain cands",     c->cands);
	free_if_valid ("nlchain members",   c->members);
	}

//----------
//
// 二叉堆
//
//----------

// 待释放堆：qe 小者在上（同 qe 时下标小者在上）

static int pending_above
   (const nlctx*	c,
	u32				a,
	u32				b)
	{
	if (c->hsp[a].qe != c->hsp[b].qe) return (c->hsp[a].qe < c->hsp[b].qe);
	return (a < b);
	}

// 跳转池：R 大者在上（同 R 时下标小者在上）

static int pool_above
   (const nlctx*	c,
	u32				a,
	u32				b)
	{
	if (c->st[a].R.score != c->st[b].R.score) return (c->st[a].R.score > c->st[b].R.score);
	return (a < b);
	}

// 遍历前沿：元素是 pool.v[] 中的位置，按所指 HSP 的池次序

static int frontier_above
   (const nlctx*	c,
	u32				a,
	u32				b)
	{
	return pool_above (c, c->pool.v[a], c->pool.v[b]);
	}

static void heap_push
   (nlctx*	c,
	nlheap*	h,
	u32		x,
	int		(*above) (const nlctx*, u32, u32))
	{
	u32		k, parent;

	k = h->len++;
	while (k > 0)
		{
		parent = (k-1) / 2;
		if (!above (c, x, h->v[parent])) break;
		h->v[k] = h->v[parent];
		k = parent;
		}
	h->v[k] = x;
	}

static u32 heap_pop
   (nlctx*	c,
	nlheap*	h,
	int		(*above) (const nlctx*, u32, u32))
	{
	u32		top, x, k, child;

	top = h->v[0];
	x   = h->v[--h->len];
	k   = 0;
	while (true)
		{
		child = 2*k + 1;
		if (child >= h->len) break;
		if ((child+1 < h->len) && (above (c, h->v[child+1], h->v[child])))
			child++;
		if (!above (c, h->v[child], x)) break;
		h->v[k] = h->v[child];
		k = child;
		}
	if (h->len > 0) h->v[k] = x;
	return top;
	}

//----------
//
// K-d 树（§6.1）
//
// 建树和划分沿用 chain.c 的 build_kd_tree / partition_segments，改动是：每个
// (p,σ) 一棵、坐标轴为 j 的出口 (qe, σ·tout)、节点存两个 double 最大值、查询
// 为 box 查询。
//
//----------

static sgnpos projection
   (const nlctx*	c,
	u32				k,
	int				axis)
	{
	const nlhsp* h = &c->hsp[c->perm[k]];

	if (axis == 0) return (sgnpos) h->qe;
	return h->sigma * (sgnpos) h->tout;
	}

static void perm_swap
   (nlctx*	c,
	u32		a,
	u32		b)
	{
	u32 t = c->perm[a];  c->perm[a] = c->perm[b];  c->perm[b] = t;
	}

// partition_points--
//	把 perm[lo..hi] 按 axis 坐标划分为 lo..m-1 ≤ m ≤ m+1..hi，返回 m（不会是 hi）。
//	（算法与 chain.c 的 partition_segments 相同）

static u32 partition_points
   (nlctx*	c,
	u32		lo,
	u32		hi,
	int		axis)
	{
	u32		m, i, j;
	sgnpos	a, b, d, pivot;

	if (hi - lo < 2)
		suicidef ("nlchain partition: cannot happen (%u,%u)", lo, hi);

	while (true)
		{
		// 取 lo、中点、hi 三者的中位数作 pivot，移到最前

		m = (lo+hi)/2;
		a = projection (c, lo, axis);
		b = projection (c, m,  axis);
		d = projection (c, hi, axis);

		if (((a <= b) && (b <= d)) || ((d <= b) && (b <= a)))
			{ perm_swap (c,lo,m);  pivot = b; }
		else if (((a <= d) && (d <= b)) || ((b <= d) && (d <= a)))
			{ perm_swap (c,lo,hi);  pivot = d; }
		else
			pivot = a;

		// 小的移到前面，大的移到后面

		i = lo;
		j = hi+1;
		while (i < j)
			{
			for (i++ ; (i<=hi)&&(projection(c,i,axis)<=pivot) ; i++)
				;
			for (j-- ; (j>=lo)&&(projection(c,j,axis)>pivot) ; j--)
				;
			perm_swap (c,i,j);
			}

		perm_swap (c,i,j);		// 撤销最后一次交换
		perm_swap (c,lo,j);		// pivot 放到最终位置

		// 不能返回 j==hi（建树会无限递归）；j==hi 说明 pivot 是最大值

		if      (j < hi)     return j;
		else if (hi-lo == 2) return hi-1;

		hi--;					// 去掉 hi 再划分一次（尾递归）
		}
	}

static nlkdnode* build_kd_tree
   (nlctx*	c,
	u32		lo,
	u32		hi,
	int		axis)
	{
	nlkdnode*	p;
	u32			m;

	p = &c->nodes[c->nodesLen++];
	p->maxN  = unreached;
	p->maxL  = unreached;
	p->loSon = p->hiSon = NULL;

	if (hi+1-lo <= bucketSize)
		{
		p->isBucket = true;
		p->loIx     = lo;
		p->hiIx     = hi;
		}
	else
		{
		p->isBucket = false;
		m = partition_points (c, lo, hi, axis);
		p->cutVal = projection (c, m, axis);
		p->loIx   = lo;
		p->hiIx   = m;
		p->loSon  = build_kd_tree (c, lo,  m,  1-axis);
		p->hiSon  = build_kd_tree (c, m+1, hi, 1-axis);
		}

	return p;
	}

// build_trees--
//	为 [a,b) 中的 HSP 按 (p,σ) 各建一棵树。

static const nlctx* groupSortCtx;

static int group_order
   (const void*	va,
	const void*	vb)
	{
	u32				a  = *((const u32*) va);
	u32				b  = *((const u32*) vb);
	const nlhsp*	ha = &groupSortCtx->hsp[a];
	const nlhsp*	hb = &groupSortCtx->hsp[b];

	if (ha->p     != hb->p)     return (ha->p     < hb->p)?     -1 : 1;
	if (ha->sigma != hb->sigma) return (ha->sigma < hb->sigma)? -1 : 1;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

static void build_trees
   (nlctx*	c,
	u32		a,
	u32		b)
	{
	u32			k, lo, hi;
	nlkdnode*	root;

	c->nodesLen = 0;
	if (b <= a) return;

	for (k=a ; k<b ; k++) c->perm[k] = k;
	groupSortCtx = c;
	qsort (&c->perm[a], b-a, sizeof(u32), group_order);

	for (lo=a ; lo<b ; lo=hi+1)
		{
		for (hi=lo ;
		     (hi+1 < b)
		  && (c->hsp[c->perm[hi+1]].p     == c->hsp[c->perm[lo]].p)
		  && (c->hsp[c->perm[hi+1]].sigma == c->hsp[c->perm[lo]].sigma) ;
		     hi++)
			;
		root = build_kd_tree (c, lo, hi, 0);
		for (k=lo ; k<=hi ; k++) c->treeOf[c->perm[k]] = root;
		}

	for (k=a ; k<b ; k++) c->treePos[c->perm[k]] = k;
	}

// kd_propagate--
//	i 的状态算完后，沿根到 i 所在叶子更新两个子树最大值（对应 chain.c 的
//	propagate_max_score）。

static void kd_propagate
   (nlctx*	c,
	u32		i)
	{
	const nlstate*	st  = &c->st[i];
	double			vN  = (st->U.score > st->R.score)? st->U.score : st->R.score;
	double			vL  = st->L.score;
	u32				pos = c->treePos[i];
	nlkdnode*		node;

	for (node=c->treeOf[i] ; node!=NULL ; )
		{
		if (vN > node->maxN) node->maxN = vN;
		if (vL > node->maxL) node->maxL = vL;
		if (node->isBucket) break;
		node = (pos <= node->hiIx)? node->loSon : node->hiSon;
		}
	}

// search_normal--
//	在一棵 (p,σ) 树里做 box 查询，对 box 内属于 Normal(i) 的 j 调用 E3、E6。

static void search_normal
   (nlctx*			c,
	const nlkdnode*	node,
	int				axis,
	u32				i,
	const nlbox*	box)
	{
	const nlstate*	st = &c->st[i];
	double			s  = c->hsp[i].s;
	u32				k, j;

	// 子树里还没有算完的 HSP

	if ((node->maxN == unreached) && (node->maxL == unreached)) return;

	// gain ≤ s_i、C_N ≥ 0（§6.1），所以子树内任一 j 的 E3 候选 ≤ maxN + s_i，
	// E6 候选 ≤ maxL + s_i；两者都严格低于当前最优时整棵跳过（用严格 <：等于
	// 时候选仍可能按平局规则胜出）

	if ((node->maxN + s < st->R.score) && (node->maxL + s < st->L.score)) return;

	if (node->isBucket)
		{
		for (k=node->loIx ; k<=node->hiIx ; k++)
			{
			j = c->perm[k];
			if (!is_normal_edge (&c->hsp[j], &c->hsp[i], c->params)) continue;
			relax_normal (c, j, i);
			relax_local  (c, j, i);
			}
		return;
		}

	// lo 子树坐标 ≤ cutVal，hi 子树坐标 ≥ cutVal

	if (box->lo[axis] <= (double) node->cutVal)
		search_normal (c, node->loSon, 1-axis, i, box);
	if (box->hi[axis] >= (double) node->cutVal)
		search_normal (c, node->hiSon, 1-axis, i, box);
	}

// search_jump--
//	从跳转池顶部按 R 降序扫描，对每个 j 调用 E2（§6.2）。

static void search_jump
   (nlctx*	c,
	u32		i)
	{
	const nlstate*	st = &c->st[i];
	const nlhsp*	hi = &c->hsp[i];
	u32				scanned = 0;
	u32				pos, j;

	if (c->pool.len == 0) return;

	c->frontier.len = 0;
	heap_push (c, &c->frontier, 0, frontier_above);
	while (c->frontier.len > 0)
		{
		pos = heap_pop (c, &c->frontier, frontier_above);
		j   = c->pool.v[pos];

		// 跳转边 ov = 0 ⇒ gain = s_i；P_J ≥ J ⇒ 候选 ≤ R[j] + s_i - J（§6.2）。
		// 按 R 降序遍历，上界严格低于目前最优时后面的 j 都不可能更好（用严格
		// <：等于时候选仍可能按平局规则胜出）

		if (c->st[j].R.score + hi->s - c->params->jumpPen < st->U.score) break;

		// 安全阀 H_max：只有这里会造成有损

		if (scanned >= c->params->maxJumpScan)
			{ c->jumpTruncations++;  break; }
		scanned++;

		if (is_jump_edge (&c->hsp[j], hi, c->params))
			relax_jump (c, j, i);

		if (2*pos+1 < c->pool.len) heap_push (c, &c->frontier, 2*pos+1, frontier_above);
		if (2*pos+2 < c->pool.len) heap_push (c, &c->frontier, 2*pos+2, frontier_above);
		}
	}

// search_candidates--
//	正式搜索层：K-d 树枚举 Normal(i)，跳转池枚举 Jump(i)。

static void search_candidates
   (nlctx*	c,
	u32		i)
	{
	const nlhsp*	hi       = &c->hsp[i];
	double			ov_max_i = c->params->overlapFrac * hsp_len (hi);
	double			sTin     = hi->sigma * (double) hi->tin;
	nlbox			box;

	// box 恰好对应 -ov_max_i ≤ dq, dt ≤ maxGap（§6.1）

	box.lo[0] = ((double) hi->qb) - c->params->maxGap;
	box.hi[0] = ((double) hi->qb) + ov_max_i;
	box.lo[1] = sTin - c->params->maxGap;
	box.hi[1] = sTin + ov_max_i;
	search_normal (c, c->treeOf[i], 0, i, &box);

	search_jump (c, i);
	}

// release_pending--
//	新的 qb 批次开始前，把 qe_j ≤ qb 的 j 从待释放堆移入跳转池（§6.2）；不满足
//	E5 的 j 不入池。

static void release_pending
   (nlctx*	c,
	unspos	qb)
	{
	u32		j;

	while ((c->pending.len > 0) && (c->hsp[c->pending.v[0]].qe <= qb))
		{
		j = heap_pop (c, &c->pending, pending_above);
		if (is_jump_source (c, j))
			heap_push (c, &c->pool, j, pool_above);
		}
	}

// run_dp--
//	对 [a,b) 中的 HSP 做 DP。同一 qb 批次的 HSP 读取批次开始前的池快照，整批
//	算完再进待释放堆（§6.2 同批冻结）。search 为搜索层，单元测试可换成暴力
//	枚举。

static void run_dp
   (nlctx*		c,
	u32			a,
	u32			b,
	nlsearch	search)
	{
	u32			k, e, i;
	unspos		qb;

	c->rangeLo      = a;
	c->rangeHi      = b;
	c->pending.len  = 0;
	c->pool.len     = 0;
	build_trees (c, a, b);

	for (k=a ; k<b ; k=e)
		{
		qb = c->hsp[k].qb;
		release_pending (c, qb);

		for (e=k ; (e<b)&&(c->hsp[e].qb==qb) ; e++)
			{
			init_state   (c, e);
			search       (c, e);
			settle_state (c, e);
			kd_propagate (c, e);
			}

		for (i=k ; i<e ; i++)
			{
			if (c->st[i].R.score != unreached)
				heap_push (c, &c->pending, i, pending_above);
			}
		}
	}

//----------
//
// 回溯（§7）
//
//----------

// emit_segment--
//	把一个共线段加入输出。hsps 按路径次序；段分数按 §7 第 3 步重新累加：段首
//	s 加段内每条正常边的 gain − C_N。first 暂存段在 members[] 中的起点，压缩
//	segtable 时再改成表内下标。

static void emit_segment
   (nlctx*		c,
	const u32*	hsps,
	u32			count,
	char		jumpIn,
	int			isSecondary,
	u32			chainId,
	u32			part2)
	{
	nlchainresult*	out = c->out;
	nlsegment*		seg;
	const nlhsp*	h, *prev;
	u32				k;

	if (out->len >= out->size)
		{
		out->size = (out->size == 0)? 64 : 2*out->size;
		out->seg  = (nlsegment*) realloc_or_die ("nlchain result", out->seg,
		                                         out->size * sizeof(nlsegment));
		}

	seg = &out->seg[out->len++];
	h   = &c->hsp[hsps[0]];
	seg->first       = c->membersLen;
	seg->count       = count;
	seg->part1       = h->p;
	seg->part2       = part2;
	seg->strand      = h->sigma;
	seg->qBeg        = h->qb;
	seg->qEnd        = h->qe;
	seg->tBeg        = hsp_tb (h);
	seg->tEnd        = hsp_te (h);
	seg->support     = h->s;
	seg->jumpIn      = jumpIn;
	seg->isSecondary = isSecondary;
	seg->chainId     = chainId;

	for (k=0 ; k<count ; k++)
		{
		h = &c->hsp[hsps[k]];
		c->members[c->membersLen++] = hsps[k];
		if (h->qe        > seg->qEnd) seg->qEnd = h->qe;
		if (hsp_tb (h) < seg->tBeg)   seg->tBeg = hsp_tb (h);
		if (hsp_te (h) > seg->tEnd)   seg->tEnd = hsp_te (h);
		if (k > 0)
			{
			prev = &c->hsp[hsps[k-1]];
			seg->support += edge_gain (prev, h) - normal_cost (prev, h, c->params);
			}
		}
	}

// trace_primary--
//	§7 第 1–3 步：选终点、回溯主路径、在跳转边处切段。返回主路径分。

static double trace_primary
   (nlctx*	c,
	u32		a,
	u32		b,
	u32		part2)
	{
	const nlpath*	path;
	u32				i, e, k, k0, t;
	u8				state, ts;

	c->pathLen = 0;

	// 1. 终点：所有 R[i] > 0 中按平局规则取最大者（再平则下标小者）

	e = noPred;
	for (i=a ; i<b ; i++)
		{
		if (c->st[i].R.score <= 0) continue;
		if ((e == noPred) || (better (&c->st[i].R, &c->st[e].R))) e = i;
		}
	if (e == noPred) return 0;

	// 2. 从 (e, R) 沿 pred 回溯，再反转成 query 升序

	i = e;  state = stateR;
	while (i != noPred)
		{
		c->pathIx   [c->pathLen] = i;
		c->pathState[c->pathLen] = state;
		c->pathLen++;
		path  = (state == stateU)? &c->st[i].U : &c->st[i].R;
		i     = path->pred;
		state = path->predState;
		}

	for (k=0 ; k<c->pathLen/2 ; k++)
		{
		t  = c->pathIx   [k];  c->pathIx   [k] = c->pathIx   [c->pathLen-1-k];  c->pathIx   [c->pathLen-1-k] = t;
		ts = c->pathState[k];  c->pathState[k] = c->pathState[c->pathLen-1-k];  c->pathState[c->pathLen-1-k] = ts;
		}

	for (k=0 ; k<c->pathLen ; k++) c->visited[c->pathIx[k]] = true;

	// 3. 在 edge = J 处切段

	for (k0=0,k=1 ; k<=c->pathLen ; k++)
		{
		if (k < c->pathLen)
			{
			path = (c->pathState[k] == stateU)? &c->st[c->pathIx[k]].U
			                                  : &c->st[c->pathIx[k]].R;
			if (path->edge != edgeJump) continue;
			}

		emit_segment (c, &c->pathIx[k0], k-k0,
		              (k0 == 0)? 'S' : jump_kind (&c->hsp[c->pathIx[k0-1]],
		                                          &c->hsp[c->pathIx[k0]]),
		              /* isSecondary */ false, /* chainId */ 0, part2);
		k0 = k;
		}

	return c->st[e].R.score;
	}

// primary_score_within--
//	主路径在 query 区间 [qa,qz) 内的分数（§5.4 第 6 步）：取 qb 落在区间内的
//	主路径 HSP，用回溯时的累计分之差求出。区间内没有主路径时返回 false。

static int primary_score_within
   (const nlctx*	c,
	unspos			qa,
	unspos			qz,
	double*			score)
	{
	u32		lo, hi, mid, first, last;
	double	before, through;

	// first = 第一个 qb ≥ qa 的路径下标

	lo = 0;  hi = c->pathLen;
	while (lo < hi)
		{
		mid = (lo+hi) / 2;
		if (c->hsp[c->pathIx[mid]].qb < qa) lo = mid+1; else hi = mid;
		}
	first = lo;

	// last+1 = 第一个 qb ≥ qz 的路径下标

	lo = first;  hi = c->pathLen;
	while (lo < hi)
		{
		mid = (lo+hi) / 2;
		if (c->hsp[c->pathIx[mid]].qb < qz) lo = mid+1; else hi = mid;
		}
	if (lo == first) return false;
	last = lo - 1;

	through = (c->pathState[last] == stateU)? c->st[c->pathIx[last]].U.score
	                                        : c->st[c->pathIx[last]].R.score;
	before  = 0;
	if (first > 0)
		before = (c->pathState[first-1] == stateU)? c->st[c->pathIx[first-1]].U.score
		                                          : c->st[c->pathIx[first-1]].R.score;
	*score = through - before;
	return true;
	}

// trace_secondary--
//	§5.4 第 2–6 步：从未访问、L ≥ T_sec 的 HSP 按 L 降序沿 predL 回溯次级链。

static const nlctx* candSortCtx;

static int cand_order
   (const void*	va,
	const void*	vb)
	{
	u32				a  = *((const u32*) va);
	u32				b  = *((const u32*) vb);
	const nlpath*	la = &candSortCtx->st[a].L;
	const nlpath*	lb = &candSortCtx->st[b].L;

	if (better (la, lb)) return -1;
	if (better (lb, la)) return  1;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

static void trace_secondary
   (nlctx*	c,
	u32		a,
	u32		b,
	u32		part2)
	{
	const nlchainparams* p = c->params;
	u32			numCands, ix, i, k, len, t, numOverlapping, firstSecondary, chainId;
	double		stopScore, secScore, primScore;
	unspos		qa, qz, ka, kz, ovBeg, ovEnd;
	nlsegment*	kept;

	if (p->secondaryScore <= 0) return;		// T_sec = 0 关闭次级链

	// 2. 候选终点

	numCands = 0;
	for (i=a ; i<b ; i++)
		{
		if ((!c->visited[i]) && (c->st[i].L.score >= p->secondaryScore))
			c->cands[numCands++] = i;
		}
	candSortCtx = c;
	qsort (c->cands, numCands, sizeof(u32), cand_order);

	firstSecondary = c->out->len;
	chainId        = 1;

	for (ix=0 ; ix<numCands ; ix++)
		{
		i = c->cands[ix];
		if (c->visited[i]) continue;

		// 3. 沿 predL 回溯，遇到已访问的 HSP 就停

		len = 0;
		stopScore = 0;
		while (true)
			{
			c->chainBuf[len++] = i;
			if (c->st[i].L.pred == noPred) break;
			if (c->visited[c->st[i].L.pred])
				{ stopScore = c->st[c->st[i].L.pred].L.score;  break; }
			i = c->st[i].L.pred;
			}

		// 4. 次级链分数 = L[end] − L[stop]

		secScore = c->st[c->cands[ix]].L.score - stopScore;
		if (secScore < p->secondaryScore) continue;

		for (k=0 ; k<len/2 ; k++)
			{ t = c->chainBuf[k];  c->chainBuf[k] = c->chainBuf[len-1-k];  c->chainBuf[len-1-k] = t; }

		qa = c->hsp[c->chainBuf[0]].qb;
		qz = qa;
		for (k=0 ; k<len ; k++)
			{ if (c->hsp[c->chainBuf[k]].qe > qz) qz = c->hsp[c->chainBuf[k]].qe; }

		// 5. 每段 query 区间限额：已保留的次级链中，与它的 query 重叠超过它自身
		//    query 长度一半的，达到 N_sec 条就丢弃

		numOverlapping = 0;
		for (k=firstSecondary ; k<c->out->len ; k++)
			{
			kept  = &c->out->seg[k];
			ka    = kept->qBeg;  kz = kept->qEnd;
			ovBeg = (ka > qa)? ka : qa;
			ovEnd = (kz < qz)? kz : qz;
			if ((ovEnd > ovBeg) && ((double) (ovEnd - ovBeg) > 0.5 * (qz - qa)))
				numOverlapping++;
			}
		if (numOverlapping >= p->secondaryPerRegion) continue;

		// 6. 相对同一 query 区间内主路径分数的比例

		if ((primary_score_within (c, qa, qz, &primScore))
		 && (secScore < p->secondaryRatio * primScore))
			continue;

		// 保留

		for (k=0 ; k<len ; k++) c->visited[c->chainBuf[k]] = true;
		emit_segment (c, c->chainBuf, len, 'S',
		              /* isSecondary */ true, chainId++, part2);
		}
	}

//----------
//
// 对外接口
//
//----------

// nlsortrec--
//	换算后的 HSP 及其来源，用于排序。

typedef struct nlsortrec
	{
	nlhsp		h;
	u32			part2;
	u8			table;
	u32			ix;
	} nlsortrec;

// DP 次序：query 分区，再按 (qb, p, tb, σ, len)（§4），最后按来源保证确定

static int sortrec_order
   (const void*	va,
	const void*	vb)
	{
	const nlsortrec* a = (const nlsortrec*) va;
	const nlsortrec* b = (const nlsortrec*) vb;

	if (a->part2      != b->part2)      return (a->part2      < b->part2)?      -1 : 1;
	if (a->h.qb       != b->h.qb)       return (a->h.qb       < b->h.qb)?       -1 : 1;
	if (a->h.p        != b->h.p)        return (a->h.p        < b->h.p)?        -1 : 1;
	if (hsp_tb(&a->h) != hsp_tb(&b->h)) return (hsp_tb(&a->h) < hsp_tb(&b->h))? -1 : 1;
	if (a->h.sigma    != b->h.sigma)    return (a->h.sigma    < b->h.sigma)?    -1 : 1;
	if (hsp_len(&a->h)!= hsp_len(&b->h))return (hsp_len(&a->h)< hsp_len(&b->h))?-1 : 1;
	if (a->table      != b->table)      return (a->table      < b->table)?      -1 : 1;
	return (a->ix < b->ix)? -1 : (a->ix > b->ix)? 1 : 0;
	}

// convert_hsp--
//	§4 的坐标换算；正反链的分支只出现在这里。反链 pos2 是反向互补 query 上的
//	坐标；query 有分区时每个分区各自就地反向互补，所以按分区镜像。

static void convert_hsp
   (seq*			target,
	seq*			query,
	const segment*	seg,
	int				sigma,
	nlhsp*			h,
	u32*			part2)
	{
	seqpartition*	sp1 = &target->partition;
	seqpartition*	sp2 = &query->partition;
	partition*		part;
	unspos			qFrom, qTo, tb, te;

	if (sp2->p == NULL)
		{ *part2 = 0;  qFrom = 0;  qTo = query->len; }
	else
		{
		part   = lookup_partition (query, seg->pos2);
		*part2 = (u32) (part - sp2->p);
		qFrom  = part->sepBefore + 1;
		qTo    = part->sepAfter;
		}

	if (sigma > 0) h->qb = seg->pos2;
	          else h->qb = qFrom + qTo - seg->pos2 - seg->length;
	h->qe = h->qb + seg->length;

	tb = seg->pos1;
	te = seg->pos1 + seg->length;
	h->tin   = (sigma > 0)? tb : te;
	h->tout  = (sigma > 0)? te : tb;
	h->sigma = sigma;
	h->p     = (sp1->p == NULL)? 0 : (u32) (lookup_partition (target, seg->pos1) - sp1->p);
	h->s     = (double) seg->s;
	h->rho   = h->s / seg->length;
	}

// compact_tables--
//	§7 第 5 步：两张表各自只留下主路径和次级链的 HSP，并按输出段的次序排列，
//	使每段在表内连续；同时把段的 first 改为表内下标。

static void compact_tables
   (nlctx*		c,
	segtable*	stMinus,
	segtable*	stPlus)
	{
	segtable*	tables[2] = { stMinus, stPlus };
	segment*	copies[2];
	u32			cur[2] = { 0, 0 };
	u32			s, k, t, m;
	nlsegment*	seg;

	for (t=0 ; t<2 ; t++)
		{
		copies[t] = (segment*) malloc_or_die ("nlchain compact",
		                                      (tables[t]->len+1) * sizeof(segment));
		memcpy (copies[t], tables[t]->seg, tables[t]->len * sizeof(segment));
		}

	for (s=0 ; s<c->out->len ; s++)
		{
		seg = &c->out->seg[s];
		t   = (seg->strand < 0)? 0 : 1;
		m   = seg->first;
		seg->first = cur[t];
		for (k=0 ; k<seg->count ; k++)
			{
			tables[t]->seg[cur[t]] = copies[t][c->srcIx[c->members[m+k]]];
			tables[t]->seg[cur[t]].filter = false;
			cur[t]++;
			}
		}

	for (t=0 ; t<2 ; t++)
		{
		tables[t]->len = cur[t];
		free_if_valid ("nlchain compact", copies[t]);
		}
	}

//----------
//
// nl_reduce_to_chains--
//	对一条 query 两条链的 HSP 做非线性串联，只保留主路径和次级链的 HSP。
//
//----------
//
// Arguments:
//	seq*		target:		target 序列（可带 [multiple] 分区）。
//	seq*		query:		query 序列；调用时处于反向互补状态也可以（只用坐标）。
//	segtable*	stMinus:	反链 HSP（单链表，已 merge、已带分数）。
//	segtable*	stPlus:		正链 HSP（同上）。
//	nlchainparams* params:	参数，hspthresh 相关的默认值已由调用者解析好。
//	nlchainresult* out:		输出各段，按 (query 分区, chainId, qBeg) 排列。
//
// Returns:
//	主路径分（各 query 分区中的最大者）；两张表原地压缩。
//
//----------

double nl_reduce_to_chains
   (seq*					target,
	seq*					query,
	segtable*				stMinus,
	segtable*				stPlus,
	const nlchainparams*	params,
	nlchainresult*			out)
	{
	nlctx			c;
	nlsortrec*		recs;
	segtable*		tables[2] = { stMinus, stPlus };
	u32				n, k, t, ix, a, b;
	double			best, primScore;
#ifdef collect_stats
	u32				s;
#endif // collect_stats

	out->len = 0;
	n = stMinus->len + stPlus->len;
	if (n == 0) return 0;

	ctx_alloc (&c, n);
	c.params = params;
	c.out    = out;

	// 统一坐标（§4）并排序

	recs = (nlsortrec*) malloc_or_die ("nlchain recs", n * sizeof(nlsortrec));
	k = 0;
	for (t=0 ; t<2 ; t++)
		{
		for (ix=0 ; ix<tables[t]->len ; ix++,k++)
			{
			convert_hsp (target, query, &tables[t]->seg[ix], (t == 0)? -1 : +1,
			             &recs[k].h, &recs[k].part2);
			recs[k].table = (u8) t;
			recs[k].ix    = ix;
			}
		}
	qsort (recs, n, sizeof(nlsortrec), sortrec_order);
	for (k=0 ; k<n ; k++)
		{
		c.hsp[k]      = recs[k].h;
		c.part2[k]    = recs[k].part2;
		c.srcTable[k] = recs[k].table;
		c.srcIx[k]    = recs[k].ix;
		}
	free_if_valid ("nlchain recs", recs);

	// 每个 query 分区分别做 DP 与回溯

	best = 0;
	for (a=0 ; a<n ; a=b)
		{
		for (b=a ; (b<n)&&(c.part2[b]==c.part2[a]) ; b++)
			;
		run_dp (&c, a, b, search_candidates);
#ifdef collect_stats
		s = out->len;
#endif // collect_stats
		primScore = trace_primary (&c, a, b, c.part2[a]);
		if (primScore > best) best = primScore;

#ifdef collect_stats
		for ( ; s<out->len ; s++)
			{
			nlchain_count_stat (numPrimarySegments);
			nlchain_add_stat   (numPrimaryHsps, out->seg[s].count);
			switch (out->seg[s].jumpIn)
				{
				case 'S': nlchain_count_stat (numJumpIn[0]);  break;
				case 'I': nlchain_count_stat (numJumpIn[1]);  break;
				case 'T': nlchain_count_stat (numJumpIn[2]);  break;
				case 'D': nlchain_count_stat (numJumpIn[3]);  break;
				case 'G': nlchain_count_stat (numJumpIn[4]);  break;
				}
			}
#endif // collect_stats

		trace_secondary (&c, a, b, c.part2[a]);

#ifdef collect_stats
		for ( ; s<out->len ; s++)
			{
			nlchain_count_stat (numSecondaryChains);
			nlchain_add_stat   (numSecondaryHsps, out->seg[s].count);
			}
#endif // collect_stats
		}

	nlchain_add_stat (numHsps,            n);
	nlchain_add_stat (numJumpTruncations, c.jumpTruncations);

	compact_tables (&c, stMinus, stPlus);
	ctx_free (&c);

	return best;
	}

//----------
//
// nl_free_result--
//	释放 nlchainresult 的段数组。
//
//----------

void nl_free_result
   (nlchainresult*	out)
	{
	if (out->seg != NULL) free_if_valid ("nlchain result", out->seg);
	out->seg  = NULL;
	out->len  = 0;
	out->size = 0;
	}

//----------
//
// nl_write_report_header, nl_write_report--
//	段表 TSV（§7 第 7 步）。坐标为各自序列内从 0 起的半开区间；query 坐标为
//	正链坐标。
//
//----------

void nl_write_report_header
   (FILE*	f)
	{
	fprintf (f, "#query\tkind\tchain\tseg\tjumpIn\ttarget\tstrand"
	            "\tqBeg\tqEnd\ttBeg\ttEnd\tscore\thsps\n");
	}

static void seq_part_info
   (seq*		s,
	u32			partIx,
	const char*	defaultName,
	const char** name,
	unspos*		offset,
	unspos*		startLoc)
	{
	seqpartition*	sp = &s->partition;
	partition*		part;

	if (sp->p == NULL)
		{
		*name = (s->useFullNames)? s->header : s->shortHeader;
		if ((*name == NULL) || ((*name)[0] == 0)) *name = defaultName;
		*offset   = 0;
		*startLoc = s->startLoc;
		}
	else
		{
		part      = &sp->p[partIx];
		*name     = &sp->pool[part->header];
		*offset   = part->sepBefore + 1;
		*startLoc = part->startLoc;
		}
	}

void nl_write_report
   (FILE*					f,
	seq*					target,
	seq*					query,
	const nlchainresult*	r)
	{
	const nlsegment*	seg;
	const char*			name1, *name2;
	unspos				offset1, offset2, startLoc1, startLoc2;
	u32					s, segNum;

	if (f == NULL) return;

	segNum = 0;
	for (s=0 ; s<r->len ; s++)
		{
		seg = &r->seg[s];
		if ((s > 0)
		 && ((seg->part2 != r->seg[s-1].part2) || (seg->chainId != r->seg[s-1].chainId)))
			segNum = 0;

		seq_part_info (target, seg->part1, "seq1", &name1, &offset1, &startLoc1);
		seq_part_info (query,  seg->part2, "seq2", &name2, &offset2, &startLoc2);

		fprintf (f, "%s\t%s\t%u\t%u\t%c\t%s\t%c"
		            "\t" unsposFmt "\t" unsposFmt "\t" unsposFmt "\t" unsposFmt
		            "\t%.1f\t%u\n",
		         name2,
		         (seg->isSecondary)? "secondary" : "primary",
		         seg->chainId, segNum++, seg->jumpIn,
		         name1, (seg->strand > 0)? '+' : '-',
		         seg->qBeg - offset2 + startLoc2 - 1,
		         seg->qEnd - offset2 + startLoc2 - 1,
		         seg->tBeg - offset1 + startLoc1 - 1,
		         seg->tEnd - offset1 + startLoc1 - 1,
		         seg->support, seg->count);
		}
	}

//----------
//
// nlchain_zero_stats, nlchain_show_stats, nlchain_generic_stats--
//	统计（§8.2 第 7 条）。
//
//----------

void nlchain_zero_stats
   (void)
	{
#ifdef collect_stats
	memset (&nlchainStats, 0, sizeof(nlchainStats));
#endif // collect_stats
	}

void nlchain_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	if (nlchainStats.numHsps == 0) return;
	fprintf (f, "      nlchain HSPs: %s\n", ucommatize(nlchainStats.numHsps));
	fprintf (f, "      primary HSPs: %s\n", ucommatize(nlchainStats.numPrimaryHsps));
	fprintf (f, "  primary segments: %s\n", ucommatize(nlchainStats.numPrimarySegments));
	fprintf (f, "  jumpIn S/I/T/D/G: %s/%s/%s/%s/%s\n",
	            ucommatize(nlchainStats.numJumpIn[0]), ucommatize(nlchainStats.numJumpIn[1]),
	            ucommatize(nlchainStats.numJumpIn[2]), ucommatize(nlchainStats.numJumpIn[3]),
	            ucommatize(nlchainStats.numJumpIn[4]));
	fprintf (f, "  secondary chains: %s\n", ucommatize(nlchainStats.numSecondaryChains));
	fprintf (f, "    secondary HSPs: %s\n", ucommatize(nlchainStats.numSecondaryHsps));
	fprintf (f, "  jump truncations: %s\n", ucommatize(nlchainStats.numJumpTruncations));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void nlchain_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "nlchain_hsps=%" PRIu64 "\n",              nlchainStats.numHsps);
	(*func) (f, "nlchain_primary_hsps=%" PRIu64 "\n",      nlchainStats.numPrimaryHsps);
	(*func) (f, "nlchain_primary_segments=%" PRIu64 "\n",  nlchainStats.numPrimarySegments);
	(*func) (f, "nlchain_secondary_chains=%" PRIu64 "\n",  nlchainStats.numSecondaryChains);
	(*func) (f, "nlchain_secondary_hsps=%" PRIu64 "\n",    nlchainStats.numSecondaryHsps);
	(*func) (f, "nlchain_jump_truncations=%" PRIu64 "\n",  nlchainStats.numJumpTruncations);
#endif // collect_stats
	}
