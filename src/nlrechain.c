//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlrechain.c
//
//----------
//
// nlrechain--
//	nlchain 的 gapped 后处理（方案 A）：把一条链方向上的全部 gapped 比对拆成无
//	gap 块，按 axtChain 的模型重新串链，删掉重算后链分低于 T_chain 的链，再把保
//	留的块放回各自的来源比对、在被删或被裁处切开输出。
//
// 算法与记号见 docs/nlchain_rechain_plan.md，下文的 §x 均指该文档章节。行为以
// kent axtChain 为准（按 §2 的描述独立实现）。代码分层：
//	公式层：gap 代价表、gapCalcCost、交叉点、connectCost、块分，纯算术。
//	转移层：offer() 是唯一的状态更新点，平局只在这里判断。
//	搜索层：2-d 树枚举候选前驱，剪枝条件只出现在这一层。
//	其余：拆块、合并完全重复块、剥链、链内整理、链分、输出重建、报告。
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
#include <float.h>				// standard C floating point limits
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "gapped_extend.h"		// score_alignment

#define  nlrechain_owner		// (make this the owner of its globals)
#include "nlrechain.h"			// interface to this module

//----------
//
// 内置代价表（§2 第 4 步）
//
// 数值为 kent gapCalc 的 medium（originalGapCosts）和 loose（defaultGapCosts）
// 两张表；qGap 与 tGap 两行相同。
//
//----------

static const unspos nlrGapPositions[11]
	= { 1, 2, 3, 11, 111, 2111, 12111, 32111, 72111, 152111, 252111 };
static const s64 nlrMediumGap[11]
	= { 350, 425, 450, 600, 900, 2900, 22900, 57900, 117900, 217900, 317900 };
static const s64 nlrMediumBoth[11]
	= { 750, 825, 850, 1000, 1300, 3300, 23300, 58300, 118300, 218300, 318300 };
static const s64 nlrLooseGap[11]
	= { 325, 360, 400, 450, 600, 1100, 3600, 7600, 15600, 31600, 56600 };
static const s64 nlrLooseBoth[11]
	= { 625, 660, 700, 750, 900, 1400, 4000, 8000, 16000, 32000, 57000 };

const nlrgaptable nlrGapCostMedium
	= { 11, nlrGapPositions, nlrMediumGap, nlrMediumGap, nlrMediumBoth };
const nlrgaptable nlrGapCostLoose
	= { 11, nlrGapPositions, nlrLooseGap,  nlrLooseGap,  nlrLooseBoth  };

//----------
//
// 内部数据结构
//
//----------

#define noIx ((u32) -1)

// 两侧完全被包含时的连接代价（axtChain 的"禁止"值）
#define enclosedCost 100000000.0

// nlrspan--
//	一个无 gap 块的范围，target/query 都是 seq->v 上从 0 起的半开区间（反链时
//	query 为反向互补后的坐标）。

typedef struct nlrspan
	{
	unspos		tStart, tEnd;
	unspos		qStart, qEnd;
	} nlrspan;

// nlrblock--
//	拆出的块（§2 第 1 步）及其 DP 状态（§2 第 4 步）。

typedef struct nlrblock
	{
	nlrspan		span;			// 块的范围
	u32			part1, part2;	// target / query 分区下标
	u32			src;			// 来源比对在输入列表中的下标
	u32			ord;			// 在来源比对中的次序
	u32			opBegin;		// 来源 script 中本块第一个 sub 操作的下标
	u32			opAfter;		// 本块之后第一个非 sub 操作的下标
	double		score;			// 块分：替换分之和（§2 第 2 步）
	double		minCol;			// 块内最低的单列替换分（搜索层的下界用）
	double		total;			// 以本块结尾的最优链分
	u32			pred;			// 最优前驱；noIx 表示链首
	} nlrblock;

// nlrpiece--
//	链中的一块，范围可能已被链内整理裁短（§2 第 6 步）。

typedef struct nlrpiece
	{
	nlrspan		span;
	u32			blk;			// 对应的 nlrblock 下标
	} nlrpiece;

// nlrbest--
//	某块的最优前驱候选（转移层）。

typedef struct nlrbest
	{
	double		score;
	u32			pred;
	} nlrbest;

// nlrnode--
//	2-d 树结点；每个结点存一个块，聚合值覆盖整棵子树（搜索层）。

typedef struct nlrnode
	{
	u32			blk;			// 本结点的块
	u32			lo, hi;			// 子结点；noIx 表示无
	u32			parent;			// 父结点；noIx 表示根
	unspos		minQ, minT;		// 子树内 qStart、tStart 的最小值
	unspos		maxQE, maxTE;	// 子树内 qEnd、tEnd 的最大值
	double		minCol;			// 子树内块的 minCol 最小值
	double		maxTotal;		// 子树内已处理块的 total 最大值；无则 -DBL_MAX
	} nlrnode;

// nlrctx--
//	一次 nl_rechain_aligns 调用的上下文。

typedef struct nlrctx
	{
	seq*		target;
	seq*		query;
	scoreset*	scoring;
	const nlrgaptable* gapCost;
	} nlrctx;

// nlrrun--
//	输出的一段：同一来源比对中次序连续、在同一条链中相邻的若干块（X2）。

typedef struct nlrrun
	{
	u32			first, count;	// 在 kept[] 中的下标区间
	u32			src, ord;		// 来源比对、首块次序（输出排序用，X3）
	} nlrrun;

//----------
//
// 公式层
//
//----------

// 代价表一行：d = 0 时为 0；表内线性插值后向下取整；超出表尾按最后一段的斜率
// 外推   (§2 第 4 步 gapCalcCost)

static s64 gap_row_cost
   (const nlrgaptable*	t,
	const s64*			cost,
	u64					d)
	{
	u32		k, last;
	u64		ds;
	s64		dv;

	if (d == 0) return 0;
	if (d <= t->position[0]) return cost[0];	// （表首为 1）

	last = t->len - 1;
	if (d >= t->position[last])
		k = last;
	else
		for (k=1 ; t->position[k]<d ; k++)		// position[k-1] < d ≤ position[k]
			;

	ds = t->position[k] - t->position[k-1];
	dv = cost[k] - cost[k-1];
	if (d >= t->position[last])
		return cost[last] + (s64) ((dv * (s64) (d - t->position[last])) / (s64) ds);
	return cost[k-1] + (s64) ((dv * (s64) (d - t->position[k-1])) / (s64) ds);
	}

// gapCalcCost(dq,dt)：负值归零；dt = 0 查 qGap(dq)，dq = 0 查 tGap(dt)，
// 否则查 bothGap(dq+dt)   (§2 第 4 步)

static s64 gap_calc_cost
   (const nlrgaptable*	t,
	s64					dq,
	s64					dt)
	{
	if (dq < 0) dq = 0;
	if (dt < 0) dt = 0;
	if (dt == 0) return gap_row_cost (t, t->qGap, (u64) dq);
	if (dq == 0) return gap_row_cost (t, t->tGap, (u64) dt);
	return gap_row_cost (t, t->bothGap, (u64) (dq + dt));
	}

// 单列替换分 sub[target][query]（与 score_alignment 相同）

static double column_score
   (const nlrctx*	c,
	unspos			tPos,
	unspos			qPos)
	{
	return (double) c->scoring->sub[c->target->v[tPos]][c->query->v[qPos]];
	}

// 块分 = 块内替换分之和   (§2 第 2 步)

static double span_score
   (const nlrctx*	c,
	const nlrspan*	b)
	{
	unspos	k, len = b->tEnd - b->tStart;
	double	s = 0;

	for (k=0 ; k<len ; k++)
		s += column_score (c, b->tStart+k, b->qStart+k);
	return s;
	}

// 交叉点：在重叠的 ov 列上，先全用右块，再从左往右逐列改用左块，取总分最高的
// 位置 pos（平局取靠前）；返回调整量 rScore + lScore - bestScore   (§2 第 4 步)

static double find_crossover
   (const nlrctx*	c,
	const nlrspan*	a,				// 左块
	const nlrspan*	b,				// 右块
	unspos			ov,
	unspos*			pos)
	{
	double	l, r, rScore = 0, lScore = 0, sc, bestScore;
	unspos	k, bestPos = 0;

	for (k=0 ; k<ov ; k++)
		{
		rScore += column_score (c, b->tStart+k,      b->qStart+k);
		lScore += column_score (c, a->tEnd-ov+k,     a->qEnd-ov+k);
		}

	sc = bestScore = rScore;
	for (k=0 ; k<ov ; k++)
		{
		l = column_score (c, a->tEnd-ov+k, a->qEnd-ov+k);
		r = column_score (c, b->tStart+k,  b->qStart+k);
		sc += l - r;
		if (sc > bestScore) { bestScore = sc;  bestPos = k+1; }
		}

	*pos = bestPos;
	return rScore + lScore - bestScore;
	}

// connectCost(a,b)：不重叠时为 gapCalcCost(dq,dt)；重叠 ov 列时，任一块被
// 完全包含则加 enclosedCost，否则加交叉点调整并令 dq、dt 各加 ov   (§2 第 4 步)

static double connect_cost
   (const nlrctx*	c,
	const nlrspan*	a,
	const nlrspan*	b)
	{
	s64		dq = ((s64) b->qStart) - ((s64) a->qEnd);
	s64		dt = ((s64) b->tStart) - ((s64) a->tEnd);
	s64		ov;
	double	adjust = 0;
	unspos	pos;

	if ((dq < 0) || (dt < 0))
		{
		ov = (dq < dt)? -dq : -dt;
		if ((ov >= (s64) (b->qEnd - b->qStart)) || (ov >= (s64) (a->qEnd - a->qStart)))
			adjust = enclosedCost;
		else
			{
			adjust = find_crossover (c, a, b, (unspos) ov, &pos);
			dq += ov;
			dt += ov;
			}
		}

	return adjust + (double) gap_calc_cost (c->gapCost, dq, dt);
	}

//----------
//
// 转移层
//
//----------

// 候选：total_j + s_i - connectCost(j,i)   (§2 第 4 步)

static double chain_candidate
   (const nlrctx*		c,
	const nlrblock*		blk,
	u32					j,
	u32					i)
	{
	return blk[j].total + blk[i].score - connect_cost (c, &blk[j].span, &blk[i].span);
	}

// is_predecessor--
//	合法前驱：qStart_j < qStart_i 且 tStart_j < tStart_i（都严格）。

static int is_predecessor
   (const nlrblock*	j,
	const nlrblock*	i)
	{
	return (j->span.qStart < i->span.qStart) && (j->span.tStart < i->span.tStart);
	}

// offer--
//	唯一的状态更新点。初值 {0, noIx}：只接受严格大于 0 的候选（axtChain 的搜索
//	初值为 0）；同分取下标较小的 j，使结果与枚举次序无关。

static void offer
   (nlrbest*	best,
	double		cand,
	u32			j)
	{
	if ((cand > best->score)
	 || ((cand == best->score) && (best->pred != noIx) && (j < best->pred)))
		{ best->score = cand;  best->pred = j; }
	}

// settle--
//	所有候选处理完后：最优候选高于块分本身时才接上前驱，否则自成链首。

static void settle
   (nlrblock*		b,
	const nlrbest*	best)
	{
	b->total = b->score;
	b->pred  = noIx;
	if ((best->pred != noIx) && (best->score > b->score))
		{ b->total = best->score;  b->pred = best->pred; }
	}

//----------
//
// 搜索层：2-d 树
//
// 块按 DP 次序（tStart、qStart 升序）编号，依次处理；处理完的块"激活"，其 total
// 沿父链更新各结点的 maxTotal。对块 i 的查询剪去三类子树：没有 qStart、tStart
// 都更小的块；没有已处理块；上界 maxTotal + s_i - 代价下界 严格小于当前最优
// （等于时不剪，平局交给 offer()）。
//
//----------

static const nlrblock* sortBlocks;		// qsort 比较函数用

static int cmp_by_q
   (const void*	va,
	const void*	vb)
	{
	u32 a = *(const u32*) va, b = *(const u32*) vb;
	if (sortBlocks[a].span.qStart != sortBlocks[b].span.qStart)
		return (sortBlocks[a].span.qStart < sortBlocks[b].span.qStart)? -1 : 1;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

static int cmp_by_t
   (const void*	va,
	const void*	vb)
	{
	u32 a = *(const u32*) va, b = *(const u32*) vb;
	if (sortBlocks[a].span.tStart != sortBlocks[b].span.tStart)
		return (sortBlocks[a].span.tStart < sortBlocks[b].span.tStart)? -1 : 1;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

// build_tree--
//	在 perm[lo,hi) 上按 axis（0 = qStart，1 = tStart）取中位数递归建树，返回
//	结点下标。

static u32 build_tree
   (const nlrblock*	blk,
	u32*			perm,
	u32				lo,
	u32				hi,
	int				axis,
	u32				parent,
	nlrnode*		node,
	u32*			numNodes,
	u32*			nodeOf)
	{
	u32			m, id, ch[2], k;
	nlrnode*	n;
	nlrnode*	child;

	if (lo >= hi) return noIx;

	sortBlocks = blk;
	qsort (perm+lo, hi-lo, sizeof(u32), (axis == 0)? cmp_by_q : cmp_by_t);
	m  = lo + (hi-lo)/2;
	id = (*numNodes)++;

	n = &node[id];
	n->blk      = perm[m];
	n->parent   = parent;
	n->minQ     = blk[perm[m]].span.qStart;
	n->minT     = blk[perm[m]].span.tStart;
	n->maxQE    = blk[perm[m]].span.qEnd;
	n->maxTE    = blk[perm[m]].span.tEnd;
	n->minCol   = blk[perm[m]].minCol;
	n->maxTotal = -DBL_MAX;
	nodeOf[perm[m]] = id;

	ch[0] = build_tree (blk, perm, lo,  m,  1-axis, id, node, numNodes, nodeOf);
	ch[1] = build_tree (blk, perm, m+1, hi, 1-axis, id, node, numNodes, nodeOf);

	n = &node[id];						// （node 不会重新分配，这里只为清楚）
	n->lo = ch[0];
	n->hi = ch[1];
	for (k=0 ; k<2 ; k++)
		{
		if (ch[k] == noIx) continue;
		child = &node[ch[k]];
		if (child->minQ   < n->minQ)   n->minQ   = child->minQ;
		if (child->minT   < n->minT)   n->minT   = child->minT;
		if (child->maxQE  > n->maxQE)  n->maxQE  = child->maxQE;
		if (child->maxTE  > n->maxTE)  n->maxTE  = child->maxTE;
		if (child->minCol < n->minCol) n->minCol = child->minCol;
		}
	return id;
	}

// activate--
//	块 i 处理完毕，把 total_i 并入沿途各结点的 maxTotal。

static void activate
   (nlrnode*		node,
	const u32*		nodeOf,
	const nlrblock*	blk,
	u32				i)
	{
	u32 id;

	for (id=nodeOf[i] ; id!=noIx ; id=node[id].parent)
		if (blk[i].total > node[id].maxTotal) node[id].maxTotal = blk[i].total;
	}

// cost_lower_bound--
//	子树内任一块 j 到块 i 的 connectCost 下界。dq_j ≥ qStart_i - maxQE、
//	dt_j ≥ tStart_i - maxTE，gapCalcCost 对两个参数都单调不减，所以 gap 部分
//	≥ gapCalcCost(qStart_i - maxQE, tStart_i - maxTE)。交叉点调整等于被舍弃
//	各列的替换分之和，可能为负：列数 ov 不超过 max(maxQE - qStart_i,
//	maxTE - tStart_i)，且不包含时 ov < |i|，每列 ≥ min(minCol_i, 子树 minCol)。

static double cost_lower_bound
   (const nlrctx*	c,
	const nlrnode*	n,
	const nlrblock*	bi)
	{
	s64		dq = ((s64) bi->span.qStart) - ((s64) n->maxQE);
	s64		dt = ((s64) bi->span.tStart) - ((s64) n->maxTE);
	s64		ovMax = (dq < dt)? -dq : -dt;
	s64		lenI  = (s64) (bi->span.qEnd - bi->span.qStart);
	double	minCol, lb;

	lb = (double) gap_calc_cost (c->gapCost, dq, dt);
	if (ovMax > 0)
		{
		if (ovMax > lenI - 1) ovMax = lenI - 1;
		minCol = (n->minCol < bi->minCol)? n->minCol : bi->minCol;
		if (minCol < 0) lb += ((double) ovMax) * minCol;
		}
	return lb;
	}

static void search_tree
   (const nlrctx*	c,
	const nlrnode*	node,
	u32				id,
	const nlrblock*	blk,
	u32				i,
	nlrbest*		best)
	{
	const nlrnode*	n;

	if (id == noIx) return;
	n = &node[id];

	// 子树里没有 qStart、tStart 都比 i 小的块，或没有已处理的块

	if ((n->minQ >= blk[i].span.qStart) || (n->minT >= blk[i].span.tStart)) return;
	if (n->maxTotal == -DBL_MAX) return;

	// 上界：total_j ≤ maxTotal，connectCost ≥ cost_lower_bound；严格小于才剪

	if (n->maxTotal + blk[i].score - cost_lower_bound (c, n, &blk[i]) < best->score)
		return;

	if (is_predecessor (&blk[n->blk], &blk[i]))
		offer (best, chain_candidate (c, blk, n->blk, i), n->blk);

	search_tree (c, node, n->hi, blk, i, best);
	search_tree (c, node, n->lo, blk, i, best);
	}

// chain_dp--
//	对一组块（同一 target 分区、query 分区，已按 DP 次序排好）求每块的 total
//	和 pred。

static void chain_dp
   (const nlrctx*	c,
	nlrblock*		blk,
	u32				n)
	{
	nlrnode*	node;
	u32*		perm;
	u32*		nodeOf;
	u32			i, numNodes, root;
	nlrbest		best;

	if (n == 0) return;

	node   = (nlrnode*) malloc_or_die ("nlrechain node",   n * sizeof(nlrnode));
	perm   = (u32*)     malloc_or_die ("nlrechain perm",   n * sizeof(u32));
	nodeOf = (u32*)     malloc_or_die ("nlrechain nodeOf", n * sizeof(u32));

	for (i=0 ; i<n ; i++) perm[i] = i;
	numNodes = 0;
	root = build_tree (blk, perm, 0, n, 0, noIx, node, &numNodes, nodeOf);

	for (i=0 ; i<n ; i++)
		{
		best.score = 0;
		best.pred  = noIx;
		search_tree (c, node, root, blk, i, &best);
		settle   (&blk[i], &best);
		activate (node, nodeOf, blk, i);
		}

	free_if_valid ("nlrechain node",   node);
	free_if_valid ("nlrechain perm",   perm);
	free_if_valid ("nlrechain nodeOf", nodeOf);
	}

//----------
//
// 拆块、合并完全重复块
//
//----------

// split_alignment--
//	§2 第 1 步：edit script 中每段连续的 sub 是一个块。追加到 *blk。

static void split_alignment
   (seq*		target,
	seq*		query,
	alignel*	a,
	u32			src,
	nlrblock**	blk,
	u32*		len,
	u32*		size)
	{
	editscript*	script = a->script;
	unspos		i = 0, j = 0, run;
	u32			opIx = 0, opBegin, ord = 0, rpt;
	nlrblock*	b;

	while (opIx < script->len)
		{
		opBegin = opIx;
		run = edit_script_run_of_subs (script, &opIx);
		if (run > 0)
			{
			if (*len >= *size)
				{
				*size = (*size == 0)? 1024 : 2 * (*size);
				*blk  = (nlrblock*) realloc_or_die ("nlrechain blocks", *blk,
				                                   (*size) * sizeof(nlrblock));
				}
			b = &(*blk)[(*len)++];
			b->span.tStart = a->beg1 - 1 + i;
			b->span.tEnd   = b->span.tStart + run;
			b->span.qStart = a->beg2 - 1 + j;
			b->span.qEnd   = b->span.qStart + run;
			b->part1   = (target->partition.p == NULL)? 0
			           : (u32) (lookup_partition (target, b->span.tStart) - target->partition.p);
			b->part2   = (query->partition.p == NULL)? 0
			           : (u32) (lookup_partition (query, b->span.qStart) - query->partition.p);
			b->src     = src;
			b->ord     = ord++;
			b->opBegin = opBegin;
			b->opAfter = opIx;
			}
		i += run;  j += run;

		for ( ; opIx<script->len ; opIx++)
			{
			rpt = edit_op_repeat (script->op[opIx]);
			switch (edit_op_operation (script->op[opIx]))
				{
				case editopIns: j += rpt;  continue;
				case editopDel: i += rpt;  continue;
				}
			break;						// （sub：回到外层）
			}
		}

	if ((i != a->end1 - a->beg1 + 1) || (j != a->end2 - a->beg2 + 1))
		suicidef ("nlrechain: edit script length mismatch (" unsposFmt "," unsposFmt ")", i, j);
	}

// cmp_blocks--
//	按 (target 分区, query 分区, tStart, qStart) 排序；起点相同时长的在前，再按
//	来源下标（X1）。合并后这也就是 DP 次序。

static int cmp_blocks
   (const void*	va,
	const void*	vb)
	{
	const nlrblock* a = (const nlrblock*) va;
	const nlrblock* b = (const nlrblock*) vb;

	if (a->part1 != b->part1) return (a->part1 < b->part1)? -1 : 1;
	if (a->part2 != b->part2) return (a->part2 < b->part2)? -1 : 1;
	if (a->span.tStart != b->span.tStart) return (a->span.tStart < b->span.tStart)? -1 : 1;
	if (a->span.qStart != b->span.qStart) return (a->span.qStart < b->span.qStart)? -1 : 1;
	if (a->span.tEnd   != b->span.tEnd)   return (a->span.tEnd   > b->span.tEnd)?   -1 : 1;
	if (a->src != b->src) return (a->src < b->src)? -1 : 1;
	return (a->ord < b->ord)? -1 : (a->ord > b->ord)? 1 : 0;
	}

// merge_exact_overlaps--
//	§2 第 3 步：起点完全相同的块并成一个，终点取较大者。块已按 cmp_blocks 排好，
//	所以每组的第一个就是最长的（等长取来源下标小的），保留它即可（X1）。返回合并
//	后的块数。

static u32 merge_exact_overlaps
   (nlrblock*	blk,
	u32			n)
	{
	u32	k, m;

	if (n == 0) return 0;
	for (k=1,m=1 ; k<n ; k++)
		{
		if ((blk[k].part1 == blk[m-1].part1) && (blk[k].part2 == blk[m-1].part2)
		 && (blk[k].span.tStart == blk[m-1].span.tStart)
		 && (blk[k].span.qStart == blk[m-1].span.qStart))
			{
			nlrechain_add_stat (bpMerged, blk[k].span.tEnd - blk[k].span.tStart);
			continue;
			}
		blk[m++] = blk[k];
		}
	return m;
	}

// score_blocks--
//	块分和块内最低单列分。

static void score_blocks
   (const nlrctx*	c,
	nlrblock*		blk,
	u32				n)
	{
	u32		k;
	unspos	p, len;
	double	col;

	for (k=0 ; k<n ; k++)
		{
		blk[k].score  = 0;
		blk[k].minCol = DBL_MAX;
		len = blk[k].span.tEnd - blk[k].span.tStart;
		for (p=0 ; p<len ; p++)
			{
			col = column_score (c, blk[k].span.tStart+p, blk[k].span.qStart+p);
			blk[k].score += col;
			if (col < blk[k].minCol) blk[k].minCol = col;
			}
		}
	}

//----------
//
// 剥链、链内整理、链分
//
//----------

// cmp_by_total--
//	剥链次序：total 降序，平局按 DP 次序。

static int cmp_by_total
   (const void*	va,
	const void*	vb)
	{
	u32 a = *(const u32*) va, b = *(const u32*) vb;
	if (sortBlocks[a].total != sortBlocks[b].total)
		return (sortBlocks[a].total > sortBlocks[b].total)? -1 : 1;
	return (a < b)? -1 : (a > b)? 1 : 0;
	}

// peel_chains--
//	§2 第 5 步：按 total 降序，从每个未用的块沿前驱回溯，遇到无前驱或前驱已被用
//	掉时停止。chainBlk 依次存放各链的块（链内按升序），chainStart[k] 为第 k 条链
//	的起始下标，chainStart[numChains] = n。返回链数。

static u32 peel_chains
   (const nlrblock*	blk,
	u32				n,
	u32*			chainBlk,
	u32*			chainStart)
	{
	u32*	order;
	char*	hit;
	u32		k, i, b, numChains = 0, len = 0, lo, hi, tmp;

	order = (u32*)  malloc_or_die ("nlrechain order", n * sizeof(u32));
	hit   = (char*) zalloc_or_die ("nlrechain hit",   n);
	for (k=0 ; k<n ; k++) order[k] = k;
	sortBlocks = blk;
	qsort (order, n, sizeof(u32), cmp_by_total);

	for (k=0 ; k<n ; k++)
		{
		i = order[k];
		if (hit[i]) continue;
		chainStart[numChains++] = len;
		lo = len;
		for (b=i ; ; b=blk[b].pred)
			{
			hit[b] = true;
			chainBlk[len++] = b;
			if ((blk[b].pred == noIx) || (hit[blk[b].pred])) break;
			}
		for (hi=len-1 ; lo<hi ; lo++,hi--)		// 回溯得到的是降序，反转
			{ tmp = chainBlk[lo];  chainBlk[lo] = chainBlk[hi];  chainBlk[hi] = tmp; }
		}
	chainStart[numChains] = len;

	free_if_valid ("nlrechain order", order);
	free_if_valid ("nlrechain hit",   hit);
	return numChains;
	}

// remove_partial_overlaps--
//	§2 第 6.1 步：链内相邻块 a、b 重叠 ov 列时，任一块被完全包含则删掉 b；否则在
//	交叉点 x 切开（b 起点 +x，a 终点 -(ov-x)），b 裁空则删掉 b，a 裁空则删掉所有
//	空块后从头再扫。返回剩下的块数。

static u32 remove_partial_overlaps
   (const nlrctx*	c,
	nlrpiece*		p,
	u32				n)
	{
	u32		a, b, k, m;
	s64		dq, dt, ov;
	unspos	x;
	int		trimmedOutA;

	if (n <= 1) return n;

	for (;;)
		{
		trimmedOutA = false;
		for (a=0,b=1 ; b<n ; b++)
			{
			dq = ((s64) p[b].span.qStart) - ((s64) p[a].span.qEnd);
			dt = ((s64) p[b].span.tStart) - ((s64) p[a].span.tEnd);
			if ((dq < 0) || (dt < 0))
				{
				ov = (dq < dt)? -dq : -dt;
				if ((ov >= (s64) (p[a].span.qEnd - p[a].span.qStart))
				 || (ov >= (s64) (p[b].span.qEnd - p[b].span.qStart)))
					{ p[b].span.qEnd = p[b].span.qStart;  continue; }	// 删掉 b

				find_crossover (c, &p[a].span, &p[b].span, (unspos) ov, &x);
				p[b].span.qStart += x;           p[b].span.tStart += x;
				p[a].span.qEnd   -= (unspos) ov - x;  p[a].span.tEnd -= (unspos) ov - x;
				if (p[b].span.qEnd <= p[b].span.qStart) continue;		// 删掉 b
				if (p[a].span.qEnd <= p[a].span.qStart) { trimmedOutA = true;  break; }
				}
			if (b != a+1) p[a+1] = p[b];	// 压掉 a、b 之间已删的块
			a = a+1;
			}
		if (!trimmedOutA) { n = a+1;  break; }

		// a 被裁空：把 [0,a] 之后未扫的块接上，删掉所有空块，从头再扫

		for (k=b ; k<n ; k++) p[a+1+k-b] = p[k];
		n = a+1 + n-b;
		for (k=0,m=0 ; k<n ; k++)
			if (p[k].span.qEnd > p[k].span.qStart) p[m++] = p[k];
		n = m;
		if (n <= 1) break;
		}
	return n;
	}

// chain_score--
//	§2 第 6.3 步：块分之和减去相邻块之间的 gapCalcCost。首尾正好相接的块合并
//	（第 6.2 步）不改变链分，因为 gapCalcCost(0,0) = 0，块分可加。

static double chain_score
   (const nlrctx*	c,
	const nlrpiece*	p,
	u32				n)
	{
	u32		k;
	double	s = 0;

	for (k=0 ; k<n ; k++)
		{
		s += span_score (c, &p[k].span);
		if (k > 0)
			s -= (double) gap_calc_cost (c->gapCost,
			                 ((s64) p[k].span.qStart) - ((s64) p[k-1].span.qEnd),
			                 ((s64) p[k].span.tStart) - ((s64) p[k-1].span.tEnd));
		}
	return s;
	}

//----------
//
// 报告（kent chain 格式，X6）
//
//----------

static u32 nextChainId = 1;

// part_info--
//	分区（或整条序列）的名字、在 seq->v 中的起点、长度、1 起的起始位置和真实
//	长度。

static void part_info
   (seq*			s,
	u32				partIx,
	const char*		defaultName,
	const char**	name,
	unspos*			offset,
	unspos*			len,
	unspos*			startLoc,
	unspos*			trueLen)
	{
	seqpartition*	sp = &s->partition;
	partition*		part;

	if (sp->p == NULL)
		{
		*name = (s->useFullNames)? s->header : s->shortHeader;
		if ((*name == NULL) || ((*name)[0] == 0)) *name = defaultName;
		*offset   = 0;
		*len      = s->len;
		*startLoc = s->startLoc;
		*trueLen  = s->trueLen;
		}
	else
		{
		part      = &sp->p[partIx];
		*name     = &sp->pool[part->header];
		*offset   = part->sepBefore + 1;
		*len      = part->sepAfter - part->sepBefore - 1;
		*startLoc = part->startLoc;
		*trueLen  = part->trueLen;
		}
	}

// write_chain--
//	t 为正链坐标；q 在反链时为整条序列反向互补后的坐标（同 kent）。首尾正好相接
//	的块合并后输出。

static void write_chain
   (FILE*			f,
	const nlrctx*	c,
	const nlrpiece*	p,
	u32				n,
	double			chainScore,
	u32				part1,
	u32				part2)
	{
	const char*		name1, *name2;
	unspos			off1, len1, loc1, size1, off2, len2, loc2, size2;
	unspos			tb, te, qb, qe, size;
	int				minus = ((c->query->revCompFlags & rcf_rev) != 0);
	u32				k;

	if (n == 0) return;
	part_info (c->target, part1, "seq1", &name1, &off1, &len1, &loc1, &size1);
	part_info (c->query,  part2, "seq2", &name2, &off2, &len2, &loc2, &size2);

	tb = p[0].span.tStart - off1 + loc1 - 1;
	te = p[n-1].span.tEnd - off1 + loc1 - 1;
	qb = p[0].span.qStart - off2;
	qe = p[n-1].span.qEnd - off2;
	if (!minus) { qb += loc2 - 1;  qe += loc2 - 1; }
	else        { qb += size2 - (loc2-1+len2);  qe += size2 - (loc2-1+len2); }

	fprintf (f, "chain %.0f %s " unsposFmt " + " unsposFmt " " unsposFmt
	            " %s " unsposFmt " %c " unsposFmt " " unsposFmt " %u\n",
	         chainScore, name1, size1, tb, te,
	         name2, size2, (minus)? '-' : '+', qb, qe, nextChainId++);

	size = p[0].span.tEnd - p[0].span.tStart;
	for (k=1 ; k<n ; k++)
		{
		if ((p[k].span.tStart == p[k-1].span.tEnd) && (p[k].span.qStart == p[k-1].span.qEnd))
			{ size += p[k].span.tEnd - p[k].span.tStart;  continue; }
		fprintf (f, unsposFmt "\t" unsposFmt "\t" unsposFmt "\n", size,
		         p[k].span.tStart - p[k-1].span.tEnd, p[k].span.qStart - p[k-1].span.qEnd);
		size = p[k].span.tEnd - p[k].span.tStart;
		}
	fprintf (f, unsposFmt "\n\n", size);
	}

//----------
//
// 输出重建（D2(a)、X2、X3）
//
//----------

static int cmp_runs
   (const void*	va,
	const void*	vb)
	{
	const nlrrun* a = (const nlrrun*) va;
	const nlrrun* b = (const nlrrun*) vb;
	if (a->src != b->src) return (a->src < b->src)? -1 : 1;
	return (a->ord < b->ord)? -1 : (a->ord > b->ord)? 1 : 0;
	}

// joins_previous--
//	X2：piece k 与链中前一块 k-1 来自同一比对、次序相邻，且接缝两侧都没被裁。
//	（同一比对中相邻的块互不重叠，链内整理不会裁这个接缝，后两个条件实际总是
//	成立；保留它们是为了让"沿用原比对的 gap"这一前提写在明处。）

static int joins_previous
   (const nlrblock*	blk,
	const nlrpiece*	prev,
	const nlrpiece*	cur)
	{
	const nlrblock* a = &blk[prev->blk];
	const nlrblock* b = &blk[cur->blk];
	return (a->src == b->src) && (b->ord == a->ord+1)
	    && (prev->span.tEnd   == a->span.tEnd)
	    && (cur->span.tStart  == b->span.tStart);
	}

// build_alignment--
//	把一段（kept[first..first+count)）构造成新的 alignel：块之间沿用来源比对的
//	indel 操作，分数用 score_alignment 重算。

static alignel* build_alignment
   (const nlrctx*	c,
	alignel**		srcAlign,
	const nlrblock*	blk,
	const nlrpiece*	kept,
	const nlrrun*	run)
	{
	const nlrpiece*	p = &kept[run->first];
	const alignel*	srcA = srcAlign[run->src];
	editscript*		script;
	u32				k, opIx;
	alignel*		a;

	script = edit_script_new ();
	edit_script_sub (&script, p[0].span.tEnd - p[0].span.tStart);
	for (k=1 ; k<run->count ; k++)
		{
		for (opIx=blk[p[k-1].blk].opAfter ; opIx<blk[p[k].blk].opBegin ; opIx++)
			if (edit_op_repeat (srcA->script->op[opIx]) > 0)
				edit_script_add (&script, edit_op_operation (srcA->script->op[opIx]),
				                 edit_op_repeat (srcA->script->op[opIx]));
		edit_script_sub (&script, p[k].span.tEnd - p[k].span.tStart);
		}

	a = (alignel*) malloc_or_die ("nlrechain alignel", sizeof(alignel));
	a->next      = NULL;
	a->isTrivial = srcA->isTrivial;
	a->beg1      = p[0].span.tStart + 1;
	a->beg2      = p[0].span.qStart + 1;
	a->end1      = p[run->count-1].span.tEnd;
	a->end2      = p[run->count-1].span.qEnd;
	a->seq1      = srcA->seq1;
	a->seq2      = srcA->seq2;
	a->script    = script;
	a->hspId     = srcA->hspId;
	a->s         = score_alignment (c->scoring, c->target->v, a->beg1-1,
	                                c->query->v, a->beg2-1, script);
	return a;
	}

//----------
//
// nl_rechain_aligns--
//	§2 第 1–8 步。
//
//----------
//
// Arguments:
//	seq*		target:		The target sequence.
//	seq*		query:		The query sequence (reverse-complemented for the
//							.. minus strand).
//	scoreset*	scoring:	The substitution scores (same as gapped extension).
//	alignel*	alignList:	The gapped alignments for this strand;  freed here.
//	const nlrechainparams* params:	T_chain, gap table, report file.
//
// Returns:
//	The replacement alignment list (possibly NULL).
//
//----------

alignel* nl_rechain_aligns
   (seq*					target,
	seq*					query,
	scoreset*				scoring,
	alignel*				alignList,
	const nlrechainparams*	params)
	{
	nlrctx		c;
	alignel**	srcAlign = NULL;
	alignel*	a, *head = NULL, *tail = NULL;
	nlrblock*	blk = NULL;
	nlrpiece*	kept = NULL, *p = NULL;
	nlrrun*		runs = NULL;
	u32*		chainBlk = NULL, *chainStart = NULL;
	u32			numAligns, numBlocks, blkSize, numKept, numRuns;
	u32			g0, g1, ch, numChains, k, n, m;
	double		s;

	if (alignList == NULL) return NULL;

	c.target  = target;
	c.query   = query;
	c.scoring = scoring;
	c.gapCost = params->gapCost;

	// 第 1 步：拆块

	for (numAligns=0,a=alignList ; a!=NULL ; a=a->next) numAligns++;
	srcAlign = (alignel**) malloc_or_die ("nlrechain srcAlign", numAligns * sizeof(alignel*));
	numBlocks = blkSize = 0;
	for (k=0,a=alignList ; a!=NULL ; k++,a=a->next)
		{
		srcAlign[k] = a;
		split_alignment (target, query, a, k, &blk, &numBlocks, &blkSize);
		}
	nlrechain_add_stat (numAlignsIn, numAligns);
#ifdef collect_stats
	for (k=0 ; k<numBlocks ; k++)
		nlrechain_add_stat (bpIn, blk[k].span.tEnd - blk[k].span.tStart);
#endif // collect_stats

	// 第 3 步：合并完全重复的块；第 2 步：块分

	if (numBlocks > 0)
		qsort (blk, numBlocks, sizeof(nlrblock), cmp_blocks);
	numBlocks = merge_exact_overlaps (blk, numBlocks);
	score_blocks (&c, blk, numBlocks);
	nlrechain_add_stat (numBlocks, numBlocks);

	chainBlk   = (u32*)      malloc_or_die ("nlrechain chainBlk",   (numBlocks+1) * sizeof(u32));
	chainStart = (u32*)      malloc_or_die ("nlrechain chainStart", (numBlocks+1) * sizeof(u32));
	p          = (nlrpiece*) malloc_or_die ("nlrechain pieces",     (numBlocks+1) * sizeof(nlrpiece));
	kept       = (nlrpiece*) malloc_or_die ("nlrechain kept",       (numBlocks+1) * sizeof(nlrpiece));
	runs       = (nlrrun*)   malloc_or_die ("nlrechain runs",       (numBlocks+1) * sizeof(nlrrun));
	numKept = numRuns = 0;

	// 按 (target 分区, query 分区) 分组，每组做第 4–7 步

	for (g0=0 ; g0<numBlocks ; g0=g1)
		{
		for (g1=g0+1 ; g1<numBlocks ; g1++)
			if ((blk[g1].part1 != blk[g0].part1) || (blk[g1].part2 != blk[g0].part2)) break;

		chain_dp (&c, blk+g0, g1-g0);
		numChains = peel_chains (blk+g0, g1-g0, chainBlk, chainStart);
		nlrechain_add_stat (numChains, numChains);

		for (ch=0 ; ch<numChains ; ch++)
			{
			n = chainStart[ch+1] - chainStart[ch];
			for (k=0 ; k<n ; k++)
				{
				p[k].blk  = g0 + chainBlk[chainStart[ch]+k];
				p[k].span = blk[p[k].blk].span;
				nlrechain_add_stat (bpTrimmed, p[k].span.tEnd - p[k].span.tStart);
				}
			n = remove_partial_overlaps (&c, p, n);
			s = chain_score (&c, p, n);

			if (params->reportFile != NULL)
				write_chain (params->reportFile, &c, p, n, s, blk[g0].part1, blk[g0].part2);

#ifdef collect_stats
			for (k=0 ; k<n ; k++)
				{
				nlrechainStats.bpTrimmed -= p[k].span.tEnd - p[k].span.tStart;
				if (s >= params->minChainScore)
					nlrechain_add_stat (bpKept,    p[k].span.tEnd - p[k].span.tStart);
				else
					nlrechain_add_stat (bpRemoved, p[k].span.tEnd - p[k].span.tStart);
				}
#endif // collect_stats

			// 第 7 步：下限

			if (s < params->minChainScore)
				{ nlrechain_count_stat (numChainsRemoved);  continue; }

			// 第 8 步：按 X2 分段

			for (k=0 ; k<n ; k++)
				{
				kept[numKept] = p[k];
				if ((k > 0) && (joins_previous (blk, &p[k-1], &p[k])))
					runs[numRuns-1].count++;
				else
					{
					runs[numRuns].first = numKept;
					runs[numRuns].count = 1;
					runs[numRuns].src   = blk[p[k].blk].src;
					runs[numRuns].ord   = blk[p[k].blk].ord;
					numRuns++;
					}
				numKept++;
				}
			}
		}

	// 第 8 步：按来源比对和次序排序（X3），构造新列表

	if (numRuns > 0)
		qsort (runs, numRuns, sizeof(nlrrun), cmp_runs);
	for (m=0 ; m<numRuns ; m++)
		{
		a = build_alignment (&c, srcAlign, blk, kept, &runs[m]);
		if (head == NULL) head = a;
		             else tail->next = a;
		tail = a;
		}
	nlrechain_add_stat (numAlignsOut, numRuns);

	free_align_list (alignList);
	free_if_valid ("nlrechain srcAlign",   srcAlign);
	free_if_valid ("nlrechain blocks",     blk);
	free_if_valid ("nlrechain chainBlk",   chainBlk);
	free_if_valid ("nlrechain chainStart", chainStart);
	free_if_valid ("nlrechain pieces",     p);
	free_if_valid ("nlrechain kept",       kept);
	free_if_valid ("nlrechain runs",       runs);
	return head;
	}

//----------
//
// nlrechain_zero_stats, nlrechain_show_stats, nlrechain_generic_stats--
//	统计（§3）。
//
//----------

void nlrechain_zero_stats
   (void)
	{
#ifdef collect_stats
	memset (&nlrechainStats, 0, sizeof(nlrechainStats));
#endif // collect_stats
	}

void nlrechain_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	if (nlrechainStats.numAlignsIn == 0) return;
	fprintf (f, "  rechain aligns in: %s\n", ucommatize(nlrechainStats.numAlignsIn));
	fprintf (f, " rechain aligns out: %s\n", ucommatize(nlrechainStats.numAlignsOut));
	fprintf (f, "     rechain blocks: %s\n", ucommatize(nlrechainStats.numBlocks));
	fprintf (f, "     rechain chains: %s\n", ucommatize(nlrechainStats.numChains));
	fprintf (f, "     chains removed: %s\n", ucommatize(nlrechainStats.numChainsRemoved));
	fprintf (f, "   rechain bp input: %s\n", ucommatize(nlrechainStats.bpIn));
	fprintf (f, "  rechain bp merged: %s\n", ucommatize(nlrechainStats.bpMerged));
	fprintf (f, " rechain bp trimmed: %s\n", ucommatize(nlrechainStats.bpTrimmed));
	fprintf (f, " rechain bp removed: %s\n", ucommatize(nlrechainStats.bpRemoved));
	fprintf (f, "    rechain bp kept: %s\n", ucommatize(nlrechainStats.bpKept));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void nlrechain_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "nlrechain_aligns_in=%" PRIu64 "\n",      nlrechainStats.numAlignsIn);
	(*func) (f, "nlrechain_aligns_out=%" PRIu64 "\n",     nlrechainStats.numAlignsOut);
	(*func) (f, "nlrechain_blocks=%" PRIu64 "\n",         nlrechainStats.numBlocks);
	(*func) (f, "nlrechain_chains=%" PRIu64 "\n",         nlrechainStats.numChains);
	(*func) (f, "nlrechain_chains_removed=%" PRIu64 "\n", nlrechainStats.numChainsRemoved);
	(*func) (f, "nlrechain_bp_in=%" PRIu64 "\n",          nlrechainStats.bpIn);
	(*func) (f, "nlrechain_bp_merged=%" PRIu64 "\n",      nlrechainStats.bpMerged);
	(*func) (f, "nlrechain_bp_trimmed=%" PRIu64 "\n",     nlrechainStats.bpTrimmed);
	(*func) (f, "nlrechain_bp_removed=%" PRIu64 "\n",     nlrechainStats.bpRemoved);
	(*func) (f, "nlrechain_bp_kept=%" PRIu64 "\n",        nlrechainStats.bpKept);
#endif // collect_stats
	}
