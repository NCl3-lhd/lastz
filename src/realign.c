//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: realign.c
//
//----------
//
// realign--
//	输出前的框内全局重比对：每条 gapped 比对沿原路径切段，段的两端固定在原路
//	径上，段内用单仿射 gap 做带约束的全局比对（Gotoh），再拼回 edit script。比
//	对的首尾坐标和参与的碱基不变，只改 gap 的摆放。
//
// 设计与依据见 docs/realign_plan.md（下文 §x 指该文档章节）与
// docs/realign_ksw2_eval.md（服务器上 Python + ksw2 原型的结果）。
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
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "gapped_extend.h"		// score_alignment

#define  realign_owner			// (make this the owner of its globals)
#include "realign.h"			// interface to this module

//----------
//
// 切段与带宽的常数（§2 D6，与原型相同）
//
//----------

#define realignChunkCols  1000	// 相邻切点至少相隔的列数
#define realignCutRun     11	// 切点所在无 gap 段的最短长度
#define realignBandMargin 16	// 带宽 = 原路径最大对角线偏移 + 此余量
#define realignMaxCells   (((u64) 1) << 30)	// 一段的回溯表超过此格数时保留原路径

// 每列的操作（展开后的 edit script）

#define colSub 0
#define colIns 1				// query（seq2）多出的碱基
#define colDel 2				// target（seq1）多出的碱基

// 回溯表每格的位

#define tbHFromDiag 0
#define tbHFromE    1
#define tbHFromF    2
#define tbHMask     3
#define tbEExtend   4			// E（插入）由 E 延伸而来，否则由 H 开 gap
#define tbFExtend   8			// F（删除）由 F 延伸而来

#define negInf ((score) -1000000000)

//----------
//
// 模块内的工作区（反复使用，按需扩大）
//
//----------

static u8*		colOp    = NULL;	// 一条比对展开后的每列操作
static u64		colOpCap = 0;
static u8*		newOp    = NULL;	// 一段重比对结果（逆序）
static u64		newOpCap = 0;
static u8*		tb       = NULL;	// 回溯表，(lt+1) 行 × (2w+1) 列
static u64		tbCap    = 0;
static score*	hPrev    = NULL;	// 按 query 下标的整行数组
static score*	hCur     = NULL;
static score*	fPrev    = NULL;
static score*	fCur     = NULL;
static u64		rowCap   = 0;

static void grow_u8
   (u8**		buf,
	u64*		cap,
	u64			need,
	char*		name)
	{
	if (need <= *cap) return;
	*cap = need + need/2 + 16;
	*buf = (u8*) realloc_or_die (name, *buf, *cap);
	}

static void grow_rows
   (u64			need)
	{
	if (need <= rowCap) return;
	rowCap = need + need/2 + 16;
	hPrev = (score*) realloc_or_die ("realign hPrev", hPrev, rowCap * sizeof(score));
	hCur  = (score*) realloc_or_die ("realign hCur",  hCur,  rowCap * sizeof(score));
	fPrev = (score*) realloc_or_die ("realign fPrev", fPrev, rowCap * sizeof(score));
	fCur  = (score*) realloc_or_die ("realign fCur",  fCur,  rowCap * sizeof(score));
	}

//----------
//
// path_score--
//	一段列操作在重比对打分下的分：替换分取 sub[target][query]，每个连续的插入
//	或删除段罚 gapOpen + l*gapExtend；相邻的插入段与删除段各算一个 gap。
//
//----------

static score path_score
   (const scoreset*		scoring,
	const realignparams* p,
	const u8*			t,			// 段内 target 起点
	const u8*			q,			// 段内 query 起点
	const u8*			op,
	u64					n)
	{
	score	s = 0;
	u64		k, run;

	for (k=0 ; k<n ; )
		{
		if (op[k] == colSub)
			{ s += scoring->sub[*(t++)][*(q++)];  k++;  continue; }
		run = 1;
		while ((k+run < n) && (op[k+run] == op[k])) run++;
		s -= p->gapOpen + ((score) run) * p->gapExtend;
		if (op[k] == colIns) q += run;
		                else t += run;
		k += run;
		}
	return s;
	}

//----------
//
// max_drift--
//	一段列操作走过的最大对角线偏移 |i-j|（i、j 为段内已消耗的 target、query
//	碱基数）。
//
//----------

static u64 max_drift
   (const u8*	op,
	u64			n)
	{
	s64		d = 0, maxD = 0;
	u64		k;

	for (k=0 ; k<n ; k++)
		{
		if      (op[k] == colIns) d--;
		else if (op[k] == colDel) d++;
		if ( d > maxD) maxD =  d;
		if (-d > maxD) maxD = -d;
		}
	return (u64) maxD;
	}

//----------
//
// banded_global--
//	t[0..lt) 对 q[0..lq) 的全局比对，单仿射 gap，带约束 |i-j| <= w。结果（逆
//	序的列操作）写入 newOp，返回列数；*best 为最优分。
//
// 平局（§2 D2）：rightAlign 时回溯中优先走 gap、gap 内优先"开"，indel 因此靠
// 右；否则优先走对角线、gap 内优先"延伸"，indel 靠左。
//
//----------

static u64 banded_global
   (const scoreset*		scoring,
	const realignparams* p,
	const u8*			t,
	u64					lt,
	const u8*			q,
	u64					lq,
	u64					w,
	score*				best)
	{
	const score	open = p->gapOpen + p->gapExtend;
	const score	ext  = p->gapExtend;
	const int	right = p->rightAlign;
	u64			width = 2*w + 1;
	u64			i, j, lo, hi, n;
	score*		swap;
	score		h, e, f, d, hLeft, eLeft, x, y;
	u8			bits;
	int			state;

	grow_u8 (&tb, &tbCap, (lt+1) * width, "realign traceback");
	grow_rows (lq + 2);
	for (j=0 ; j<=lq+1 ; j++)
		hPrev[j] = hCur[j] = fPrev[j] = fCur[j] = negInf;

	for (i=0 ; i<=lt ; i++)
		{
		lo = (i > w)? i - w : 0;
		hi = (i + w < lq)? i + w : lq;
		hLeft = eLeft = negInf;
		for (j=lo ; j<=hi ; j++)
			{
			bits = 0;

			// E(i,j)：消耗 query[j-1]，来自 (i,j-1)

			e = negInf;
			if (j > lo)
				{
				x = hLeft - open;
				y = eLeft - ext;
				if ((right)? (y > x) : (y >= x)) { e = y;  bits |= tbEExtend; }
				                            else   e = x;
				if (e < negInf) e = negInf;
				}

			// F(i,j)：消耗 target[i-1]，来自 (i-1,j)；j = hi(i-1)+1 时 (i-1,j) 在
			// 上一行的带外，读到的是上一行写在 hi+1 处的 negInf

			f = negInf;
			if (i > 0)
				{
				x = hPrev[j] - open;
				y = fPrev[j] - ext;
				if ((right)? (y > x) : (y >= x)) { f = y;  bits |= tbFExtend; }
				                            else   f = x;
				if (f < negInf) f = negInf;
				}

			// H(i,j)

			if ((i == 0) && (j == 0))
				h = 0;
			else
				{
				d = ((i > 0) && (j > 0))? hPrev[j-1] + scoring->sub[t[i-1]][q[j-1]]
				                        : negInf;
				if (d < negInf) d = negInf;
				h = d;  bits |= tbHFromDiag;
				if ((right)? (e >= h) : (e > h)) { h = e;  bits = (bits & ~tbHMask) | tbHFromE; }
				if ((right)? (f >= h) : (f > h)) { h = f;  bits = (bits & ~tbHMask) | tbHFromF; }
				}

			hCur[j] = h;  fCur[j] = f;
			hLeft = h;    eLeft = e;
			tb[i*width + (j + w - i)] = bits;
			}
		if (hi + 1 <= lq) { hCur[hi+1] = negInf;  fCur[hi+1] = negInf; }
		swap = hPrev;  hPrev = hCur;  hCur = swap;
		swap = fPrev;  fPrev = fCur;  fCur = swap;
		}
	*best = hPrev[lq];

	// 回溯（从 (lt,lq) 的 H 出发），列操作逆序写入 newOp

	grow_u8 (&newOp, &newOpCap, lt + lq + 1, "realign newOp");
	n = 0;
	i = lt;  j = lq;  state = 0;		// 0 = H, 1 = E, 2 = F
	while ((i > 0) || (j > 0))
		{
		bits = tb[i*width + (j + w - i)];
		if (state == 0)
			{
			switch (bits & tbHMask)
				{
				case tbHFromDiag: newOp[n++] = colSub;  i--;  j--;  break;
				case tbHFromE:    state = 1;                       break;
				default:          state = 2;                       break;
				}
			}
		else if (state == 1)
			{
			newOp[n++] = colIns;  j--;
			if (!(bits & tbEExtend)) state = 0;
			}
		else
			{
			newOp[n++] = colDel;  i--;
			if (!(bits & tbFExtend)) state = 0;
			}
		}
	return n;
	}

//----------
//
// realign_one--
//	重比对一条比对：展开、切段、逐段带约束全局比对、拼回新的 edit script。
//
//----------

static void realign_one
   (const seq*			target,
	const seq*			query,
	const scoreset*		scoring,
	alignel*			a,
	const realignparams* p)
	{
	editscript*	script = a->script;
	editscript*	newScript;
	u64			numCols, k, rpt, c, start, last, mid, c0, c1, n, lt, lq, w, wMax, drift;
	u64			ti, qi, ti0, qi0, tiCut, qiCut;
	u32			opIx;
	u8			code;
	score		orig, best;
	int			kept;

	// 展开 edit script

	numCols = 0;
	for (opIx=0 ; opIx<script->len ; opIx++)
		numCols += edit_op_repeat (script->op[opIx]);
	grow_u8 (&colOp, &colOpCap, numCols, "realign colOp");
	numCols = 0;
	for (opIx=0 ; opIx<script->len ; opIx++)
		{
		rpt = edit_op_repeat (script->op[opIx]);
		switch (edit_op_operation (script->op[opIx]))
			{
			case editopIns: code = colIns;  break;
			case editopDel: code = colDel;  break;
			default:        code = colSub;  break;
			}
		for (k=0 ; k<rpt ; k++) colOp[numCols++] = code;
		}
	realign_add_stat (colsIn, numCols);

	// 逐段处理：切点是长度 >= realignCutRun 的无 gap 段中点，相邻切点至少相隔
	// realignChunkCols 列（与原型 realign_maf.py 的 cut_points 相同）

	newScript = edit_script_new ();
	ti = a->beg1 - 1;  qi = a->beg2 - 1;	// 扫描位置
	ti0 = ti;  qi0 = qi;  c0 = 0;			// 当前段的起点
	last = 0;
	c = 0;
	for (;;)
		{
		// 找下一个切点（或比对末尾）

		c1 = numCols;  tiCut = qiCut = 0;
		while (c < numCols)
			{
			if (colOp[c] != colSub)
				{
				if (colOp[c] == colIns) qi++; else ti++;
				c++;
				continue;
				}
			start = c;
			while ((c < numCols) && (colOp[c] == colSub)) c++;
			ti += c - start;  qi += c - start;
			if (c - start < realignCutRun) continue;
			mid = (start + c) / 2;
			if ((mid - last >= realignChunkCols) && (numCols - mid >= realignCutRun/2))
				{
				c1 = mid;  last = mid;
				tiCut = ti - (c - mid);  qiCut = qi - (c - mid);
				break;
				}
			}
		if (c1 == numCols) { tiCut = a->end1;  qiCut = a->end2; }

		// 段 [c0,c1)：target [ti0,tiCut)，query [qi0,qiCut)

		lt = tiCut - ti0;  lq = qiCut - qi0;
		orig = path_score (scoring, p, target->v + ti0, query->v + qi0, colOp + c0, c1 - c0);
		drift = max_drift (colOp + c0, c1 - c0);
		w = ((lt > lq)? lt - lq : lq - lt);
		if (drift > w) w = drift;
		w += realignBandMargin;
		wMax = (lt > lq)? lt : lq;
		if (w > wMax) w = wMax;
		realign_count_stat (numChunks);

		kept = false;
		for (;;)
			{
			if ((lt+1) * (2*w+1) > realignMaxCells)
				{ kept = true;  break; }
			n = banded_global (scoring, p, target->v + ti0, lt, query->v + qi0, lq, w, &best);
			if (w >= wMax) break;
			drift = max_drift (newOp, n);	// (逆序不影响最大偏移)
			if (drift < w) break;
			w = 2*w + 16;
			if (w > wMax) w = wMax;
			realign_count_stat (numWidened);
			}

		if (kept)
			{
			for (k=c0 ; k<c1 ; k++)
				edit_script_add (&newScript, (colOp[k]==colSub)? editopSub
				                           : (colOp[k]==colIns)? editopIns : editopDel, 1);
			realign_add_stat (colsOut, c1 - c0);
			}
		else
			{
			if (best < orig)
				suicidef ("realign: chunk score " scoreFmtSimple " below original path score " scoreFmtSimple
				          " (target " unsposFmt ".." unsposFmt ", query " unsposFmt ".." unsposFmt ")",
				          best, orig, (unspos) ti0, (unspos) tiCut, (unspos) qi0, (unspos) qiCut);
			if (best > orig) { realign_count_stat (numImproved); }
			for (k=n ; k>0 ; )
				{
				code = newOp[--k];
				rpt = 1;
				while ((k > 0) && (newOp[k-1] == code)) { k--;  rpt++; }
				edit_script_add (&newScript, (code==colSub)? editopSub
				                           : (code==colIns)? editopIns : editopDel, rpt);
				}
			realign_add_stat (colsOut, n);
			}

		if (c1 == numCols) break;
		c0 = c1;  ti0 = tiCut;  qi0 = qiCut;
		c = c1;   ti  = tiCut;  qi  = qiCut;
		}

	// 新 script 必须恰好覆盖原比对的两侧

	ti = qi = 0;
	for (opIx=0 ; opIx<newScript->len ; opIx++)
		{
		rpt = edit_op_repeat (newScript->op[opIx]);
		switch (edit_op_operation (newScript->op[opIx]))
			{
			case editopIns: qi += rpt;              break;
			case editopDel: ti += rpt;              break;
			default:        ti += rpt;  qi += rpt;  break;
			}
		}
	if ((ti != a->end1 - a->beg1 + 1) || (qi != a->end2 - a->beg2 + 1))
		suicidef ("realign: new script length mismatch (target " unsposFmt ".." unsposFmt ")",
		          a->beg1, a->end1);

	free_if_valid ("realign old script", a->script);
	a->script = newScript;
	a->s      = score_alignment ((scoreset*) scoring, target->v, a->beg1-1,
	                             query->v, a->beg2-1, newScript);
	}

//----------
//
// realign_aligns--
//	对列表中的每条比对做框内重比对（§1）。
//
//----------
//
// Arguments:
//	seq*		target:		The target sequence.
//	seq*		query:		The query sequence (reverse-complemented for the
//							.. minus strand).
//	scoreset*	scoring:	lastz 自身的打分；替换分用于重比对，整体用于重算
//							.. 输出分数（§2 D4）。
//	alignel*	alignList:	即将输出的比对，原地修改。
//	const realignparams* params:	gap 罚分与平局方向。
//
// Returns:
//	alignList（同一个列表）。
//
//----------

alignel* realign_aligns
   (seq*					target,
	seq*					query,
	scoreset*				scoring,
	alignel*				alignList,
	const realignparams*	params)
	{
	alignel*	a;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		realign_count_stat (numAligns);
		realign_one (target, query, scoring, a, params);
		}
	return alignList;
	}

//----------
//
// realign_zero_stats, realign_show_stats, realign_generic_stats--
//
//----------

void realign_zero_stats
   (void)
	{
#ifdef collect_stats
	memset (&realignStats, 0, sizeof(realignStats));
#endif // collect_stats
	}

void realign_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	if (realignStats.numAligns == 0) return;
	fprintf (f, "     realign aligns: %s\n", ucommatize(realignStats.numAligns));
	fprintf (f, "     realign chunks: %s\n", ucommatize(realignStats.numChunks));
	fprintf (f, "    realign widened: %s\n", ucommatize(realignStats.numWidened));
	fprintf (f, "   realign improved: %s\n", ucommatize(realignStats.numImproved));
	fprintf (f, "    realign cols in: %s\n", ucommatize(realignStats.colsIn));
	fprintf (f, "   realign cols out: %s\n", ucommatize(realignStats.colsOut));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void realign_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "realign_aligns=%" PRIu64 "\n",   realignStats.numAligns);
	(*func) (f, "realign_chunks=%" PRIu64 "\n",   realignStats.numChunks);
	(*func) (f, "realign_widened=%" PRIu64 "\n",  realignStats.numWidened);
	(*func) (f, "realign_improved=%" PRIu64 "\n", realignStats.numImproved);
	(*func) (f, "realign_cols_in=%" PRIu64 "\n",  realignStats.colsIn);
	(*func) (f, "realign_cols_out=%" PRIu64 "\n", realignStats.colsOut);
#endif // collect_stats
	}
