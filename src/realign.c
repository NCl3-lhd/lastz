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
#include <math.h>				// exp, log（MEA）
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
// MEA（后验解码，docs/mea_plan.md）
//
// 段内路径分布取 Boltzmann 形式：P(path) ∝ exp(λ·S(path))，S 与 banded_global
// 的打分相同（替换分 + 单仿射 gap）。λ 很大时分布集中到 Viterbi 路径上。
//
// 数值：每条全局路径的对角线列数 d 与 gap 列数满足 2d + #gap = lt + lq，所以把
// 每个对角线列的权重乘 e^{-λ smax}、每个 gap 列乘 e^{-λ smax/2}，所有路径同乘
// 一个常数，分布不变；这样所有单步权重都 ≤ 1，不会上溢。前向、后向各按行缩放，
// 每行的对数缩放因子单独记下。
//
// 状态按最后一列区分：M（替换）、X（插入，消耗 query）、Y（删除，消耗 target）。
// 段首 (0,0) 视为 M 状态，权重 1。X 只能由 M/Y 开（×go）或由 X 延伸（×ge），Y
// 同理，与 banded_global 的路径集合和 path_score 的计分一致。
//
//----------

#define meaMaxCells (((u64) 1) << 26)	// 一段的 double 矩阵超过此格数时保留 Viterbi 路径
#define meaNegInf   (-1.0e300)

static double*	post      = NULL;	// 前向 M 值，随后原地改为匹配后验；(lt+1)×(2w+1)
static u64		postCap   = 0;
static double*	rowScale  = NULL;	// 前向每行的对数缩放（累计）
static u64		rowScaleCap = 0;
static double*	meaRow[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
static u64		meaRowCap = 0;
static u8*		vitOp     = NULL;	// Viterbi 路径（逆序）的副本
static u64		vitOpCap  = 0;

static double*	expSub    = NULL;	// exp(λ(sub[a][b] − smax))，256×256
static double	meaGo, meaGe;		// 开 gap 列、延伸 gap 列的权重（含 smax/2 的平移）
static const scoreset* expSubFor = NULL;	// expSub 对应的打分与参数
static double	expSubLambda = -1, expSubOpen = -1, expSubExtend = -1;

static void grow_double
   (double**	buf,
	u64*		cap,
	u64			need,
	char*		name)
	{
	if (need <= *cap) return;
	*cap = need + need/2 + 16;
	*buf = (double*) realloc_or_die (name, *buf, *cap * sizeof(double));
	}

static void grow_mea_rows
   (u64			need)
	{
	int			k;

	if (need <= meaRowCap) return;
	meaRowCap = need + need/2 + 16;
	for (k=0 ; k<6 ; k++)
		meaRow[k] = (double*) realloc_or_die ("realign meaRow", meaRow[k], meaRowCap * sizeof(double));
	}

//----------
//
// mea_lambda_ka--
//	打分矩阵（ACGT 部分）在均匀背景下的 Karlin-Altschul λ，即
//	Σ_{a,b} e^{λ s(a,b)} / 16 = 1 的正根。期望分不为负或没有正分时无解，退出。
//
//----------

static double mea_ka_f
   (const scoreset*	scoring,
	double			lambda)
	{
	static const u8 nts[4] = { 'A', 'C', 'G', 'T' };
	double		sum = 0;
	int			a, b;

	for (a=0 ; a<4 ; a++) for (b=0 ; b<4 ; b++)
		sum += exp (lambda * (double) scoring->sub[nts[a]][nts[b]]);
	return sum / 16 - 1;
	}

static double mea_lambda_ka
   (const scoreset*	scoring)
	{
	static const u8 nts[4] = { 'A', 'C', 'G', 'T' };
	double		mean = 0, lo, hi, mid;
	int			a, b, k, anyPos = false;

	for (a=0 ; a<4 ; a++) for (b=0 ; b<4 ; b++)
		{
		mean += (double) scoring->sub[nts[a]][nts[b]] / 16;
		if (scoring->sub[nts[a]][nts[b]] > 0) anyPos = true;
		}
	if ((mean >= 0) || (!anyPos))
		suicide ("--mea needs a scoring matrix with negative expected score and some positive score");

	lo = 0;  hi = 1e-6;
	for (k=0 ; (k<2000) && (mea_ka_f (scoring, hi) <= 0) ; k++)
		{ lo = hi;  hi *= 2; }
	for (k=0 ; k<200 ; k++)
		{
		mid = (lo + hi) / 2;
		if (mea_ka_f (scoring, mid) > 0) hi = mid;
		                            else lo = mid;
		}
	return (lo + hi) / 2;
	}

//----------
//
// mea_prepare--
//	按打分和参数建立 expSub、meaGo、meaGe（与上次相同时直接复用）。
//
//----------

static void mea_prepare
   (const scoreset*		scoring,
	const realignparams* p)
	{
	double		lambda, smax;
	int			a, b;

	lambda = p->meaLambdaScale * mea_lambda_ka (scoring);
	if ((expSubFor == scoring) && (expSubLambda == lambda)
	 && (expSubOpen == (double) p->gapOpen) && (expSubExtend == (double) p->gapExtend))
		return;

	smax = (double) scoring->sub[0][0];
	for (a=0 ; a<256 ; a++) for (b=0 ; b<256 ; b++)
		if ((double) scoring->sub[a][b] > smax) smax = (double) scoring->sub[a][b];

	if (expSub == NULL)
		expSub = (double*) malloc_or_die ("realign expSub", 256*256*sizeof(double));
	for (a=0 ; a<256 ; a++) for (b=0 ; b<256 ; b++)
		expSub[a*256+b] = exp (lambda * ((double) scoring->sub[a][b] - smax));
	meaGo = exp (-lambda * ((double) p->gapOpen + (double) p->gapExtend + smax/2));
	meaGe = exp (-lambda * ((double) p->gapExtend + smax/2));

	expSubFor    = scoring;
	expSubLambda = lambda;
	expSubOpen   = (double) p->gapOpen;
	expSubExtend = (double) p->gapExtend;
	}

//----------
//
// mea_posteriors--
//	带约束 |i-j| <= w 的前向-后向，结果 post[i*(2w+1) + (j+w-i)] = 格 (i,j)
//	（target 第 i 个碱基与 query 第 j 个碱基对齐，1 起）的后验；i 或 j 为 0 的格
//	记 0。带外的路径不计入（带外权重视为 0）。
//
// 带外的格不必写哨兵：带随 i 单调移动，行数组在开头整体清零，相邻行越出本行带
// 的那一格（前向的 hi+1、后向的 lo-1）此前从未被写过，读到的仍是 0。
//
// 返回 false 表示某行全部下溢、缩放溢出，或前向与后向的总权重相差超过 1e-6（相
// 对）——λ 很大时，与行内最大值相差 e^-709 以上的格会下溢成 0，丢掉的质量在两
// 个方向上不同。调用方保留 Viterbi 路径。正确性由 realign_test 的穷举对照保证。
//
//----------

#define meaIx(i,j) ((i)*width + ((j)+w-(i)))

static int mea_posteriors
   (const u8*	t,
	u64			lt,
	const u8*	q,
	u64			lq,
	u64			w)
	{
	const double go = meaGo, ge = meaGe;
	u64			width = 2*w + 1;
	u64			i, j, lo, hi;
	double		*mP, *mC, *xP, *xC, *yP, *yC, *swap;
	double		m, x, y, tot, mx, inv, mLeft, xLeft, yLeft;
	double		logZ, logB, logZb, r, diag, right, down, bm, bx, by, bxRight;
	const double* es;

	grow_double (&post, &postCap, (lt+1) * width, "realign post");
	grow_double (&rowScale, &rowScaleCap, lt+1, "realign rowScale");
	grow_mea_rows (lq + 2);
	mP = meaRow[0];  mC = meaRow[1];
	xP = meaRow[2];  xC = meaRow[3];
	yP = meaRow[4];  yC = meaRow[5];
	for (j=0 ; j<=lq+1 ; j++)
		mP[j] = mC[j] = xP[j] = xC[j] = yP[j] = yC[j] = 0;

	// 前向

	for (i=0 ; i<=lt ; i++)
		{
		lo = (i > w)? i - w : 0;
		hi = (i + w < lq)? i + w : lq;
		es = (i > 0)? expSub + 256*t[i-1] : NULL;
		mLeft = xLeft = yLeft = 0;
		mx = 0;
		for (j=lo ; j<=hi ; j++)
			{
			if ((i == 0) && (j == 0))
				{ m = 1;  x = y = 0; }
			else
				{
				m = ((i > 0) && (j > 0))? (mP[j-1] + xP[j-1] + yP[j-1]) * es[q[j-1]] : 0;
				x = (j > lo)? (mLeft + yLeft) * go + xLeft * ge : 0;
				y = (i > 0)?  (mP[j] + xP[j]) * go + yP[j] * ge : 0;
				}
			mC[j] = mLeft = m;
			xC[j] = xLeft = x;
			yC[j] = yLeft = y;
			tot = m + x + y;
			if (tot > mx) mx = tot;
			}
		if ((!(mx > 0)) || (!isfinite (mx))) return false;
		inv = 1 / mx;
		for (j=lo ; j<=hi ; j++)
			{
			mC[j] *= inv;  xC[j] *= inv;  yC[j] *= inv;
			post[meaIx(i,j)] = mC[j];
			}
		rowScale[i] = ((i > 0)? rowScale[i-1] : 0) + log (mx);
		swap = mP;  mP = mC;  mC = swap;
		swap = xP;  xP = xC;  xC = swap;
		swap = yP;  yP = yC;  yC = swap;
		}
	tot = mP[lq] + xP[lq] + yP[lq];
	if (!(tot > 0)) return false;
	logZ = log (tot) + rowScale[lt];

	// 后向（复用同一组行数组：mP/xP/yP 存第 i+1 行的 B_M/B_X/B_Y），同时算后验

	for (j=0 ; j<=lq+1 ; j++)
		mP[j] = mC[j] = xP[j] = xC[j] = yP[j] = yC[j] = 0;
	logB = 0;
	for (i=lt ; ; i--)
		{
		lo = (i > w)? i - w : 0;
		hi = (i + w < lq)? i + w : lq;
		es = (i < lt)? expSub + 256*t[i] : NULL;
		bxRight = 0;
		mx = 0;
		for (j=hi ; ; j--)
			{
			if ((i == lt) && (j == lq))
				bm = bx = by = 1;
			else
				{
				diag  = ((i < lt) && (j < lq))? es[q[j]] * mP[j+1] : 0;
				right = (j < hi)? bxRight : 0;
				down  = (i < lt)? yP[j] : 0;
				bm = diag + go*right + go*down;
				bx = diag + ge*right + go*down;
				by = diag + go*right + ge*down;
				}
			mC[j] = bm;  xC[j] = bx;  yC[j] = by;
			bxRight = bx;
			if (bm > mx) mx = bm;
			if (bx > mx) mx = bx;
			if (by > mx) mx = by;
			if (j == lo) break;
			}
		if ((!(mx > 0)) || (!isfinite (mx))) return false;
		inv = 1 / mx;
		for (j=lo ; j<=hi ; j++)
			{ mC[j] *= inv;  xC[j] *= inv;  yC[j] *= inv; }
		logB += log (mx);

		// 后验 = 前向 M × 后向 B_M / Z

		r = exp (rowScale[i] + logB - logZ);
		if (!isfinite (r)) return false;
		for (j=lo ; j<=hi ; j++)
			post[meaIx(i,j)] = ((i > 0) && (j > 0))? post[meaIx(i,j)] * mC[j] * r : 0;

		swap = mP;  mP = mC;  mC = swap;
		swap = xP;  xP = xC;  xC = swap;
		swap = yP;  yP = yC;  yC = swap;
		if (i == 0) break;
		}

	// 段首按 M 状态出发，后向总权重应与前向一致

	if (!(mP[0] > 0)) return false;
	logZb = log (mP[0]) + logB;
	if (fabs (logZb - logZ) > 1e-6) return false;
	return true;
	}

//----------
//
// mea_path_obj--
//	一条路径（逆序列操作）在 MEA 目标 Σ_匹配列 (post − θ) 下的值。路径须在带内。
//
//----------

static double mea_path_obj
   (const u8*	opRev,
	u64			n,
	u64			w,
	double		theta)
	{
	u64			width = 2*w + 1;
	u64			i = 0, j = 0, k;
	double		obj = 0;

	for (k=n ; k>0 ; )
		{
		switch (opRev[--k])
			{
			case colSub: i++;  j++;  obj += post[meaIx(i,j)] - theta;  break;
			case colIns: j++;                                          break;
			default:     i++;                                          break;
			}
		}
	return obj;
	}

//----------
//
// mea_decode--
//	带内最大化 Σ_匹配列 (post − θ)（gap 列记 0）的全局路径。结果（逆序列操作）
//	写入 newOp，返回列数；*obj 为目标值。
//
// 平局与 banded_global 相同：rightAlign 时优先删除、其次插入、最后对角线（indel
// 靠右）；否则优先对角线、其次插入、最后删除。带外哨兵同 mea_posteriors，靠开
// 头整体置为 meaNegInf。
//
//----------

static u64 mea_decode
   (u64			lt,
	u64			lq,
	u64			w,
	double		theta,
	int			right,
	double*		obj)
	{
	u64			width = 2*w + 1;
	u64			i, j, lo, hi, n;
	double		*vP, *vC, *swap;
	double		v, d, e, f, vLeft;
	u8			bits;

	grow_u8 (&tb, &tbCap, (lt+1) * width, "realign traceback");
	grow_mea_rows (lq + 2);
	vP = meaRow[0];  vC = meaRow[1];
	for (j=0 ; j<=lq+1 ; j++) vP[j] = vC[j] = meaNegInf;

	for (i=0 ; i<=lt ; i++)
		{
		lo = (i > w)? i - w : 0;
		hi = (i + w < lq)? i + w : lq;
		vLeft = meaNegInf;
		for (j=lo ; j<=hi ; j++)
			{
			if ((i == 0) && (j == 0))
				{ v = 0;  bits = colSub; }
			else
				{
				d = ((i > 0) && (j > 0))? vP[j-1] + post[meaIx(i,j)] - theta : meaNegInf;
				e = (j > lo)? vLeft : meaNegInf;
				f = (i > 0)?  vP[j] : meaNegInf;
				v = d;  bits = colSub;
				if ((right)? (e >= v) : (e > v)) { v = e;  bits = colIns; }
				if ((right)? (f >= v) : (f > v)) { v = f;  bits = colDel; }
				}
			vC[j] = vLeft = v;
			tb[meaIx(i,j)] = bits;
			}
		swap = vP;  vP = vC;  vC = swap;
		}
	*obj = vP[lq];

	grow_u8 (&newOp, &newOpCap, lt + lq + 1, "realign newOp");
	n = 0;
	i = lt;  j = lq;
	while ((i > 0) || (j > 0))
		{
		bits = tb[meaIx(i,j)];
		newOp[n++] = bits;
		if      (bits == colSub) { i--;  j--; }
		else if (bits == colIns) j--;
		else                     i--;
		}
	return n;
	}

//----------
//
// mea_chunk--
//	对一段做 MEA。进入时 newOp 中是带宽 w 下的 Viterbi 路径（nVit 列，逆序）；
//	返回时 newOp 中是 MEA 路径，返回列数。MEA 路径碰到带边时按 banded_global 的
//	规则加宽重算；矩阵过大或数值下溢时保留 Viterbi 路径。
//
// 校验：Viterbi 路径在带内，是可行解，所以 MEA 目标值不低于它（违反即退出）。
//
//----------

static u64 mea_chunk
   (const u8*			t,
	u64					lt,
	const u8*			q,
	u64					lq,
	u64					w,
	u64					wMax,
	u64					nVit,
	const realignparams* p)
	{
	u64			n = 0;
	double		obj, vitObj;

	grow_u8 (&vitOp, &vitOpCap, nVit + 1, "realign vitOp");
	memcpy (vitOp, newOp, nVit);
	realign_count_stat (meaChunks);

	for (;;)
		{
		if ((lt+1) * (2*w+1) > meaMaxCells) goto keep_viterbi;
		if (!mea_posteriors (t, lt, q, lq, w)) goto keep_viterbi;
		vitObj = mea_path_obj (vitOp, nVit, w, p->meaTheta);
		n = mea_decode (lt, lq, w, p->meaTheta, p->rightAlign, &obj);
		if (obj < vitObj - 1e-9 * (double) (lt + lq + 1))
			suicidef ("realign: MEA objective %.12g below Viterbi path's %.12g", obj, vitObj);
		if (w >= wMax) break;
		if (max_drift (newOp, n) < w) break;
		w = 2*w + 16;
		if (w > wMax) w = wMax;
		realign_count_stat (meaWidened);
		}

	if ((n != nVit) || (memcmp (newOp, vitOp, n) != 0))
		{ realign_count_stat (meaChanged); }
	return n;

keep_viterbi:
	memcpy (newOp, vitOp, nVit);
	realign_count_stat (meaFallback);
	return nVit;
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
			if (p->mea)
				n = mea_chunk (target->v + ti0, lt, query->v + qi0, lq, w, wMax, n, p);
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

	if (params->mea) mea_prepare (scoring, params);
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
	if (realignStats.meaChunks > 0)
		{
		fprintf (f, "         mea chunks: %s\n", ucommatize(realignStats.meaChunks));
		fprintf (f, "        mea widened: %s\n", ucommatize(realignStats.meaWidened));
		fprintf (f, "       mea fallback: %s\n", ucommatize(realignStats.meaFallback));
		fprintf (f, "        mea changed: %s\n", ucommatize(realignStats.meaChanged));
		}
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
	if (realignStats.meaChunks > 0)
		{
		(*func) (f, "realign_mea_chunks=%" PRIu64 "\n",   realignStats.meaChunks);
		(*func) (f, "realign_mea_widened=%" PRIu64 "\n",  realignStats.meaWidened);
		(*func) (f, "realign_mea_fallback=%" PRIu64 "\n", realignStats.meaFallback);
		(*func) (f, "realign_mea_changed=%" PRIu64 "\n",  realignStats.meaChanged);
		}
#endif // collect_stats
	}
