//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: realign_test.c
//
//----------
//
// realign 的单元测试（docs/realign_plan.md §4 第 1 步）。
//
// 直接包含 realign.c 以访问其内部函数。
//	(1) banded_global 对独立写的 O(nm) Gotoh 参考解：分数相等；回溯路径重算的
//	    分等于 DP 分、恰好消耗两条序列、不出带。
//	(2) 加宽：原路径走对角线而最优路径偏移远大于初始带宽时，加宽后的结果等于
//	    不带约束的最优。
//	(3) 平局方向：同聚物与串联重复中的 indel，right 靠右、left 靠左，分数相同。
//	(4) realign_one：随机比对（含多段），首尾坐标与两侧长度不变、新 script 在重
//	    比对打分下的分不低于原 script、输出分数等于 score_alignment。
//
// 用法:
//	realign_test [轮数]
//	realign_test --tsv <mat25> <O> <E> <right>   （读 stdin 的段，格式同服务器
//	    原型 ksw_realign：id \t target \t query \t w；输出 id \t score \t 最终带宽；
//	    mat25 为 C A G T N 次序的 5x5 矩阵，用于与原型逐段对照分数）
//
//----------

#define collect_stats			// (统计计数用于确认确实切出了多段)
#include "realign.c"
#include "identity_dist.h"
#include "coverage_dist.h"
#include "continuity_dist.h"
#include "output.h"

//----------
//
// 链接用的桩函数：gapped_extend.o 中 gappily_extend_hsps 引用了它们，本测试
// 不会调到（同 nlrechain_test.c）
//
//----------

static void not_in_test (const char* name)
	{ suicidef ("realign_test: %s should not be called", name); }

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
		if (numFailures <= 30)                                               \
			{                                                                \
			fprintf (stderr, "FAIL %s:%d: ", __FILE__, __LINE__);            \
			fprintf (stderr, __VA_ARGS__);                                   \
			fprintf (stderr, "\n");                                          \
			}                                                                \
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

static double rnd01 (void)
	{ return rnd (1000000) / 1000000.0; }

// HOXD55（UCSC 远缘用）与 lastz 默认的 HOXD70

static scoreset* hoxd55 (void)
	{
	score t[4][4] = { {   91,  -90,  -25, -100 },
	                  {  -90,  100, -100,  -25 },
	                  {  -25, -100,  100,  -90 },
	                  { -100,  -25,  -90,   91 } };
	return new_dna_score_set (t, -100, -100, 400, 30);
	}

static scoreset* hoxd70 (void)
	{
	score t[4][4] = { {   91, -114,  -31, -123 },
	                  { -114,  100, -125,  -31 },
	                  {  -31, -125,  100, -114 },
	                  { -123,  -31, -114,   91 } };
	return new_dna_score_set (t, -100, -100, 400, 30);
	}

static const char randomChars[] = "ACGTACGTACGTACGTacgtN";

static void random_bases (u8* s, u32 n, u32 alphabet)
	{
	u32 k;
	for (k=0 ; k<n ; k++)
		{
		if (alphabet < 4) s[k] = "ACGT"[rnd(alphabet)];
		             else s[k] = randomChars[rnd(sizeof(randomChars)-1)];
		}
	}

//----------
//
// ref_global--
//	独立的参考实现：完整的 (lt+1)×(lq+1) 三个矩阵，带约束 |i-j| <= w（w 足够大
//	即不带约束），返回全局最优分。
//
//----------

static s64 ref_global
   (scoreset*	scoring,
	s64			gapOpen,
	s64			gapExtend,
	const u8*	t,
	u32			lt,
	const u8*	q,
	u32			lq,
	u64			w)
	{
	const s64	NEG = -((s64) 1 << 50);
	u32			n = (lq+1), i, j;
	s64*		H = (s64*) malloc_or_die ("ref H", (lt+1) * n * sizeof(s64));
	s64*		E = (s64*) malloc_or_die ("ref E", (lt+1) * n * sizeof(s64));
	s64*		F = (s64*) malloc_or_die ("ref F", (lt+1) * n * sizeof(s64));
	s64			best, v;

	for (i=0 ; i<=lt ; i++)
		for (j=0 ; j<=lq ; j++)
			{
			s64 h = NEG, e = NEG, f = NEG;
			u64 off = (i > j)? i - j : j - i;
			if (off <= w)
				{
				if (j > 0)
					{
					e = H[i*n + j-1] - gapOpen - gapExtend;
					v = E[i*n + j-1] - gapExtend;
					if (v > e) e = v;
					}
				if (i > 0)
					{
					f = H[(i-1)*n + j] - gapOpen - gapExtend;
					v = F[(i-1)*n + j] - gapExtend;
					if (v > f) f = v;
					}
				if ((i == 0) && (j == 0)) h = 0;
				else
					{
					if ((i > 0) && (j > 0))
						h = H[(i-1)*n + j-1] + (s64) scoring->sub[t[i-1]][q[j-1]];
					if (e > h) h = e;
					if (f > h) h = f;
					}
				}
			if (h < NEG) h = NEG;
			if (e < NEG) e = NEG;
			if (f < NEG) f = NEG;
			H[i*n + j] = h;  E[i*n + j] = e;  F[i*n + j] = f;
			}
	best = H[lt*n + lq];
	free (H);  free (E);  free (F);
	return best;
	}

// 逆序列操作 → 正序

static void reverse_ops (const u8* rev, u64 n, u8* out)
	{
	u64 k;
	for (k=0 ; k<n ; k++) out[k] = rev[n-1-k];
	}

// 列操作消耗的两侧长度

static void op_lengths (const u8* op, u64 n, u64* lt, u64* lq)
	{
	u64 k;
	*lt = *lq = 0;
	for (k=0 ; k<n ; k++)
		{
		if (op[k] != colIns) (*lt)++;
		if (op[k] != colDel) (*lq)++;
		}
	}

static void op_string (const u8* op, u64 n, char* buf, u32 bufLen)
	{
	u64 k;
	u32 m = 0;
	for (k=0 ; (k<n) && (m+1<bufLen) ; k++)
		buf[m++] = "SID"[op[k]];
	buf[m] = 0;
	}

//----------
//
// (1) banded_global 对参考解
//
//----------

static void test_banded_vs_ref (u32 rounds)
	{
	scoreset*	sets[2];
	score		gaps[4][2] = { { 400, 30 }, { 400, 8 }, { 150, 10 }, { 500, 0 } };
	u8			t[200], q[200], fwd[400];
	u32			r, lt, lq, alphabet;
	u64			w, n, ct, cq, drift;
	score		best;
	s64			ref;
	realignparams p;

	sets[0] = hoxd55 ();  sets[1] = hoxd70 ();
	for (r=0 ; r<rounds ; r++)
		{
		scoreset* sc = sets[rnd(2)];
		u32 g = rnd(4);
		p.enabled    = true;
		p.gapOpen    = gaps[g][0];
		p.gapExtend  = gaps[g][1];
		p.rightAlign = rnd(2);
		alphabet = 1 + rnd(5);			// 1..4 只用前几个碱基（制造平局），5 含小写与 N
		lt = rnd(60);  lq = rnd(60);
		random_bases (t, lt, alphabet);
		if ((lq > 0) && (lt > 0) && (rnd(2)))
			{									// q 由 t 变异得到
			u32 k, m = 0;
			for (k=0 ; (k<lt) && (m<lq) ; k++)
				{
				double x = rnd01();
				if (x < 0.08) continue;
				if (x < 0.16) { q[m++] = "ACGT"[rnd(4)];  if (m < lq) q[m++] = t[k];  continue; }
				q[m++] = (rnd01() < 0.25)? "ACGT"[rnd(4)] : t[k];
				}
			lq = m;
			}
		else
			random_bases (q, lq, alphabet);
		w = ((lt > lq)? lt - lq : lq - lt) + rnd(20);
		if (rnd(4) == 0) w = 1000;

		n = banded_global (sc, &p, t, lt, q, lq, w, &best);
		ref = ref_global (sc, p.gapOpen, p.gapExtend, t, lt, q, lq, w);
		check ((s64) best == ref, "banded_global score " scoreFmtSimple " vs ref %" PRId64
		       " (lt=%u lq=%u w=%" PRIu64 " O=" scoreFmtSimple " E=" scoreFmtSimple " right=%d)",
		       best, ref, lt, lq, w, p.gapOpen, p.gapExtend, p.rightAlign);
		reverse_ops (newOp, n, fwd);
		op_lengths (fwd, n, &ct, &cq);
		check ((ct == lt) && (cq == lq), "path consumes %" PRIu64 ",%" PRIu64 " of %u,%u", ct, cq, lt, lq);
		check (path_score (sc, &p, t, q, fwd, n) == best, "path rescore != DP score");
		drift = max_drift (fwd, n);
		check (drift <= w, "path leaves band (drift %" PRIu64 " > w %" PRIu64 ")", drift, w);
		}
	free_score_set ("hoxd55", sets[0]);
	free_score_set ("hoxd70", sets[1]);
	}

//----------
//
// 测 realign_one 用的比对构造
//
//----------

typedef struct testaln
	{
	seq*		target;
	seq*		query;
	alignel*	a;
	} testaln;

static seq* make_seq_from (const u8* s, u32 len)
	{
	seq* sq = (seq*) zalloc_or_die ("test seq", sizeof(seq));
	sq->v = (u8*) malloc_or_die ("test seq v", len+1);
	memcpy (sq->v, s, len);
	sq->v[len] = 0;
	sq->len = sq->trueLen = len;
	sq->startLoc = 1;
	return sq;
	}

static void free_test_seq (seq* sq)
	{
	free_if_valid ("test seq v", sq->v);
	free_if_valid ("test seq",   sq);
	}

// make_aln--
//	随机构造一条比对：target 核心长 coreLen，两侧各有 flank 个无关碱基；query
//	由 target 核心按 sub/indel 率变异得到，变异过程同时给出原 edit script。

static testaln make_aln (u32 coreLen, u32 flank, double subRate, double indelRate, u32 alphabet)
	{
	testaln		x;
	u8*			t = (u8*) malloc_or_die ("t", coreLen + 2*flank + 1);
	u8*			q = (u8*) malloc_or_die ("q", 10*coreLen + 2*flank + 16);	// 每个 target 碱基最多带 8 个插入
	u32			pos, m, d, ins;			// pos：核心内已消耗的 target 碱基数
	editscript*	script = edit_script_new ();

	if (coreLen < 2) coreLen = 2;
	random_bases (t, coreLen + 2*flank, alphabet);
	random_bases (q, flank, alphabet);
	m = flank;
	edit_script_sub (&script, 1);  q[m++] = t[flank];  pos = 1;	// 以替换开头
	while (pos + 1 < coreLen)
		{
		double x = rnd01();
		if (x < indelRate/2)
			{
			d = 1 + rnd(6);
			if (pos + d + 1 < coreLen) { edit_script_del (&script, d);  pos += d; }
			continue;
			}
		if (x < indelRate)
			{
			ins = 1 + rnd(8);
			edit_script_ins (&script, ins);
			random_bases (q + m, ins, alphabet);  m += ins;
			}
		edit_script_sub (&script, 1);
		q[m++] = (rnd01() < subRate)? "ACGT"[rnd(4)] : t[flank+pos];
		pos++;
		}
	edit_script_sub (&script, 1);  q[m++] = t[flank+pos];  pos++;	// 以替换结尾
	random_bases (q + m, flank, alphabet);
	x.target = make_seq_from (t, coreLen + 2*flank);
	x.query  = make_seq_from (q, m + flank);
	x.a = (alignel*) zalloc_or_die ("test alignel", sizeof(alignel));
	x.a->beg1 = flank + 1;  x.a->end1 = flank + pos;
	x.a->beg2 = flank + 1;  x.a->end2 = m;
	x.a->script = script;
	free (t);  free (q);
	return x;
	}

static void free_aln (testaln* x)
	{
	free_align_list (x->a);
	free_test_seq (x->target);
	free_test_seq (x->query);
	}

// script 在重比对打分下的分（展开后用 path_score）

static score script_realign_score (scoreset* sc, const realignparams* p, testaln* x)
	{
	editscript*	s = x->a->script;
	u64			n = 0, k;
	u32			opIx;
	u8*			op;
	score		v;

	for (opIx=0 ; opIx<s->len ; opIx++) n += edit_op_repeat (s->op[opIx]);
	op = (u8*) malloc_or_die ("ops", n + 1);
	n = 0;
	for (opIx=0 ; opIx<s->len ; opIx++)
		{
		u32 rpt = edit_op_repeat (s->op[opIx]);
		u8 code = (edit_op_operation(s->op[opIx]) == editopIns)? colIns
		        : (edit_op_operation(s->op[opIx]) == editopDel)? colDel : colSub;
		for (k=0 ; k<rpt ; k++) op[n++] = code;
		}
	v = path_score (sc, p, x->target->v + x->a->beg1 - 1, x->query->v + x->a->beg2 - 1, op, n);
	free (op);
	return v;
	}

static void script_lengths (editscript* s, u64* ct, u64* cq)
	{
	u32 opIx;
	*ct = *cq = 0;
	for (opIx=0 ; opIx<s->len ; opIx++)
		{
		u32 rpt = edit_op_repeat (s->op[opIx]);
		switch (edit_op_operation (s->op[opIx]))
			{
			case editopIns: *cq += rpt;  break;
			case editopDel: *ct += rpt;  break;
			default:        *ct += rpt;  *cq += rpt;  break;
			}
		}
	}

//----------
//
// (2) 加宽：q = t[40:] + 40 个随机碱基，原比对是一条对角线（全替换）。最优路
//	径要先删 40 再插 40，偏移 40 远大于初始带宽 16，必须加宽。
//
//----------

static void test_widening (void)
	{
	scoreset*		sc = hoxd70 ();
	realignparams	p = { true, 400, 8, true };
	u8				t[300], q[300];
	u32				k, len = 300;
	testaln			x;
	editscript*		script = edit_script_new ();
	u64				w0 = realignStats.numWidened;
	s64				ref;

	random_bases (t, len, 4);
	for (k=0 ; k<len-40 ; k++) q[k] = t[k+40];
	random_bases (q + len - 40, 40, 4);
	edit_script_sub (&script, len);
	x.target = make_seq_from (t, len);
	x.query  = make_seq_from (q, len);
	x.a = (alignel*) zalloc_or_die ("test alignel", sizeof(alignel));
	x.a->beg1 = x.a->beg2 = 1;
	x.a->end1 = x.a->end2 = len;
	x.a->script = script;

	ref = ref_global (sc, p.gapOpen, p.gapExtend, t, len, q, len, 100000);
	realign_one (x.target, x.query, sc, x.a, &p);
	check (realignStats.numWidened > w0, "widening was not triggered");
	check ((s64) script_realign_score (sc, &p, &x) == ref,
	       "widened result " scoreFmtSimple " != unbanded optimum %" PRId64,
	       script_realign_score (sc, &p, &x), ref);
	free_aln (&x);
	free_score_set ("hoxd70", sc);
	}

//----------
//
// (3) 平局方向
//
//----------

static void tie_case (scoreset* sc, const char* ts, const char* qs, int right, const char* expect)
	{
	realignparams	p = { true, 400, 30, right };
	u8				fwd[64];
	char			buf[64];
	score			best;
	u64				n;

	n = banded_global (sc, &p, (const u8*) ts, strlen(ts), (const u8*) qs, strlen(qs), 20, &best);
	reverse_ops (newOp, n, fwd);
	op_string (fwd, n, buf, sizeof(buf));
	check (strcmp (buf, expect) == 0, "tie %s vs %s right=%d: got %s, expected %s", ts, qs, right, buf, expect);
	}

static void test_ties (void)
	{
	scoreset* sc = hoxd70 ();
	tie_case (sc, "AAAA",     "AAA",      true,  "SSSD");
	tie_case (sc, "AAAA",     "AAA",      false, "DSSS");
	tie_case (sc, "AAA",      "AAAA",     true,  "SSSI");
	tie_case (sc, "AAA",      "AAAA",     false, "ISSS");
	tie_case (sc, "GACACACT", "GACACT",   true,  "SSSSSDDS");
	tie_case (sc, "GACACACT", "GACACT",   false, "SDDSSSSS");
	tie_case (sc, "GACACT",   "GACACACT", true,  "SSSSSIIS");
	tie_case (sc, "GACACT",   "GACACACT", false, "SIISSSSS");
	free_score_set ("hoxd70", sc);
	}

//----------
//
// (4) realign_one / realign_aligns 的整体性质
//
//----------

static void test_realign_one (u32 rounds)
	{
	scoreset*		sets[2];
	score			gaps[3][2] = { { 400, 8 }, { 400, 30 }, { 200, 10 } };
	u32				r;
	u64				chunks0 = realignStats.numChunks, aligns0 = realignStats.numAligns;
	u64				multi = 0;

	sets[0] = hoxd55 ();  sets[1] = hoxd70 ();
	for (r=0 ; r<rounds ; r++)
		{
		scoreset*		sc = sets[rnd(2)];
		u32				g = rnd(3);
		realignparams	p;
		u32				coreLen = (rnd(4) == 0)? 2000 + rnd(6000) : 20 + rnd(400);
		testaln			x = make_aln (coreLen, rnd(30), 0.05 + 0.4*rnd01(), 0.25*rnd01(), 1 + rnd(5));
		unspos			b1 = x.a->beg1, e1 = x.a->end1, b2 = x.a->beg2, e2 = x.a->end2;
		score			before, after;
		u64				ct, cq, c0 = realignStats.numChunks;

		p.enabled = true;  p.gapOpen = gaps[g][0];  p.gapExtend = gaps[g][1];  p.rightAlign = rnd(2);
		before = script_realign_score (sc, &p, &x);
		realign_aligns (x.target, x.query, sc, x.a, &p);
		after = script_realign_score (sc, &p, &x);
		script_lengths (x.a->script, &ct, &cq);
		check ((x.a->beg1 == b1) && (x.a->end1 == e1) && (x.a->beg2 == b2) && (x.a->end2 == e2),
		       "realign changed alignment bounds");
		check ((ct == e1 - b1 + 1) && (cq == e2 - b2 + 1), "new script lengths %" PRIu64 ",%" PRIu64
		       " vs %u,%u", ct, cq, e1 - b1 + 1, e2 - b2 + 1);
		check (after >= before, "realigned score " scoreFmtSimple " < original " scoreFmtSimple, after, before);
		check (x.a->s == score_alignment (sc, x.target->v, x.a->beg1-1, x.query->v, x.a->beg2-1, x.a->script),
		       "alignment score is not score_alignment of the new script");
		if (realignStats.numChunks - c0 > 1) multi++;
		free_aln (&x);
		}
	check (realignStats.numAligns - aligns0 == rounds, "numAligns stat");
	check (multi > 0, "no alignment was split into more than one chunk");
	fprintf (stderr, "realign_one: %u alignments, %" PRIu64 " chunks, %" PRIu64 " split into >1 chunk\n",
	         rounds, realignStats.numChunks - chunks0, multi);
	free_score_set ("hoxd55", sets[0]);
	free_score_set ("hoxd70", sets[1]);
	}

//----------
//
// --tsv 模式：与服务器原型 ksw_realign 逐段对照分数
//
//----------

static int tsv_mode (int argc, char** argv)
	{
	score			mat[25], tmpl[4][4];
	const char*		order = "CAGT";
	const char*		acgt  = "ACGT";
	realignparams	p;
	scoreset*		sc;
	char*			line = NULL;
	size_t			cap = 0;
	ssize_t			len;
	int				i, a, b;
	char*			s;

	if (argc != 6) suicide ("usage: realign_test --tsv <mat25> <O> <E> <right>");
	s = argv[2];
	for (i=0 ; i<25 ; i++)
		{
		mat[i] = (score) strtol (s, &s, 10);
		if (*s == ',') s++;
		}
	for (a=0 ; a<4 ; a++)
		for (b=0 ; b<4 ; b++)
			tmpl[a][b] = mat[(strchr(order, acgt[a]) - order) * 5 + (strchr(order, acgt[b]) - order)];
	sc = new_dna_score_set (tmpl, mat[24], mat[24], 0, 0);
	p.enabled = true;
	p.gapOpen = atoi (argv[3]);  p.gapExtend = atoi (argv[4]);  p.rightAlign = atoi (argv[5]);

	while ((len = getline (&line, &cap, stdin)) > 0)
		{
		char *id = line, *t, *q, *ws;
		u64  lt, lq, w, wMax, n, drift;
		score best;
		if (line[len-1] == '\n') line[--len] = 0;
		t  = strchr (id, '\t');  *t++  = 0;
		q  = strchr (t,  '\t');  *q++  = 0;
		ws = strchr (q,  '\t');  *ws++ = 0;
		lt = strlen (t);  lq = strlen (q);  w = strtoull (ws, NULL, 10);
		wMax = (lt > lq)? lt : lq;
		if (w > wMax) w = wMax;
		for (;;)					// 与 realign_one 相同的加宽规则
			{
			n = banded_global (sc, &p, (u8*) t, lt, (u8*) q, lq, w, &best);
			if (w >= wMax) break;
			drift = max_drift (newOp, n);
			if (drift < w) break;
			w = 2*w + 16;
			if (w > wMax) w = wMax;
			}
		printf ("%s\t" scoreFmtSimple "\t%" PRIu64 "\n", id, best, w);
		}
	free (line);
	free_score_set ("tsv scores", sc);
	return 0;
	}

//----------
//
// main
//
//----------

int main (int argc, char** argv)
	{
	u32 rounds = 3000;

	if ((argc > 1) && (strcmp (argv[1], "--tsv") == 0))
		return tsv_mode (argc, argv);
	if (argc > 1) rounds = (u32) atoi (argv[1]);

	test_banded_vs_ref (rounds);
	test_widening ();
	test_ties ();
	test_realign_one (rounds / 10 + 1);

	fprintf (stderr, "realign_test: %d checks, %d failures\n", numChecks, numFailures);
	return (numFailures == 0)? 0 : 1;
	}
