//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlrechain.h
//
//----------

#ifndef nlrechain_H				// (prevent multiple inclusion)
#define nlrechain_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef nlrechain_owner
#define global
#else
#define global extern
#endif

//----------
//
// data structures and types
//
// 设计见 docs/nlchain_rechain_plan.md；本文件中的 §x 均指该文档章节。
//
//----------

// nlrgaptable--
//	axtChain linearGap 的三行代价表（§2 第 4 步 gapCalcCost）。表首位置为 1。

typedef struct nlrgaptable
	{
	u32			len;
	const unspos* position;		// 严格递增的 d
	const s64*	qGap;			// dt = 0 时按 dq 查
	const s64*	tGap;			// dq = 0 时按 dt 查
	const s64*	bothGap;		// dq、dt 都 > 0 时按 dq+dt 查
	} nlrgaptable;

extern const nlrgaptable nlrGapCostMedium;
extern const nlrgaptable nlrGapCostLoose;

// nlrechainparams--
//	gapped 之后重新串链的参数（§3）。

typedef struct nlrechainparams
	{
	const nlrgaptable* gapCost;	// 内置 medium / loose，跟随 --nlchaingapcost
	score		minChainScore;	// T_chain：重算后的链分低于它的链整条删除；
								// .. 0 表示关闭（调用方不调用 nl_rechain_aligns）
	FILE*		reportFile;		// 剥链结果（kent chain 格式，删除前），可为 NULL
	} nlrechainparams;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	u64 numAlignsIn;			// 输入的 gapped 比对数
	u64 numAlignsOut;			// 输出的比对数
	u64 numBlocks;				// 拆出的无 gap 块数
	u64 numChains;				// 剥出的链数
	u64 numChainsRemoved;		// 低于下限而删除的链数
	u64 bpIn;					// 输入碱基对（按块累加）
	u64 bpMerged;				// 完全重复块合并掉的碱基对
	u64 bpTrimmed;				// 链内重叠裁剪（含整块删除）掉的碱基对
	u64 bpRemoved;				// 删除的链上的碱基对
	u64 bpKept;					// 保留的碱基对
	} nlrechainStats;

#define nlrechain_count_stat(field)   ++nlrechainStats.field
#define nlrechain_add_stat(field,val) (nlrechainStats.field += val)
#else
#define nlrechain_count_stat(field)
#define nlrechain_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void nlrechain_zero_stats    (void);
void nlrechain_show_stats    (FILE* f);
void nlrechain_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in nlrechain.c
//
//----------

// alignList：一次 finish_one_strand 中过滤后的 gapped 比对（单一链方向；反链时
// query 已反向互补）。返回替换后的新列表，原列表被释放。
alignel* nl_rechain_aligns (seq* target, seq* query, scoreset* scoring,
                            alignel* alignList, const nlrechainparams* params);

#undef global
#endif // nlrechain_H
