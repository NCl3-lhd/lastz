//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: realign.h
//
//----------

#ifndef realign_H				// (prevent multiple inclusion)
#define realign_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef realign_owner
#define global
#else
#define global extern
#endif

//----------
//
// data structures and types
//
// 设计见 docs/realign_plan.md；本文件中的 §x 均指该文档章节。
//
//----------

// realignparams--
//	输出前框内全局重比对的参数（§2 D2）。

typedef struct realignparams
	{
	int			enabled;		// true => 输出前做重比对
	score		gapOpen;		// 重比对用的 gap 罚分；长 l 的 gap 罚
	score		gapExtend;		// .. gapOpen + l*gapExtend
	int			rightAlign;		// true => 平局时 indel 靠右；false => 靠左
	} realignparams;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	u64 numAligns;				// 重比对的比对数
	u64 numChunks;				// 切出的段数
	u64 numWidened;				// 因路径碰到带边而加宽重算的次数
	u64 numImproved;			// 重比对分高于原路径的段数
	u64 colsIn;					// 重比对前的比对列数
	u64 colsOut;				// 重比对后的比对列数
	} realignStats;

#define realign_count_stat(field)   ++realignStats.field
#define realign_add_stat(field,val) (realignStats.field += val)
#else
#define realign_count_stat(field)
#define realign_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void realign_zero_stats    (void);
void realign_show_stats    (FILE* f);
void realign_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in realign.c
//
//----------

// alignList：一次 finish_one_strand 中即将输出的 gapped 比对（反链时 query 已反
// 向互补）。原地替换每条比对的 edit script 并按 scoring（lastz 自身的打分）重算
// 分数；比对的首尾坐标、次序不变。
alignel* realign_aligns (seq* target, seq* query, scoreset* scoring,
                         alignel* alignList, const realignparams* params);

#undef global
#endif // realign_H
