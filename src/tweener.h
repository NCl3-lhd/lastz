//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: tweener.h
//
//----------

#ifndef tweener_H				// (prevent multiple inclusion)
#define tweener_H

// other files

#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed matching stuff
#include "chain.h"				// segment chaining stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "edit_script.h"		// alignment edit script stuff
#include "nlchain.h"			// non-linear chaining stuff

// establish ownership of global variables

#ifdef tweener_owner
#define global
#else
#define global extern
#endif

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int   numTweeners;
	u64   totalArea;
	int   totalCoverage;
	int   tweenerCoverage;
	u64   numNlWindows;			// nlchain 插值：实际搜索的窗口数
	u64   numNlWindowsSkipped;	// nlchain 插值：被无损预筛跳过的窗口数
	u64   nlWindowArea;			// nlchain 插值：搜索窗口的面积和
	u64   numNlHsps[2];			// nlchain 插值：找到的 HSP 数，0 反链，1 正链
	u64   numPass1Dropped;		// 第一次主路径 HSP 在第二次串链中被丢掉的数目
	} tweenerStats;

// stats macros

#define tweener_count_stat(field)   ++tweenerStats.field
#define tweener_uncount_stat(field) --tweenerStats.field
#define tweener_set_stat(field,val) (tweenerStats.field = val)
#define tweener_add_stat(field,val) (tweenerStats.field += val)
#else
#define tweener_count_stat(field)
#define tweener_uncount_stat(field)
#define tweener_set_stat(field,val)
#define tweener_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void tweener_zero_stats    (void);
void tweener_show_stats    (FILE* f);
void tweener_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in tweener.c
//
//----------

alignel* tweener_interpolate (alignel* a, seq* seq1, seq* seq2,
                              int selfCompare, int inhibitTrivial,
                              const s8 charToBits[], seed* tweenSeed,
                              scoreset* scoring, scoreset* maskedScoring,
                              tback* tb, score xDrop, int gappedAllBounds,
                              score yDrop, int trimToPeak, score scoreThresh,
                              score diagPen, score antiPen,
                              int scale, chainer connect, u32 windowSize);
void     tweener_collect_nl_hsps (seq* target, seq* query, const u8* qToComplement,
                              const s8 upperCharToBits[], seed* innerSeed,
                              scoreset* maskedScoring, score xDrop,
                              score innerThreshold, u32 windowSize,
                              const nlchainresult* r,
                              const segtable* stMinus, const segtable* stPlus,
                              segtable** innerMinus, segtable** innerPlus);

#undef global
#endif // tweener_H
