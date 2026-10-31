//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: nlchain.h
//
//----------

#ifndef nlchain_H				// (prevent multiple inclusion)
#define nlchain_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff

// establish ownership of global variables

#ifdef nlchain_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef nlchain_owner
int nlchain_dbgChaining = false;
#else
global int nlchain_dbgChaining;
#endif

//----------
//
// data structures and types
//
// 设计见 docs/nlchain_design.md；本文件中的 §x.y 均指该文档章节。
//
//----------

// nlgaptable--
//	正常边 g(d) 的分段线性表（axtChain linearGap 的 qGap 行，§5.1）。

typedef struct nlgaptable
	{
	u32			len;
	const unspos* position;		// 严格递增的 d
	const double* cost;			// 对应的 g(d)
	} nlgaptable;

// nlchainparams--
//	非线性串联的全部参数（§8.1、§8.3）。

typedef struct nlchainparams
	{
	const nlgaptable* gapCost;	// 正常边 g(d)：内置 medium / loose
	double		antiPen;			// R：两侧都有间隔时的每碱基罚分（同 --chain 的 anti）
	unspos		maxGap;				// 正常边：dq、dt 上限
	unspos		maxDiag;			// 正常边：|dq-dt| 上限
	double		overlapFrac;		// f_ov：相邻 HSP 在 query 上的最大重叠比例
	score		jumpPen;			// 异常边基础罚分 J
	score		spanCoef;			// λ
	unspos		spanScale;			// D0
	double		spanCap;			// L（log2 单位）
	score		soloScore;			// T_solo
	score		minSourceSupport;	// T_src，0 关闭
	u32			maxJumpScan;		// H_max
	score		secondaryScore;		// T_sec；0 关闭次级链
	u32			secondaryPerRegion;	// N_sec
	double		secondaryRatio;		// r_sec
	} nlchainparams;

// nlsegment--
//	输出的一个共线段（§7 第 3、4 步）。

typedef struct nlsegment
	{
	u32			first, count;		// 在该链方向压缩后的 segtable 中的下标区间
	u32			part1;				// target 分区下标
	u32			part2;				// query 分区下标（query 无分区时为 0）
	int			strand;				// +1 / -1
	unspos		qBeg, qEnd;			// 正链 query 坐标，半开
	unspos		tBeg, tEnd;			// target 全局坐标，半开
	double		support;			// 段分数（§7 第 3 步）
	char		jumpIn;				// 'S','I','T','D','G'
	int			isSecondary;		// 0 主路径，1 次级链
	u32			chainId;			// 主路径为 0，次级链从 1 编号
	} nlsegment;

typedef struct nlchainresult
	{
	u32			len, size;
	nlsegment*	seg;
	} nlchainresult;

// 内置的 g(d) 表（§5.1）

extern const nlgaptable nlGapCostMedium;
extern const nlgaptable nlGapCostLoose;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	u64 numHsps;				// 进入 DP 的 HSP 数（两条链合计）
	u64 numPrimaryHsps;			// 主路径 HSP 数
	u64 numPrimarySegments;		// 主路径段数
	u64 numJumpIn[5];			// 主路径各段的进入方式：S,I,T,D,G
	u64 numSecondaryChains;		// 保留的次级链数
	u64 numSecondaryHsps;		// 次级链 HSP 数
	u64 numJumpTruncations;		// 跳转扫描达到 H_max 被截断的次数
	} nlchainStats;

#define nlchain_count_stat(field)   ++nlchainStats.field
#define nlchain_add_stat(field,val) (nlchainStats.field += val)
#else
#define nlchain_count_stat(field)
#define nlchain_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void nlchain_zero_stats    (void);
void nlchain_show_stats    (FILE* f);
void nlchain_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in nlchain.c
//
//----------

// stMinus / stPlus：反链 / 正链的单链 HSP 表，各自已 merge、已带分数；反链 pos2 为反向互补坐标。
// 返回主路径分（未缩放）；两张表各自原地压缩为主路径和次级链的 HSP，out 按 (chainId, qBeg) 列出各段。
double nl_reduce_to_chains (seq* target, seq* query,
                            segtable* stMinus, segtable* stPlus,
                            const nlchainparams* params, nlchainresult* out);
void   nl_free_result      (nlchainresult* out);
void   nl_write_report_header (FILE* f);
void   nl_write_report     (FILE* f, seq* target, seq* query,
                            const nlchainresult* r);

#undef global
#endif // nlchain_H
