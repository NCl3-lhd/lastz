# 实施方案：lastz 内的框内全局重比对（`--realign`）

2026-11-05。依据：`docs/realign_ksw2_eval.md`。服务器原型（Python + ksw2）在 simTest 上把 G3a + 500/30 的 F_np 从 0.7463 提到 0.7609（重比对 400/8，indel 右对齐）。本方案把它做进 lastz，用自写的整数 DP，不依赖 ksw2。

## 1. 做什么

对每条已经确定、即将输出的 gapped 比对，在它的框内重新摆放 gap。框的首尾坐标不变，参与的碱基不变，只改 edit script。

1. **切段**（同原型）：沿原比对路径，在长度 ≥ 11 的无 gap 段中点处切开，相邻切点至少相隔 1000 列。切点在原路径上，所以每段两端固定。
2. **每段做带约束的全局比对**：Gotoh 单仿射，替换分直接用 lastz 的 `scoring->sub`（与 `score_alignment` 相同，含 `--ambiguous` 和小写的处理），gap 罚分用 `--realign` 给的 O、E。
   - 带约束：以主对角线为中心，带宽 w = 原路径在段内的最大对角线偏移 + 16。
   - 回溯出的路径碰到带边（|i−j| ≥ w）时，把 w 改为 2w+16 重算，直到不碰边或带已覆盖整段。自写 DP 在带内本身是精确的，加宽只是为了不让带约束改变结果。
   - 平局按 indel 右对齐（默认）或左对齐处理。
3. **把各段拼成新的 edit script**，替换原比对的 script，然后重算 `s`。
4. **校验**（违反即 `suicide`，属于实现错误）：
   - 每段新路径在重比对打分下的分 ≥ 原路径的分（原路径落在初始带内，是可行解）；
   - 新 script 恰好消耗原比对两侧的长度。

## 2. 决策项（括号内为推荐）

- **D1 适用范围**：
  - (a) 通用选项 `--realign`，任何 gapped 输出都能用，不要求 `--nlchain`。【推荐】重比对只依赖最终的比对列表，与 nlchain 无关。不过收益只在 nlchain G3a 配置下验证过，文档里会写明。
  - (b) 只在 `--nlchain` 下可用，命名为 `--nlchainrealign`。
- **D2 参数形式**：
  - `--realign=<O>,<E>[,right|left]`，O、E 必须显式给出，单位与 `--gap` 相同，平局方向默认 right。【推荐】最优值依赖打分矩阵，原型只在 HoxD55 上测过（400,8）；没有可靠的默认值，不如强制显式给出。
  - 或：单写 `--realign` 时默认取 O = 当前 gap open、E = 8。
- **D3 位置**：`finish_one_strand` 中"print the gapped alignments"之前，即重新串链（nlchain）或 `tweener_interpolate`（非 nlchain）之后、`mirror_alignments` 之前。
  - identity/coverage 等过滤仍按重比对之前的比对判断。
  - 动态掩蔽只用 beg1/end1（`masking.c` `census_mask_aligns`），不受影响。
- **D4 输出分数**：
  - (a) 按 lastz 自身的打分重算（`score_alignment`，gap 罚分取 `--gap`）。【推荐】这样输出里的分数在所有格式下含义一致。代价是重比对后的路径在 lastz 打分下可能更低，有的会低于 `--gappedthresh`，但不会因此再被过滤（过滤已经在前面做完）。
  - (b) 按重比对的打分。
- **D5 不做双仿射**：原型里三组双仿射都不如单仿射。
- **D6 切段与带宽的常数**（1000、11、16）写成源码常量，不开选项，与原型一致；原型没有扫过这几个值。
- **D7 互斥与约束**：
  - `--realign` 要求 gapped，与 `--nogapped` 一起用时报错；O、E 不能为负。
  - `lastz_D`（浮点打分）同样可用：DP 用 `score` 类型。
  - `--format=segments`（deGapify）照常，按新 script 输出。

## 3. 文件

| 文件 | 改动 |
|---|---|
| `src/realign.h` / `src/realign.c`（新） | `realignparams`（O、E、平局方向）；`alignel* realign_aligns (seq* target, seq* query, scoreset* scoring, alignel* alignList, const realignparams*)` 原地替换 script 和分数；切段、带约束 Gotoh（3 个状态，int 或 `score`，回溯用每格 1 字节）、加宽、拼接；`-Dcollect_stats` 计数（比对数、段数、加宽次数、分数提高的段数、列数变化） |
| `src/realign_test.c`（新） | 单元测试，`#include "realign.c"`（同 `nlrechain_test.c` 的做法） |
| `src/lastz.h` / `src/lastz.c` | 参数与默认值（关闭）、`--realign=` 解析、互斥检查、调用点、`debugClockRealign`、统计输出、`--help` 中加一行 |
| `src/Makefile` | `realign` 加入构建；新增 `make realign_test` |
| `.gitignore` | 加入 `src/realign_test` |

## 4. 步骤与验证

1. **realign.c 与单元测试**（本地 VM，gcc 11 `-Werror`）：
   - DP 对暴力参考解：随机小例子（长度 1–60，替换率 0.1/0.3/0.6，indel 率 0/0.05/0.2，含 N 与小写），每段的分 = 独立写的 O(nm) Gotoh（带约束与不带约束各一份）；回溯路径重算的分 = DP 分；路径恰好消耗两条序列。
   - 加宽：构造一定会碰带边的例子，确认加宽后结果等于不带约束的最优。
   - 平局方向：同聚物（AAAA 对 AAA 这类）以及串联重复中的 indel，right 时 gap 在最右、left 时在最左；两者分数相同。
   - 切段：随机 edit script 的比对，切点都落在无 gap 段内；拼回后两侧序列与 beg/end 不变；新 script 的分（重比对打分）≥ 原 script 的分。
   - 变异检验：把加宽条件、平局比较、gap 状态转移各改一处，测试应报出失败。
2. **接入 lastz**：
   - `make base_tests` 全过。
   - 不加 `--realign` 时，与 07e0f6f 的输出逐字节一致：沿用 §11.2 的 24 组非 nlchain 选项和 8 组 nlchain 选项，比较 stdout、stderr、退出码。
   - 互斥报错各一例；`lastz_D` 加 `--realign` 能跑通。
   - 冒烟：pseudocat/pseudopig、aglobin 上加 `--realign=400,8`，比对数和 beg/end 都不变，输出格式（maf、axt、general、sam、cigar、segments）都能正常解析。
3. **与原型对照**（服务器，simTest 20 个作业）：
   - G3a + `--gap=500,30` + `--realign=400,8`（right / left 各一组），用 run_tune/eval_tune 评测。
   - 预期与 `rk_r400_8` / `rk_g400_10` 的 F 相差 ≤ 0.002：差别只来自原型的 int8 取整（s=8 时转换分 −25 → −24）。
   - 段级对照：在单元测试里加一个读 TSV 的模式，读入原型用过的真实段，按同一套缩放参数算分，与 `ksw_realign` 的分逐段一致。这一项可以排除原型取整之外的实现差异。
   - CPU：用统计版二进制记录 `realign` 阶段的耗时，预期远低于原型的 9.5%（原型是 Python）。
4. **收尾**：
   - 在 `docs/realign_ksw2_eval.md` 里补"lastz 实现"一节；
   - 在 `nlchain_design.md` §8.3 的推荐配置里加上 `--realign=400,8`（同时写明 `--gap=500,30`）；
   - `--help` 中加 `--realign` 的说明。

## 5. 未包括

- 更细的 E（1–7、9–15）扫描、切段长度和带宽余量扫描、以 400/30 的 lastz 输出为输入。实现完成后，用 lastz 内的 `--realign` 直接扫，比原型方便（不再受 int8 粒度限制）。
- pair-HMM 后验解码（MEA）。
