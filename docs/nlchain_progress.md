# nlchain 实现进展

日期：2026-10-31。对应设计：`docs/nlchain_design.md`（下文 §x.y 均指该文档）。本次完成 §9 的第 1–3 步（本地），第 0、4、5 步（远端）未做。代码尚未提交。

## 1. 已完成

| 文件 | 内容 |
|---|---|
| `src/nlchain.h` / `src/nlchain.c`（新） | §8.1 接口；公式层、转移层 E1–E6、搜索层（每个 (p,σ) 一棵 K-d 树，用出口坐标做 box 查询；跳转池为按 R 降序的堆，配合待释放堆和同批冻结；H_max 安全阀）；主路径回溯与切段（S/I/T/D/G）；次级链（§5.4 第 1–6 步）；两张单链表原地压缩；段表 TSV；统计 |
| `src/nlchain_test.c`（新） | 单元测试，直接 `#include "nlchain.c"`，参考解只替换搜索层 |
| `src/lastz.h` / `src/lastz.c` | §8.2 第 1–8 条：参数与默认值、CLI、互斥检查、`collectHspsFromBoth`、调用点 `nlchain_both_strands()`（在反链 `split_anchors` 之后；两张表分别 merge 后再调 `nl_reduce_to_chains`）、统计、`debugClockNlChaining`、报告文件 |
| `src/Makefile` | `nlchain` 加入构建；新增 `make nlchain_test` |
| `.gitignore` | 加入 `src/nlchain_test` |

### CLI

```
--nlchain
--nlchainjump=<J>[,<λ>,<D0>,<L>]          默认 J=hspthresh, 500, 1000, 8
--nlchaingap=<maxGap>,<maxDiag>[,<R>]     默认 50000, 10000, 0
--nlchaingapcost=medium|loose             默认 medium
--nlchainsolo=<T_solo>                    默认 2×hspthresh
--nlchainsource=<T_src>                   默认 0（关闭）
--nlchainoverlap=<f_ov>                   默认 0.5
--nlchainscan=<H_max>                     默认 1000
--nlchainsecondary=<T_sec>[,<N_sec>,<r_sec>]  默认 2×hspthresh, 5, 0.5；T_sec=0 关闭次级链
--nlchainreport=<file>                    段表 TSV
```

自适应 hspthresh（`--hspthresh=top...`）下，J、T_solo、T_sec 中的 hspthresh 取本次两条链收集到的最低 HSP 分（在 merge 之前统计）。

## 2. 验证结果

- **单元测试** `make nlchain_test`：默认 3000 轮，约 52 万项检查，0 失败（跑 6000 轮约 96 万项，同样 0 失败）。
  - 公式层：`gap_cost` 在 medium/loose 两张表的 11 个节点上都等于表值；验证了插值和外推；手算了 C_N、gain、P_J、D 的三种情形、jump 类型；两端同时重叠时只撤回一次（ov ≥ ov_t、ov_t > ov、dq=dt=−k 三种情况）；Normal(i) 的 maxGap、maxDiag、ov_max、target 推进四处边界。
  - 搜索层对暴力枚举：n=12 / 300 / 5000 三档随机数据，U/R/L 的分数、前驱、段内 HSP 数、前驱状态、边类型全部一致；H_max 未触发。
  - 对穷举全部合法路径：n ≤ 12，每个 HSP 的 U/R/L 最优路径与 DP 结果逐元组一致（平局按 §5.3 规则）。随机数据覆盖到了全部转移类型：R←N、U←J、R←J、L←N，以及 T_src 挡掉跳出。
  - 变异检验：把两处剪枝条件改回 `≤`，测试分别报出 4919 项和 8 项失败；正式实现为 0 失败。说明平局用例确实会走到剪枝的边界。
  - 次级链：较弱的旁系拷贝作为一整条链（6 个 HSP）被回收；长主路径（40 个 HSP）内部的 3-HSP 片段重复拷贝被 r_sec 保留；依次出现倒位和易位时，段的 jumpIn 为 S/I/T。
- **不加 `--nlchain` 的回归**：`make base_tests` 全部通过。另外 22 组选项（lav/maf/axt/general/sam/cigar/segments、`--chain`、`--chain=20,1`、`--nogapped`、自适应阈值、`--queryhspbest`、`--queryhsplimit`、`--strand=plus/minus`、`--anyornone`、`--inner`、target/query `[multiple]` 等）的输出与原版二进制逐字节一致。
- **互斥报错**：以下组合都会报错退出（exit 1）：`--strand=plus`、`--strand=minus`、`--chain`、`--segments`、`--anyornone`、`--limitperquery`、`--chores=`、query 的 `[chores=]`、`--nogfextend`；不带 `--nlchain` 单独给 `--nlchainxxx` 也会报错。
- **人造序列**（target 两条染色体用 `[multiple]`；query = A[0:8k] + rc(A[8k:16k]) + B[5k:13k] + A[16k:24k]，8% 替换）：段表依次为 `S + chrA 0–8000`、`I − chrA 8000–16000`、`T + chrB 5000–13000`、`T + chrA 16000–24000`，query 和 target 坐标都与构造一致。由此确认反链 pos2 是反向互补坐标，换算 `qb = Q − pos2 − len` 正确。
- **同对角线正反链**：正链 HSP 和反链 HSP 都在对角线 0 上，pos2 区间 [1200,1500) 重叠，两者各自保留，链方向正确，未被合并。
- **query `[multiple]`**：按 query 分区运行的段表与逐条 query 运行完全一致。
- **pseudocat vs pseudopig**（`-Dcollect_stats` 构建）：22 个 HSP，主路径 17 个，6 段；jumpIn 为 S3/I0/T0/D2/G1，加起来正好 6 段；17 等于各段 HSP 数之和，也等于 gapped 阶段的 anchors 数。`lastz_D`、`--hspthresh=top10%`、`--queryhspbest`、`--nogapped`、全参数自定义、调低 T_sec 这几组都能正常运行。在 aglobin human vs cow 上，次级链回收到了珠蛋白簇的旁系拷贝。
- **断点处 query 重叠（§7 待定问题）**：只看主路径时，pseudopig 相邻 gapped 结果在 query 上的重叠为 0；aglobin 只有 1 处、52 bp。但 y-drop 会让 gapped 结果超出所在段的 query 区间，两组数据分别超出 6348 bp 和 2779 bp。数据量太小，**是否处理仍待第 4 步大数据再定**。

## 3. 与设计文档的差异（均已确认或有理由）

1. **剪枝用严格 `<`**（§6.1、§6.2 原文为 `≤`）：按原文等号剪枝，上界恰好等于当前最优时会剪掉可能按平局规则胜出的候选，做不到 §9 第 1 步"与穷举完全一致"。已经用户确认。
2. **L 也用 `nlpath`**（原文为 `double L; u32 predL`）：增加了 L 链的 HSP 数，使 L 与 U/R 共用同一个 `better()`。已经用户确认。
3. **段表 TSV**：首列为 query 名（已经用户确认）；target 列给序列名，不给 p 下标；坐标为各序列内从 0 起的半开区间，query 用正链坐标。为此 `nlsegment` 增加了 `part2`（query 分区下标）。
4. **压缩方式**：没有调用 `filter_marked_segments`，而是直接按输出段的次序重排保留下来的 HSP。原因是 §8.1 要求 `first/count` 是表内连续区间，而 `filter_marked_segments` 只保持原顺序。之后 `finish_one_strand` 会重新排序，对后续流程没有影响。
5. **统计**：沿用 lastz 的约定，只在 `-Dcollect_stats` 构建下收集。原版 `--stats` 在默认构建下本来就直接报"not implemented"，所以这样做不影响可见性。注意：现有的 `quantum.c` 在 `collect_stats` 下有一个 `-Werror` 报错（`lastIx` 未使用），这是原有问题，与本次改动无关；构建统计版时需要去掉 `-Werror`。
6. **"target 非 `[multiple]` 且有多条序列"的警告没有加**：这种情况下 lastz 本身就会报错退出（`lastz.c` "contains more than one sequence"），走不到 nlchain。
7. **`--nlchainxxx` 参数不会隐式打开 `--nlchain`**，单独使用时报错。
8. `g(d)` 用 double 做线性插值，不像 axtChain 那样截断为 int；这一点与 §5.1 的定义一致。
9. 有多个 query 分区时，返回值取各分区主路径分的最大值（与 `try_reduce_to_chain` 的做法一致）。

## 4. 未解决 / 风险

- **次级链的区域限额检查是 O(k²)**（k 为一个 query 分区内已保留的次级链数）。人–鼠全染色体上如果次级链达到 10^5 条，这一步可能耗时明显。上远端之前应先看 `secondary chains` 的统计，必要时改成区间索引。
- 默认构建的 `unspos` 只支持 2^31 的序列长度。全基因组 `[multiple]` target（人约 3.1 Gb）要用 `lastz_32` 构建。
- §7 的断点重叠是否处理仍待定（见上）。
- `README.lastz.html` 和 `lastz --help` 还没有写入新选项（等参数在 test set 上定下来再补）。

## 5. 下一步

1. 第 0 步（远端，不写代码）：在 Human–Mouse baseline MAF 上估算 query 端 1:1 和两端 1:1 的 recall 上限，并把 `--chain` 和 baseline 放到无旁系同源版真值上重评。上限不够就停下讨论。
2. 第 4 步（远端）：用 `lastz_32` 构建，target 用全基因组 `[multiple]`，每条 query 染色体一个作业。与 baseline_plain、`--chain`、LASTZ+axtChain/chainNet 比较 recall/precision/F/CPU/内存；同时看 `--stats` 里的 H_max 截断次数、次级链数和断点重叠量。
3. 第 5 步：在 test set 上调 §8.3 的参数并冻结，再跑 16 个物种对。
