# nlchain 实现进展

日期：2026-11-01。对应设计：`docs/nlchain_design.md`（下文 §x.y 均指该文档）。已完成 §9 的第 1–3 步（本地）、第 0 步（远端，结果见 §5）、第 4 步（远端，结果见 §6）和第 5 步的调参部分（远端，结果见 §8，默认值已冻结）；16 个物种对的最终运行尚未做。nlchain 代码已在 `nlchain` 分支提交（fa30eda，gcc 构建修复 7939842，第 4 步结果 5a1684d）；第 5 步的默认值改动和文档更新尚未提交。

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
--nlchainjump=<J>[,<λ>,<D0>,<L>]          默认 J=hspthresh, 1000, 3000, 8
--nlchaingap=<maxGap>,<maxDiag>[,<R>]     默认 50000, 10000, 0
--nlchaingapcost=medium|loose             默认 medium
--nlchainsolo=<T_solo>                    默认 1.25×hspthresh
--nlchainsource=<T_src>                   默认 0（关闭）
--nlchainoverlap=<f_ov>                   默认 0.5
--nlchainscan=<H_max>                     默认 1000
--nlchainsecondary=<T_sec>[,<N_sec>,<r_sec>]  默认 0, 1, 0.9；T_sec=0 关闭次级链（默认关闭；需要旁系拷贝时推荐 T_sec=4×hspthresh）
--nlchainreport=<file>                    段表 TSV
```

自适应 hspthresh（`--hspthresh=top...`）下，J、T_solo 中的 hspthresh 取本次两条链收集到的最低 HSP 分（在 merge 之前统计）。

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
- **断点处 query 重叠（§7 待定问题）**：只看主路径时，pseudopig 相邻 gapped 结果在 query 上的重叠为 0；aglobin 只有 1 处、52 bp。但 y-drop 会让 gapped 结果超出所在段的 query 区间，两组数据分别超出 6348 bp 和 2779 bp。数据量太小，当时留待第 4 步大数据再定；第 4 步的结论见 §6。

## 3. 与设计文档的差异（均已确认或有理由）

1. **剪枝用严格 `<`**（§6.1、§6.2 原文为 `≤`）：按原文等号剪枝，上界恰好等于当前最优时会剪掉可能按平局规则胜出的候选，做不到 §9 第 1 步"与穷举完全一致"。已经用户确认。
2. **L 也用 `nlpath`**（原文为 `double L; u32 predL`）：增加了 L 链的 HSP 数，使 L 与 U/R 共用同一个 `better()`。已经用户确认。
3. **段表 TSV**：首列为 query 名（已经用户确认）；target 列给序列名，不给 p 下标；坐标为各序列内从 0 起的半开区间，query 用正链坐标。为此 `nlsegment` 增加了 `part2`（query 分区下标）。
4. **压缩方式**：没有调用 `filter_marked_segments`，而是直接按输出段的次序重排保留下来的 HSP。原因是 §8.1 要求 `first/count` 是表内连续区间，而 `filter_marked_segments` 只保持原顺序。之后 `finish_one_strand` 会重新排序，对后续流程没有影响。
5. **统计**：沿用 lastz 的约定，只在 `-Dcollect_stats` 构建下收集。原版 `--stats` 在默认构建下本来就直接报"not implemented"，所以这样做不影响可见性。注意：现有的 `quantum.c` 在 `collect_stats` 下有一个 `-Werror` 报错（`lastIx` 未使用），这是原有问题，与本次改动无关；构建统计版时需要去掉 `-Werror`（`definedForAll` 里不带 `-Werror`）。
   - fa30eda 在服务器的 gcc 13 上默认构建失败：`if (currParams->nlchain) dbg_timing_report(...)` 在未定义 `dbgTiming` 时展开为空语句，触发 `-Werror=empty-body`（本地 clang 不报）。7939842 给这条 `if` 加了花括号，本地 `base_tests`、远端 `nlchain_test` 和 `base_tests` 都已通过。
6. **"target 非 `[multiple]` 且有多条序列"的警告没有加**：这种情况下 lastz 本身就会报错退出（`lastz.c` "contains more than one sequence"），走不到 nlchain。
7. **`--nlchainxxx` 参数不会隐式打开 `--nlchain`**，单独使用时报错。
8. `g(d)` 用 double 做线性插值，不像 axtChain 那样截断为 int；这一点与 §5.1 的定义一致。
9. 有多个 query 分区时，返回值取各分区主路径分的最大值（与 `try_reduce_to_chain` 的做法一致）。

## 4. 未解决 / 风险

- **次级链的区域限额检查是 O(k²)**（k 为一个 query 分区内已保留的次级链数）。在 simHuman–simMouse 上，7 条染色体合计 39,230 条次级链，整个 nlchaining 阶段只用 8.7 s，不是问题；真实人类基因组上次级链会更多，届时再看。
- 默认构建的 `unspos` 只支持 2^31 的序列长度。全基因组 `[multiple]` target（人约 3.1 Gb）要用 `lastz_32` 构建（Alignathon 的 simHuman 只有 191 Mb，不受影响）。
- ~~次级链按当前默认参数是净负收益~~：第 5 步已改为默认关闭（§8）。
- **与 chainNet 的 F 差距**（simTest 上 0.012，全量 simHuman–simMouse 上 0.004）调参补不上，主要来自 gapped 延伸越出段边界（§8）。
- §7 的断点重叠：第 4 步量级很小，建议暂不处理（§6）。
- `README.lastz.html` 和 `lastz --help` 还没有写入新选项（参数已冻结，可以补了）。

## 5. 第 0 步结果：1:1 选择的 recall 上限（远端，2026-10-31）

服务器 `MSA_Dedicated_Server_Remote:/mnt/sda/liuhaodong/workspace/experiment/lastz_profiling`。一条命令可重跑：`scripts/step0_recall_ub.sh`（已完成的阶段会跳过）。贪心选块和碱基对统计用 `scripts/maf_1to1.py`，结果在 `results/step0_recall_ub/`（`noparalogy/accuracy.tsv`、`ancestor/accuracy.tsv`、`stats/`）。输入是 10-26 那次全量运行的 MAF（target simHuman 按染色体，query simMouse 全基因组，UCSC 远缘参数）。

**无旁系同源版真值**
- Alignathon 页面上的链接 `simMammals.noparalogyMafs.maf.gz` 是 404。实际文件以官方 `downloadMammals.sh` 为准：`simMammals.noparalogyMafs.tar.gz`（1,458,185,918 B，md5 `bffc18321f937a1eac183c946049d190`，与官方脚本一致），里面有 `ancestor.noparalogies.maf` 和 `burnin.noparalogies.maf` 两个文件。这里只用前者（对齐到 MRCA），与之前评测用的 ancestor 真值口径一致。
- 用 `project_truth.py` 投影到 simHuman–simMouse，得到 99,958,410 个碱基对，与 `--check` 的独立计数一致（ancestor 真值是 157,142,590）。
- **侧别**：按碱基对统计，两个物种的每个碱基都只出现一次（重复使用 0），所以是**两端 1:1**。作为对照，ancestor 真值里 simHuman 有 1810 万、simMouse 有 903 万个碱基对重复使用了碱基。

**方法**：输入 `baseline_plain` 的 MAF（421,047 块），按块分数降序贪心选块。候选块与已选块并集的重叠不到自身长度的 50% 就接受，接受的块整块保留、不裁剪。q1 只看 simMouse（query）一侧，both 要求两侧都满足。mafComparator 抽 10M 样本，参数与 `eval.sh` 相同。

| 组 | 块数 | 预测碱基对 | 无旁系真值 recall | precision | F | ancestor 真值 recall | precision | F |
|---|---|---|---|---|---|---|---|---|
| baseline_plain | 421,047 | 2.72 亿 | 0.714 | 0.262 | 0.383 | 0.647 | 0.373 | 0.473 |
| `--chain` | 40,275 | 0.53 亿 | **0.293** | 0.550 | 0.382 | 0.236 | 0.698 | 0.353 |
| 贪心 q1（query 端 1:1） | 113,818 | 1.43 亿 | **0.697** | 0.487 | 0.574 | 0.626 | 0.688 | 0.656 |
| 贪心 both（两端 1:1） | 102,727 | 1.37 亿 | **0.691** | 0.505 | 0.584 | 0.611 | 0.703 | 0.654 |

- baseline 换种子重评：recall 0.71425 / 0.71418，precision 0.26173 / 0.26190，抽样误差可以忽略。
- **残余重叠**（按碱基对统计的重复使用，占该组碱基对的比例）：both 组 simHuman 210 万（1.5%）、simMouse 173 万（1.3%）；q1 组 simHuman 877 万（6.1%，即 target 端重复使用）、simMouse 235 万（1.6%）。`--chain` 自身也有 230 万 / 203 万。"不裁剪"让上限偏乐观的幅度在 1–2% 量级。

**结论**
1. **上限远高于 `--chain`，继续做**：同一份无旁系真值上，1:1 选择的 recall 上限约 0.69–0.70，`--chain` 只有 0.293，约为 2.4 倍。1:1 选择本身只损失 baseline recall 的 2–3%（0.714 → 0.697）。`--chain` 丢掉的 recall 绝大部分来自"只留一条共线链"，不是来自 1:1 约束，这正是 nlchain 主路径要救回的部分。
2. **侧别**：q1 比 both 多出约 627 万个碱基对，其中只有约 60 万是真同源（recall 只高 0.006），precision 低 1.8 个百分点。也就是说，target 端重复使用（对应 nlchain 主路径经 D 回跳）在两端 1:1 的真值上大多算成假阳性。
3. 注意这个上限是"从 baseline 全部 gapped 结果里挑"。nlchain 在 gapped 之前筛 HSP，实际 recall 只会更低，不会更高。

## 6. 第 4 步结果：simHuman–simMouse 全量比较（远端，2026-11-01）

**运行方式**
- RUN_ID：`nl_simHuman-simMouse_baseline_plain+chain+nlchain+nlchain_primary+nlchain_prof+nlchain_primary_prof_all_20261031T104720Z`。
- 输入：target 为 `simHuman.fa[multiple]`，query 为单条 simMouse 染色体，7 条染色体 × 6 组，共 42 个作业，每个作业独占一个核。
- 参数：UCSC 远缘参数（`--hspthresh=2200 --gappedthresh=6000 --ydrop=3400 --masking=50 --scores=HoxD55.q --ambiguous=iupac`）。nlchain 其余参数全部取 §1 的默认值；`nlchain_primary` 加 `--nlchainsecondary=0`。
- 构建：二进制由 7939842 构建。`bin/lastz_nl` 是正式版（`-g`）；`bin/lastz_nl_prof` 额外加了 `-Dcollect_stats -DdbgTiming -DdbgTimingGappedExtend`，只用于计数和分阶段计时。两版的 MAF 已逐组核对一致。
- LASTZ+chainNet：在本次 `baseline_plain` 的 MAF 上，按 UCSC doBlastzChainNet 的标准流程，用 kent 工具跑了两组参数：
  - 主：`-linearGap=medium -minScore=3000`；
  - 对照：`-linearGap=loose -minScore=5000`。

  流程是 `axtChain -scoreScheme=HoxD55 | chainAntiRepeat → chainMergeSort → chainPreNet → chainNet -minSpace=1 | netSyntenic → netToAxt | axtSort → axtToMaf`。kent 工具来自 hgdownload 的 linux.x86_64 预编译版，sha256 记在 `bin/kent/MANIFEST.tsv`。MAF→axt 由 `scripts/maf_axt.py` 转换，7 条染色体都做了往返校验，0 差异。
- 评测：mafComparator 抽 10M 样本，种子为 1；baseline 换种子 2 重评，recall 0.70914 / 0.70906，抽样误差可以忽略。
- 脚本（服务器 `scripts/`）：
  - `bootstrap_nl.sh`：构建、测试、冒烟；
  - `run_nl.sh` / `run_one_nl.sh`：运行全量；
  - `chainnet.sh`：LASTZ+chainNet；
  - `eval_nl.sh`：评测；
  - `nl_breakpoint.py`：断点统计；
  - `nl_report.py`：汇总。

  汇总表在 `results/runs/<RUN_ID>/report_nl.md`。两组 `chainnet.sh` 并行跑时，会同时生成同一个 `.2bit` 而产生竞争；应先串行生成 2bit，或者把两组串行跑。

**准确度与资源**（CPU = user+sys，7 个作业求和；chainNet 的 CPU = baseline lastz CPU + 链化/建 net 各步 CPU）

| 组 | CPU s | 相对 baseline 提速 | 最慢作业墙钟 s | 峰值内存 GB | MAF 块 | 无旁系 R | P | F | ancestor R | P | F |
|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline_plain | 55,151 | 1.00× | 22,412 | 1.65 | 528,729 | 0.709 | 0.228 | 0.345 | 0.642 | 0.324 | 0.431 |
| `--chain` | 32,445 | 1.70× | 11,473 | 1.38 | 39,425 | 0.290 | 0.551 | 0.380 | 0.234 | 0.700 | 0.351 |
| **nlchain_primary** | 32,762 | 1.68× | 11,580 | 1.73 | 102,627 | 0.700 | 0.511 | 0.591 | 0.619 | 0.712 | 0.662 |
| nlchain（含次级链） | 32,484 | 1.70× | 11,681 | 1.73 | 139,868 | 0.701 | 0.453 | 0.550 | 0.628 | 0.637 | 0.632 |
| **chainNet medium/3000** | 55,226 | 1.00× | — | 1.65 | 131,235 | 0.702 | 0.516 | 0.595 | 0.620 | 0.718 | 0.666 |
| chainNet loose/5000 | 55,228 | 1.00× | — | 1.65 | 137,546 | 0.700 | 0.516 | 0.594 | 0.619 | 0.717 | 0.665 |

**碱基对重复使用**（被多于一个碱基对用到的碱基，占该组碱基对的比例）：

| 组 | simHuman 端 | simMouse 端 |
|---|---|---|
| baseline | 54.9% | 54.2% |
| `--chain` | 4.2% | 3.7% |
| nlchain_primary | 4.9% | 3.1% |
| nlchain | 13.9% | 11.9% |
| chainNet medium | 0（target 端的 net 保证了这一点） | 3.6% |

**分阶段计时与计数**（剖析版，7 条染色体求和）

| | seed hit search | nlchaining | gapped extension | 总计 |
|---|---|---|---|---|
| nlchain_primary | 31,937 s | 7.5 s | 528 s | 32,651 s |
| nlchain | 31,089 s | 8.7 s | 711 s | 31,974 s |

- **主路径**：HSP 5,424,209 个，其中 399,337 个（7.4%）进入主路径，分成 25,698 段。jumpIn 为 S 7 / I 12,862 / T 3,617 / D 4,108 / G 5,104。S 正好对应 7 条 query，说明每条 query 只有一条主路径；I（倒位）占一半。
- **次级链**：39,230 条，含 87,660 个 HSP。H_max 截断 0 次。
- **gapped 阶段**：DP 单元 nlchain_primary 为 5.66e10，nlchain 为 7.17e10。
- **chainNet 的链化和建 net**：合计只用 75 s CPU，其中 axtChain 占 40 s。
- **冒烟阶段（chrO）的对比**：baseline 的 gapped extension 为 254.5 s，nlchain_primary 只有 9.5 s。全量没有跑 baseline 的剖析版，但两者 CPU 差 22,389 s，与"省掉的就是 gapped 时间"相符。§9 的三个停止阈值都没有触发：
  - 次级链的追溯和检查合计不到总时间的 0.03%（阈值 2%）；
  - nlchain 的 gapped 时间约为 baseline 的 5%（阈值 50%）；
  - 截断次数为 0（阈值 1%）。

**断点处 query 重叠**（nlchain_primary，`nl_breakpoint.py`）
- 102,627 个块全部能归到主路径段（按 query、target、链方向匹配，取 query 重叠最大的段）。
- **段间**：25,691 个相邻段对中，有 5,940 个（23%）两侧的块在 query 上重叠。重叠合计 0.88 Mbp，占比对上的 query 碱基（152.8 Mbp）的 0.57%；中位数 64 bp，最大 4,891 bp。
- **全体**：被 ≥2 个块覆盖的 query 碱基共 5.05 Mbp（3.3%），大部分发生在段内，不在断点处。
- **越界**：35% 的块超出了所属段的 query 区间，超出部分合计 12.0 Mbp。原因是 y-drop 延伸越过了段的边界。

**结论**
1. **主要目标达成**：nlchain_primary 的准确度与 LASTZ+chainNet（medium/3000）基本持平。
   - 无旁系真值：F 0.591 对 0.595；ancestor 真值：F 0.662 对 0.666。precision 低约 0.005，这超出了抽样误差，但幅度很小。
   - CPU 是后者的 1/1.68，墙钟约为一半。
   - 与 `--chain` 相比，在 CPU 几乎相同的情况下，无旁系 recall 从 0.290 提到 0.700（2.4 倍）。
   - chainNet 的两组参数几乎没有差别（F 相差 0.001），所以下面只用 medium/3000 比较。
2. **加速的上限现在在 seed 阶段**：nlchain 的总时间里，seed hit search 占 97.8%，gapped 只占 1.6%，chaining 可以忽略。要继续提速只能动 seed 阶段（例如 seed 模式、步长），这超出了 nlchain 的范围。
3. **次级链按当前默认参数（T_sec = 2×hspthresh，N_sec = 5，r_sec = 0.5）是净负收益**。
   - 代价很小：gapped 时间只增加 183 s，占总时间的 0.57%。
   - 但收益为负：无旁系 recall +0.001，precision −0.058；即使在含旁系同源的 ancestor 真值上，也是 recall +0.009，precision −0.075，F 下降 0.030。
   - 重复使用：simHuman 端从 4.9% 升到 13.9%。也就是说，回收到的大多不是真值认可的同源对。
   - 第 5 步要么把默认改为关闭，要么在 test set 上收紧 T_sec / r_sec。
4. **断点重叠暂不处理**：断点处的重叠只占比对碱基的 0.57%。即使全部裁掉，precision 的上限收益也不到 1 个百分点。nlchain_primary 的 simMouse 端重复使用（3.1%）已经低于 chainNet（3.6%）。
5. **评测口径的副作用：运行方式会改变 baseline 和 `--chain` 的结果**。
   - 第 0 步用的是旧运行方式：每条人类染色体对整个小鼠基因组跑一次。在无旁系真值上，baseline 为 0.714 / 0.262，`--chain` 为 0.293 / 0.550。
   - 本次每个作业只有一条 query，baseline 变成 0.709 / 0.228，预测碱基对从 2.72 亿增到 3.11 亿；`--chain` 变成 0.290 / 0.551。
   - 原因是 `--masking` 的动态掩蔽只影响同一次运行中后续的 query（见 `lastz.c` 开头的注释），作业拆分方式不同，掩蔽的效果就不同。
   - 因此设计 §9 的"`--chain` 与逐染色体运行结果相同"忽略了掩蔽。各方法之间的比较必须在同一种运行方式下进行（本节已经这样做）。
   - nlchain_primary 的 recall 0.700 也略高于第 0 步的贪心上限 0.691，主要也是因为两者的输入不同。另外 nlchain 的 gapped 结果并不是 baseline 结果的子集：anchor 少了以后，延伸不再被已有比对挡住。所以第 0 步的上限只能作为参考，不是严格的上界。

## 7. 下一步

1. **评测口径（已定，2026-10-31）**：第 4 步三方（nlchain、`--chain`、LASTZ+axtChain/chainNet）都用**原始输出**，直接用 mafComparator 在对应真值上评，**不做任何 1:1 后过滤，也不做 netSwap**。理由：无旁系真值是两端 1:1，任何一侧的重复使用都同样算偏差，三方按各自实际输出被扣分，比较是公平的。`maf_1to1.py` 只用于第 0 步估上限，不进入评测流程。第 0 步的数据可以用来解读结果：target 端重复使用大约让 precision 低 2 个百分点。
2. ~~第 4 步~~：已完成，见 §6。
3. ~~第 5 步调参~~：已完成，默认值已冻结并写入 `lastz.c`，见 §8。
4. **提交**：`src/lastz.c`、`src/lastz.h` 的默认值改动，以及 `docs/` 的更新（尚未提交）。
5. **第 5 步后半：跑 16 个物种对**。参数按 §8 冻结，不再改动。
   - 运行方式与第 4 步相同（target `[multiple]`，每条 query 染色体一个作业），每个物种对同时跑 LASTZ+chainNet medium/3000 作为对照；
   - 服务器的 `run_nl.sh`、`chainnet.sh`、`eval_nl.sh` 写死了 simMammals 和 simHuman–simMouse，需要按 `run_tune.sh` 的方式参数化，还要接入灵长类集（数据尚未下载）；
   - 近缘物种对没有调过参，结果要单独看。
6. **候选算法改动**：按段的 query/target 区间裁剪 gapped 输出，用来补与 chainNet 的 precision 差距（§8 结论 3）。这属于算法改动，需要先评估、再决定做不做。
7. 补 `README.lastz.html` 和 `lastz --help`。

## 8. 第 5 步结果：simTest 调参（远端，2026-11-01）

详见 `docs/nlchain_tune_testset.md`。

**数据与方法**
- Alignathon simTest：5 个物种，每个约 11 Mb，树与哺乳类集相同，10 个物种对的距离为 0.176–0.606，没有近缘物种对。
- 两份真值投影都与独立计数一致；noparalogy 真值是两端 1:1。
- 运行方式与第 4 步相同，用 UCSC 远缘参数（h=2200）。
- 共跑了 79 组参数，分 4 个阶段：主路径单参数扫描、λ×T_solo×gap 表×D0 组合、次级链网格、候选组合复核。每组都显式写全全部 nlchain 参数。
- 评测：mafComparator 10M 样本，两份真值各评一次，取 10 个物种对的平均。
- 服务器脚本：`prepare_test.sh`、`run_tune.sh`（`TUNE_BIN` 可换二进制）、`chainnet_pair.sh`、`eval_tune.sh`、`tune_report.py`、`gen_tune_groups.py`、`diag_1to1.sh`；结果在 `results/tune/`。

**结论**
1. **冻结值**（用户确认方案 A）：λ 从 500 改为 1000，D0 从 1000 改为 3000，T_solo 从 2h 改为 1.25h；次级链默认关闭（T_sec=0，N_sec=1，r_sec=0.9）；其余不变。
   - noparalogy 真值上的平均 F 从 0.7060 升到 0.7085，10/10 个物种对都提高；ancestor 真值上的平均 F 从 0.7069 降到 0.7054。
   - 默认值写入 `lastz.c` 后，服务器构建的二进制只给 `--nlchain` 时，20 个作业的输出与对应调参组逐字节一致。
2. **主路径参数整体不敏感**：单参数扫描的 F 变化都不超过 ±0.005。maxGap、maxDiag 完全不影响结果。主要取舍在 λ 和 J：跳转越贵，precision 越高、recall 越低。
3. **与 chainNet 的差距调参补不上**：F_np 0.7085 对 0.7205，recall 相同，差距全在 precision，10 个物种对都是这样。
   - 诊断：两边都做两端 1:1 过滤后，precision 差距只从 0.026 缩到 0.021，碱基重复使用只解释了约 1/5。
   - 其余部分最可能来自 gapped 延伸越出段边界（第 4 步测得 35% 的块越界）。
4. **次级链与原判据的冲突**：原判据是"ancestor 真值上的 F 不降就默认开启"。按这条，T_sec ≥ 3h、N_sec=1、r_sec=0.9 的几组都达标（F_anc +0.001 ~ +0.003，10/10 个物种对）。但这些组在每个物种对上都降低 noparalogy F（−0.002 ~ −0.006），所以用户决定默认关闭。需要旁系拷贝时推荐 `--nlchainsecondary=<4h>,1,0.9`。
5. **口径提醒**：§6 结论 3 的"次级链净负收益"来自 simHuman–simMouse，它是 16 个最终测试物种对之一，不能用来定默认值。本节的决定只依据 simTest。
