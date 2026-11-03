# nlchain 实现进展

日期：2026-11-01，2026-11-02 更新。对应设计：`docs/nlchain_design.md`（下文 §x.y 均指该文档）。已完成 §9 的第 1–3 步（本地）、第 0 步（远端，结果见 §5）、第 4 步（远端，结果见 §6）和第 5 步的调参部分（远端，结果见 §8，默认值已冻结并在 06533fd 提交）；16 个物种对的最终运行尚未做。nlchain 代码已在 `nlchain` 分支提交（fa30eda，gcc 构建修复 7939842，第 4 步结果 5a1684d，默认值 06533fd）。

**2026-11-02**：方案 A（gapped 之后重新串链 + 链分下限，`docs/nlchain_rechain_plan.md`）的本地实现和验证已完成（实施方案 §4 第 1–6 步），见 §9；尚未提交，远端对拍和评测（第 7–8 步）尚未做。

**2026-11-03**：nlchain 模式下 gapped 之前做正反两链插值（`docs/nlchain_prechain_inner_plan.md`）的本地实现和验证已完成（该方案 §7 第 1–4 步），见 §11，已在 19fa9e5 提交；远端实验（第 5 步，该方案 §8）已完成，见 §12：`--inner=2000` 使 F_np/F_anc 各升约 0.01，再把 gapped 门槛降到 2000 再升约 0.009；补跑的 G3a（D6 取 a 并降门槛）F_np 0.7422、F_anc 0.7355，且没有可见开销。据此 D6 默认改为 a（第二次串链沿用外层 J、T_solo），gapped 门槛不改代码、在文档中推荐与 `--inner` 同值（§12.6）。

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
--nlchainclean=<T_chain>                  默认 5000；gapped 后重新串链，删掉链分 < T_chain 的链；0 关闭（§9）
--nlchaincleanreport=<file>               剥链结果，kent chain 格式（删除前的全部链）
```

`--nlchaingapcost` 同时决定重新串链用的 axtChain 代价表（medium / loose 的三行表）。

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
- **与 chainNet 的 F 差距**（simTest 上 0.012，全量 simHuman–simMouse 上 0.004）调参补不上。根因已定位：不是越出段边界，而是 nlchain 缺少 chainNet 流程中"gapped 后拆块重新串链、丢掉低分碎片链"这一步，见 `docs/nlchain_precision_gap.md`。
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
4. ~~提交默认值改动~~：已在 06533fd 提交。
5. **第 5 步后半：跑 16 个物种对**。参数按 §8 冻结，不再改动。
   - 运行方式与第 4 步相同（target `[multiple]`，每条 query 染色体一个作业），每个物种对同时跑 LASTZ+chainNet medium/3000 作为对照；
   - 服务器的 `run_nl.sh`、`chainnet.sh`、`eval_nl.sh` 写死了 simMammals 和 simHuman–simMouse，需要按 `run_tune.sh` 的方式参数化，还要接入灵长类集（数据尚未下载）；
   - 近缘物种对没有调过参，结果要单独看。
6. **补 precision 差距**：根因与方案见 `docs/nlchain_precision_gap.md`，实施方案见 `docs/nlchain_rechain_plan.md`。
   - 原候选"按段的 query/target 区间裁剪"已否定：模拟后 F 下降。
   - **已确认（2026-11-02）**：方案 A 在 `--nlchain` 下默认开启，即 gapped 之后按无 gap 块重新串链（照搬 axtChain 的算法），删掉链分 < 5000 的链。
   - 次级链按用途选：只要 1:1 直系同源用默认参数；需要旁系同源加 `--nlchainsecondary=<4h>,1,0.9`。
   - simTest 上的模拟（10 个物种对平均）：

     | 配置 | F_np | F_anc |
     |---|---|---|
     | 默认 + A | 0.7247 | 0.7167 |
     | 次级链 + A | 0.7208 | 0.7194 |
     | chainNet | 0.7205 | 0.7191 |

   - ~~编码、单测、本地回归~~：已完成（2026-11-02），见 §9。D1–D6 与补充细节 X1–X6 均按推荐确认。
   - ~~远端第 7 步（与 axtChain 对拍）、第 8 步（两套配置评测）~~：已完成（2026-11-02），全部达标，见 §10。
   - 步骤 5（16 个物种对）：方案 A 已通过 simTest 复核，何时跑由用户决定。
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
   - 其余部分当时推测来自 gapped 延伸越出段边界（第 4 步测得 35% 的块越界）。这一推测后来被否定，见 `docs/nlchain_precision_gap.md`：越界碱基对有 58% 是真阳性；差距来自缺少"重新串链 + 链分下限"这一步。
4. **次级链与原判据的冲突**：原判据是"ancestor 真值上的 F 不降就默认开启"。按这条，T_sec ≥ 3h、N_sec=1、r_sec=0.9 的几组都达标（F_anc +0.001 ~ +0.003，10/10 个物种对）。但这些组在每个物种对上都降低 noparalogy F（−0.002 ~ −0.006），所以用户决定默认关闭。需要旁系拷贝时推荐 `--nlchainsecondary=<4h>,1,0.9`。
5. **口径提醒**：§6 结论 3 的"次级链净负收益"来自 simHuman–simMouse，它是 16 个最终测试物种对之一，不能用来定默认值。本节的决定只依据 simTest。

## 9. 方案 A 实现：gapped 之后重新串链（本地，2026-11-02）

实施方案：`docs/nlchain_rechain_plan.md`（下文"方案 §x"指该文档）。完成了方案 §4 的第 1–6 步，代码尚未提交。

### 9.1 已完成

| 文件 | 内容 |
|---|---|
| `src/nlrechain.h` / `src/nlrechain.c`（新） | 方案 §2 第 1–8 步：拆块、合并完全重复块、块分、DP（2-d 树）、剥链、链内重叠裁剪、链分、下限、按 X2 分段重建 `alignel`；kent chain 格式报告；`-Dcollect_stats` 计数 |
| `src/nlrechain_test.c`（新） | 单元测试，直接 `#include "nlrechain.c"`；`gapped_extend.o` 引用的过滤/输出函数用桩函数代替 |
| `src/lastz.h` / `src/lastz.c` | `nlrechainparams nlr`、默认值、`--nlchainclean`、`--nlchaincleanreport`、`--nlchaingapcost` 同时选三行表、互斥检查、`debugClockNlRechain`、统计；调用点在 `finish_one_strand` 的过滤之后、`tweener_interpolate` 之前 |
| `src/Makefile` | `nlrechain` 加入构建；新增 `make nlrechain_test` |

用户确认的设计选择（2026-11-02）：D1 照搬 axtChain；D2 在被删块处切开原比对分段输出；D3 `--nlchainclean=<T_chain>` 默认 5000、0 关闭；D4 绝对值；D5 不做 chainAntiRepeat 的重复/低复杂度调整；D6 动态掩蔽用删碎片之后的比对。补充细节：
- X1：完全重复块合并后归较长块的来源比对（等长取输入次序靠前的）。
- X2：链内相邻、来源比对中次序相邻、接缝未被裁，三者都满足才合成一条比对，块之间沿用原比对的 indel 操作。
- X3：输出按来源比对的输入次序排列，同一来源内按块次序。
- X4：平局规则，见 9.3。
- X5：`--nogapped` 下跳过；显式给 `--nlchainclean`/`--nlchaincleanreport` 时报错。
- X6：做 `--nlchaincleanreport`。
- 对拍门槛：Jaccard ≥ 0.995。

### 9.2 验证结果

- **单元测试** `make nlrechain_test`：默认 3000 轮，约 292 万项检查，0 失败（macOS clang 和 ubuntu gcc 11 都通过）。
  - gap 表：medium/loose 三行的 11 个节点；插值截断、外推、bothGap 按 dq+dt、负值归零的手算值；与按 kent 描述（1..110 预计算、浮点插值后截断）独立写的参考实现在 dq、dt ∈ [0,400]² 及 4 万个随机大值上逐点一致；gapCalcCost 对两个参数单调不减（搜索层下界依赖这一点）。
  - 交叉点与 connectCost：手算样例，包括调整量为负（−228）、三个位置平局取 pos=0、包含时加 1e8。
  - 拆块：随机 edit script（含 ins 紧接 del、相邻的两个 sub 操作），块的碱基对集合等于比对的碱基对集合；块分之和减去 gap 罚分等于 `score_alignment`。
  - DP：2-d 树与暴力枚举（每个 i 枚举全部 j）逐块比较 total 和 pred，四档规模（n ≤ 12 到 n ≈ 2500）、两张表、1–4 个字母的字母表（制造平局），全部一致；共 43.7 万条前驱边，其中 7.9 万条的交叉点调整为负。
  - offer/settle 的规则单独手算（暴力解与 2-d 树共用这两个函数，对拍覆盖不到）：候选必须严格大于 0 且严格大于块分本身；同分取下标小的前驱。
  - 剥链与链内整理：分别与按 kent 写法独立实现的参考解比较，全部一致；覆盖了交叉点裁剪（1.2 万块）和包含时删块（2400 条链）。
  - 输出重建：输出比对的碱基对（按出现次数）等于保留链的块的碱基对；每条输出比对的分数等于 `score_alignment`；单条强比对原样输出（beg/end/script/分数/hspId 不变）；X1、X3 的具体例子；下限高于链分时输出为空。
  - 变异检验：把剪枝条件改成 `≤`、去掉下界中的负调整项、去掉平局取小下标、交叉点平局改取靠后、bothGap 改查 qGap、去掉"高于块分本身"、去掉"高于 0"，测试分别报出 17、4308、442、1、33 万、1、18 项失败。唯一测不出的是 X2 的"接缝未被裁"条件：同一比对中相邻的块互不重叠，链内整理不会裁这个接缝，这两个条件实际总是成立（代码注释已说明）。
- **回归**（macOS clang 与 ubuntu gcc 11 各跑一遍，对照 06533fd 构建的二进制）：
  - `make base_tests` 19 项全过；`make nlchain_test` 0 失败。
  - 不加 `--nlchain` 的 22 组选项（lav/maf/axt/general/sam/cigar/segments、`--chain`、`--chain=20,1`、`--nogapped`、`--hspthresh=top10%`、`--queryhspbest`、`--queryhsplimit`、`--strand=plus/minus`、`--anyornone`、`--inner`、target/query `[multiple]`、`--masking`、`--identity/--coverage`、UCSC 远缘参数）：stdout、stderr、退出码逐字节一致。
  - `--nlchain --nlchainclean=0` 与 06533fd 的 `--nlchain` 在 8 组选项上逐字节一致（默认、maf、次级链、自适应阈值、UCSC 远缘参数、target `[multiple]`、loose 表、`--inner`）。输出头部回显的命令行多了 ` --nlchainclean=0`，比较前去掉。
  - gcc 11 加 `-Werror` 构建 `lastz`、`lastz_D`、两个单测均无警告。
- **互斥报错**：不带 `--nlchain` 单独给 `--nlchainclean`/`--nlchaincleanreport`、`--nlchain --nogapped --nlchainclean=…`、负的 T_chain，都报错退出；`--nlchain --nogapped`（不显式给）正常运行，跳过这一步。
- **冒烟**（`-Dcollect_stats -DdbgTiming` 构建）：pseudocat vs pseudopig、aglobin human vs cow，5 组参数下计数都自洽：输入碱基对 = 完全重复合并掉的 + 链内裁掉的 + 删除链上的 + 保留的。
- **小数据与 axtChain 对拍**（本地 OrbStack ubuntu VM 直接运行远端的 kent linux 二进制，sha256 与 `bin/kent/MANIFEST.tsv` 一致）：先用 `--nlchain --nlchainclean=0` 输出 axt，再用 `axtChain -minScore=0 -linearGap=medium -scoreScheme=HoxD55.axtChain.q` 串链。与 `--nlchain --nlchaincleanreport` 比较：
  - 5 组（catpig、aglobin，各含 UCSC 参数、次级链 8800,1,0.9、默认参数）的链分 ≥ 5000 的链逐条相同（块和链分都一致），输出 axt 的碱基对与 axtChain ≥5000 的链的 Jaccard 都是 1.000000。
  - 不设下限时，链分 ≥ 0 的全部链也逐条相同。报告里多出的只有负分链，axtChain 的 `-minScore=0` 本身就会丢掉它们。

### 9.3 与实施方案的差异（均有理由）

1. **数值类型用 double，不用 s64**（X4 原定 s64）。`lastz_D` 构建的 `score` 是 double；整数矩阵下所有分数都是 2^53 以内的整数和，double 运算精确，与 s64 等价。gap 代价仍按整数截断（`gap_row_cost` 用整数算术，与 kent 的浮点插值截断逐点一致）。
2. **平局与前驱规则（X4 的具体化，均按 kent 行为）**：DP 按 (tStart, qStart) 次序处理；前驱候选必须严格大于 0（kent 搜索初值为 0）且严格大于块分本身；同分取 DP 次序靠前的前驱；剥链按 total 降序，同分按 DP 次序。kent 的同分次序取决于树的遍历顺序，这一点不保证相同。
3. **2-d 树的剪枝下界比 kent 严格**。kent 用 gap 代价作为 connectCost 的下界，但交叉点调整可能为负（等于被舍弃各列的替换分之和），所以 kent 的剪枝在极少数情况下会漏掉最优前驱。这里的下界加上了"重叠列数上界 × 最低单列分"一项，保证与暴力枚举完全一致。因此在 kent 漏解的情况下，本实现可能选到更好的前驱。小数据对拍中没有出现这种差异，第 7 步的 Jaccard 会反映它的总体影响。
4. **报告写出全部链**，包括负分链；axtChain `-minScore=0` 会丢掉负分链，对拍时按链分 ≥ 0 比较。
5. **链内整理的"裁空"分支实际不可达**：包含（ov ≥ 任一块长）时已整块删除，否则交叉点不超过 ov，裁剪后两块都非空。代码仍按 kent 的流程保留这两个分支。
6. **首尾相接块的合并（方案 §2 第 6.2 步）只在报告里做**：gapCalcCost(0,0) = 0、块分可加，合并不改变链分；输出比对按 X2 分段，不需要合并。

### 9.4 未解决 / 风险

- 性能只在小数据上看过（`nlchain rechaining` 计时为 0.000 s），全量规模要在第 8 步确认 CPU ≤ lastz 总 CPU 的 1%。
- 与模拟的已知差别仍在：N 等非 ACGT 碱基的计分（lastz `--ambiguous=iupac` 对 axtChain 的 HoxD55 矩阵），以及 D6 的动态掩蔽效应。
- `README.lastz.html` 和 `lastz --help` 尚未写入 nlchain 的任何选项（原计划在方案 §4 第 10 步一并补）。
- 单测的随机流中有一处 `rnd() * rnd()` 的求值顺序由编译器决定，clang 和 gcc 下的随机用例不同（两边都 0 失败）。

### 9.5 下一步（新会话）

1. **提交**：本地实现与文档，先给用户看文件列表、diff 摘要和 commit message。
2. **方案 §4 第 7 步（远端对拍）**：
   - 在服务器上从提交后的分支构建新二进制（gcc，`-Werror`），先跑 `make nlrechain_test`、`make nlchain_test`、`make base_tests`。
   - 用 `TUNE_BIN` 换成新二进制，在 simTest 20 个作业上跑 `--nlchain --nlchainclean=0 --nlchaincleanreport=…`；其输出应与 `runs/tune/nl_step5default/` 逐字节一致（只差命令行回显）。
   - 再与 `runs/prec_sim/nl_step5default.medium0/<物种对>/chain_raw/` 逐链比较，以及与 `runs/prec_sim/xover/<物种对>/raw5000.chain` 比较碱基对，要求 Jaccard ≥ 0.995。
   - 注意：对拍要用 `--nlchainclean=0` 的输出作为 axtChain 的输入，因为 D6 会让默认运行的后续 query 掩蔽不同。本地对拍脚本的做法见 9.2。
3. **方案 §4 第 8 步（远端评测）**：默认参数和 `--nlchainsecondary=8800,1,0.9` 各跑 20 个作业（原始输出，mafComparator `--samples=10000000 --seed=1`，两份真值），并记录 CPU。验收线：F 与模拟值相差 ≤ 0.002；默认配置 F_np > 0.7205；新增阶段 CPU ≤ 1%。
4. 第 9 步（可选，在 A 之上重扫 T_sec）、第 10 步（README、`--help`、design 文档收尾）。16 个最终测试物种对等第 8 步通过后由用户决定。


> 9.5 的第 1–3 项已完成：实现在 8e9f350 提交；第 7、8 步见 §10。第 2 项写的"`--nlchainclean=0 --nlchaincleanreport=…`"行不通：T_chain=0 时整个重新串链阶段被跳过（`lastz.c` 中 `minChainScore > 0` 才调用），报告为空。实际改用 `--nlchainclean=1`，见 10.2。

## 10. 方案 A 远端复核：第 7–8 步（远端，2026-11-02）

服务器 `MSA_Dedicated_Server_Remote:/mnt/sda/liuhaodong/workspace/experiment/lastz_profiling`，设 `PROJECT_ROOT` 为该目录。结果汇总在 `results/nlr/`。

### 10.1 构建与测试

- 源码：`git archive 8e9f350` 解到 `external/lastz-8e9f350`。构建脚本 `scripts/bootstrap_nlr.sh`，仿照 `bootstrap_nl.sh`，不覆盖已有二进制，sha256 记在 `bin/MANIFEST_nlr.tsv`。
  - `bin/lastz_nlr`：正式版，默认 Makefile，含 `-Werror`，gcc 13.1，sha256 `13271116…89ddc6`。
  - `bin/lastz_nlr_prof`：`-Dcollect_stats -DdbgTiming`，sha256 `38e8a301…dbf86`。这一版不加 `-Werror`：`-Dcollect_stats` 下上游 `pos_table.c:2193/2195` 的 `sprintf` 会触发 gcc 13 的 `-Wformat-overflow`，与本次改动无关。脚本另外检查 `nl*.c/h` 没有警告。
- 测试（在源码副本 `build/nlr_tests/` 中运行）：`make nlrechain_test` 共 2,929,359 项检查，0 失败；`make nlchain_test` 共 518,512 项检查，0 失败；`make base_tests` 19 项全部 SUCCESS；`make lastz_D` 加 `-Werror` 通过。正式版、单测、`lastz_D` 的构建日志都没有 warning。
- 冒烟：两版在 pseudocat vs pseudopig 上的 `--nlchain` 输出和清洗报告逐字节一致。

### 10.2 第 7 步：与 axtChain 对拍（通过）

- 参数组加在 `config/tune_groups.tsv` 末尾，用 `TUNE_BIN=bin/lastz_nlr` 运行：`nlr_clean0`、`nlr_clean1`、`nlr_default`、`nlr_sec`、`nlr_default_prof`。注意：重新运行 `gen_tune_groups.py` 会覆盖这些行。
- `scripts/run_tune.sh` 只加了一行：组参数含 `--nlchainclean=` 时，自动加 `--nlchaincleanreport=<作业>.clean.chain`。原版备份在 `logs/run_tune.sh.bak_nlr`。
- **`nlr_clean0`（`--nlchain --nlchainclean=0`）**：20 个作业的 MAF（去掉 `#` 行）和 `segs.tsv` 都与 `nl_step5default` 逐字节一致，说明关闭开关就回到旧行为。这一组的清洗报告为空，原因见 9.5 后的说明。
- **对拍用 `nlr_clean1`（`--nlchainclean=1`）**：报告写出全部链，只删掉链分 < 1 的链；axtChain 的 `-minScore=0` 本来也会丢掉这些链。20 个作业的 `segs.tsv` 与 `nl_step5default` 逐字节一致。
- 对拍脚本 `scripts/nlr_vs_axtchain.py <组> <物种对>`：
  - 逐作业比较 lastz 报告中链分 ≥ 0 的链与 `chain_raw/<作业>.chain`。比较内容是头部（不含 id）、块列表和链分，按多重集比较。
  - 按物种对计算链分 ≥ 5000 的碱基对集合与 `xover/<物种对>/raw5000.chain` 的 Jaccard。
  - 结果在 `results/nlr/step7.tsv`，差异链示例在 `results/nlr/step7_diff.log`。

| 物种对 | 链 ≥0 一致 / kent | 链 ≥5000 一致 / kent | Jaccard（≥5000 碱基对） |
|---|---|---|---|
| Mouse-Rat | 10671 / 10702 | 889 / 900 | 0.999894 |
| Cow-Dog | 15738 / 15770 | 968 / 979 | 0.999870 |
| Dog-Human | 16216 / 16263 | 946 / 961 | 0.999875 |
| Cow-Human | 15729 / 15761 | 950 / 960 | 0.999872 |
| Human-Mouse | 11255 / 11298 | 762 / 778 | 0.999838 |
| Human-Rat | 11181 / 11218 | 671 / 684 | 0.999725 |
| Dog-Mouse | 9012 / 9048 | 522 / 533 | 0.999869 |
| Dog-Rat | 7892 / 7921 | 455 / 464 | 0.999914 |
| Cow-Mouse | 6786 / 6807 | 472 / 481 | 0.999921 |
| Cow-Rat | 6362 / 6385 | 433 / 442 | 0.999829 |
| 合计 | 110842 / 111173（99.70%） | 7068 / 7182（98.4%） | 0.99973–0.99992 |

- **Jaccard 全部 ≥ 0.9997，远高于门槛 0.995**。另有一项自检：kent 一侧按作业统计的 ≥5000 链数与 `raw5000.chain` 的链数逐对相同。
- **差异链的性质**：没有"块完全相同、只有链分不同"的情况（n_scorediff = 0）。
  - 差异链里有正链上的链，而正链先于负链处理，不受删链经由掩蔽的影响，所以差异来自串链细节本身，不是输入不同。
  - 细看了一条：Mouse-Rat 正链主链，两边端点相同、块数都是 64856。差别只在 3 处相邻块的切分（交叉点位置不同，例如 kent 的 `129 0 3 / 3 0 137` 对 lastz 的 `120 0 137 / 12 0 3`）。lastz 的链分高 4833（279,358,693 对 279,353,860）。
  - 这与 9.3 第 3 条的预期一致：kent 的剪枝下界不严格，偶尔会漏掉最优前驱。但只细看了这一条，其余差异链的原因没有逐条核实。

### 10.3 第 8 步：两套配置评测（全部达标）

运行：`nlr_default`（`--nlchain`）、`nlr_sec`（`--nlchain --nlchainsecondary=8800,1,0.9`），各 20 个作业；评测用 `eval_tune.sh`（mafComparator `--samples=10000000 --seed=1`，两份真值）；汇总用 `tune_report.py`，结果在 `results/nlr/step8_report.txt`。运行前 `results/tune/long.tsv` 备份为 `logs/long.tsv.bak_nlr`。

| 组 | R_np | P_np | F_np | R_anc | P_anc | F_anc | 模拟 F_np / F_anc | CPU (s) |
|---|---|---|---|---|---|---|---|---|
| nl_step5default（旧） | 0.7517 | 0.6707 | 0.7085 | 0.6681 | 0.7499 | 0.7054 | — | 1421 |
| **nlr_default** | 0.7492 | 0.7022 | **0.7247** | 0.6632 | 0.7819 | **0.7167** | 0.7247 / 0.7167 | 1678* |
| **nlr_sec** | 0.7494 | 0.6947 | **0.7208** | 0.6691 | 0.7803 | **0.7194** | 0.7208 / 0.7194 | 1702* |
| chainNet medium/3000 | 0.7505 | 0.6930 | 0.7205 | 0.6696 | 0.7780 | 0.7191 | — | 2884 |

- **F 与模拟值完全相同**（四个数的差都是 0.0000，验收线是 ≤ 0.002）。D6 的动态掩蔽效应在平均值上看不出来。
- **默认配置 F_np 0.7247 > 0.7205（chainNet）**，在 9/10 个物种对上不低于 chainNet，唯一低的是 Dog-Mouse（0.6221 对 0.6222）。
  - 默认配置的 F_anc 比 chainNet 低 0.0024，在 5/10 个物种对上不低于 chainNet。
  - 次级链配置的 F_np 0.7208、F_anc 0.7194，两项都略高于 chainNet，与模拟的结论相同。
- 碱基重复使用（tMul/qMul）从 0.046/0.038 降到 0.016/0.009。
- **新增阶段 CPU 0.21%（验收线 1%）**：`nlr_default_prof` 的 20 个作业中，`nlchain rechaining` 合计 3.08 s，`total run time` 合计 1436.9 s。计时用墙钟（`gettimeofday`），单线程作业下近似 CPU 时间。
- *CPU 一列的说明：`nlr_default`/`nlr_sec` 的 1678/1702 s 比旧的 1421 s 高约 18%，但原因是服务器负载，不是新代码。证据有两条：
  - `nlr_clean0` 跳过了新阶段、输出与旧版逐字节相同，与这两组同批 60 个作业并行时为 1712 s；
  - `nlr_default_prof` 单独以 20 个作业并行时为 1437 s。
  - 服务器当时的负载约为 60/256。

### 10.4 未解决 / 风险

- 第 7 步约 0.3% 的链与 axtChain 不同，只细看了一条，原因判断为交叉点或前驱的平局与剪枝差异，没有逐条核实。对碱基对和 F 的影响可以忽略（Jaccard ≥ 0.9997，F 与模拟相同）。
- 9.4 中其余的风险仍然存在：N 碱基计分、全基因组规模下的性能、`README.lastz.html` 和 `--help` 尚未补。
- 服务器上的新脚本和改动（`bootstrap_nlr.sh`、`nlr_vs_axtchain.py`、`run_tune.sh` 加的一行、`tune_groups.tsv` 加的 5 行）只在服务器上，不在本仓库。本仓库的脚本副本在本次会话的 scratchpad 中，会话结束后就没有了。

### 10.5 下一步（新会话）

1. 提交本文档的更新（先给用户看 diff 和 commit message）。
2. 方案 §4 第 9 步（可选）：在 A 之上重扫 T_sec ∈ {3h, 4h, 6h}，N_sec=1，r_sec=0.9。用 `TUNE_BIN=bin/lastz_nlr`，参数组写进 `tune_groups.tsv`，评测流程同 10.3。
3. 方案 §4 第 10 步：更新 `docs/nlchain_design.md`（§7 回溯流程加重新串链一步、§8 接口加 `--nlchainclean`/`--nlchaincleanreport`/`--nlchaingapcost`），补写 `README.lastz.html` 和 `lastz --help`。
4. 16 个最终测试物种对（§7 第 5 项）：由用户决定何时跑。服务器脚本需要先参数化，参数按 §8 加方案 A，二进制用 `bin/lastz_nlr`。

## 11. gapped 之前插值：本地实现与验证（2026-11-03）

实施方案：`docs/nlchain_prechain_inner_plan.md`（下文"方案 §x"指该文档）。完成了方案 §7 的第 1–4 步，全部在本地 OrbStack ubuntu VM（gcc 11.4）上完成。代码尚未提交；第 5 步（远端实验，方案 §8）未做。

### 11.1 已完成

| 文件 | 内容 |
|---|---|
| `src/tweener.h` / `src/tweener.c` | 新增 `tweener_collect_nl_hsps`：遍历第一次串链的主路径，在段内相邻 HSP 之间的窗口里用 7-mer 种子正反两个方向找 HSP，放进调用者给的两张临时表；无损预筛；`tweenerStats` 增加 `numNlWindows`、`numNlWindowsSkipped`、`nlWindowArea`、`numNlHsps[2]`、`numPass1Dropped`。`tweener_interpolate`、`bounded_align`、`collect_inner_hsps` 一行未动 |
| `src/lastz.h` / `src/lastz.c` | 参数 `nlInnerJump`、`nlInnerSolo`（默认 -1）；`--nlchaininner=<J2>[,<T_solo2>]`；D8 与"`--nlchaininner` 要求 `--inner`"两项检查；`nlchain_both_strands` 按方案 §3 第 3–5 步改为两次串链；私有函数 `append_inner_hsps`（D9）；`collect_stats` 下的 numPass1Dropped 统计；nlchain 下 `finish_one_strand` 不再调用 `tweener_interpolate` |
| `test_data/nlchain_inner/make_inner_test.py`（新） | 构造测试数据生成脚本（方案 §7 第 3 步），按 x-drop 校准 X、W 的 HSP 分数，按方案 §2.1 的公式核对两档跳转代价，写出 `expected.json` |
| `test_data/nlchain_inner/check_inner_test.py`（新） | 运行并判定构造测试的全部预期，可选对照 d9db595 二进制和 `collect_stats` 构建 |

`src/nlchain.c`/`nlchain.h` 未改。决策按方案 §4 的推荐执行，D6 取 §4.1 的 b2（J2 = min(J, ⌊innerThreshold/2⌋)，T_solo2 = min(T_solo, ⌊5·innerThreshold/4⌋)，对第二次串链的全部 HSP 生效）。

用法（构造测试）：

```
python3 test_data/nlchain_inner/make_inner_test.py <数据目录>
python3 test_data/nlchain_inner/check_inner_test.py <新 lastz> <数据目录> [--base <d9db595 的 lastz>] [--stats <collect_stats 构建>]
```

### 11.2 验证结果

- **构建**：gcc 11.4 加 `-Werror` 构建 `lastz`、`lastz_D`、`nlchain_test`、`nlrechain_test`，0 警告。`-Dcollect_stats -DdbgTiming`（不加 `-Werror`）构建通过，5 条警告都在原有文件里（`gapped_extend.c` 5824/5830 行 2 条 `-Wmisleading-indentation`、`pos_table.c` 2 条、`quantum.c` 1 条），新改的文件没有警告。
- **测试**：`make base_tests` 19 项 SUCCESS；`nlchain_test` 518,512 项检查 0 失败；`nlrechain_test` 2,929,359 项检查 0 失败。
- **回归**（对照 d9db595 构建的二进制，两边用相同的相对路径运行，比较 stdout、stderr、退出码）：
  - 不加 `--nlchain` 的 24 组全部逐字节一致。数据为 pseudocat vs pseudopig，22 组为：lav（默认）、maf、axt、general、sam、cigar、`--nogapped --format=segments`、`--chain`、`--chain=20,1`、`--nogapped`、`--hspthresh=top10%`、`--queryhspbest=5`、`--queryhsplimit=10`、`--strand=plus`、`--strand=minus`、`--anyornone`、`--inner=2000`、target `[multiple]`（pseudopig 三条）、query `[multiple]`、`--masking=3`、`--identity=70 --coverage=5`、UCSC 远缘参数（HoxD55、`--hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=4000`）。另加 2 组：aglobin `--inner=2000`、aglobin `--inner=1500 --chain`。之前会话的回归脚本已经丢失，这份清单是按 §9.2 的描述重建的。
  - `--nlchain` 的 8 组全部逐字节一致，含 `--nlchainreport` 段表：默认、maf、次级链 `8800,1,0.9`（aglobin）、`--hspthresh=top50%`、UCSC 远缘参数（aglobin，不加 `--inner`）、target `[multiple]`、`--nlchaingapcost=loose`、query `[multiple]`（替换原来的 `--inner` 组，方案 §2.6）。
  - `--nlchain --inner=2000 --nogapped`：两版都报 "--inner requires --gapped"，逐字节一致（见 11.3 第 2 条）。
- **互斥**：`--nlchain --inner --self` 报 "--nlchain can't be used with both --inner and --self"；`--nlchaininner` 不带 `--inner`、取负值（J2 或 T_solo2）、给 3 个值都报错；不带 `--nlchain` 时报 "require --nlchain"。`--self --nlchain`（不带 `--inner`）和非 nlchain 的 `--self --inner` 照常运行。
- **构造测试**（数据种子 20261103）：X 长 68 bp，s_X = 6476；W 长 23 bp，s_W = 2246。W→X→B 两次跳转的代价按方案 §2.1 的公式算，J2 = 1000 时为 3984，J = 3000 时为 7984，s_X 落在两者之间。`check_inner_test.py` 21 项全部通过：
  - 推导值：段表依次为 `S + A…W`（2 个 HSP）、`I − X`、`I + B`，坐标与构造一致；输出中有 X 的反链比对，没有 W；加 `--gappedthresh=2000` 后 W 出现。
  - `--nlchaininner=3000,3750`（D6a）：只有一段 `S + A…B`（3 个 HSP，含 W），没有 X，也没有 W 的比对。
  - `--nlchaininner=1000`：与推导值逐字节一致。
  - d9db595（同样参数，去掉 `--nlchaininner`）：没有 X，W 由旧 tweener 输出。
  - 变体：query 整体反向互补后，段表的 σ 对调，段表和比对坐标都与正向运行一致（这一组让正链方向搜索找 X、反链方向搜索找 W，两条坐标换算都覆盖到了）；target `[multiple]`、query `[multiple]`（构造都放在第二条序列上）与单条运行逐字节一致。
  - 自适应 `--hspthresh=top50%`：W 在主路径第一段的末尾；加 `--inner` 的输出只比不加多出 X。
  - `collect_stats` 构建：窗口 1 个，预筛跳过 0 个，插值 HSP 反/正 = 1/1，numPass1Dropped = 0，nlchain HSPs = 4（第一次 2 个 + 插值 2 个，说明 D10 的快照恢复生效）。
- **D9 临时调试**（只在 VM 的源码副本里加打印，仓库代码不含）：在 `finish_one_strand` 打印 `lowAnchorScore`，在 `append_inner_hsps` 前后打印 `lowScore`。构造数据加 `--hspthresh=top50%`：
  - 带恢复：加不加 `--inner`，`lowAnchorScore` 都是 836（top40% 时都是 974，top30% 时都是 1157）。
  - 去掉恢复的对照构建：反链表在压缩后为空，加入 X 后 `lowScore` 变成 6476，两表取小后 `lowAnchorScore` 从 836 变成 842。说明不做 D9 时门槛确实会变。
- **冒烟**（`collect_stats` 构建，`--nlchain --inner=2000`，加 `--debug=anchors` 看每次的 HSP 数）：

  | 数据 | 第一次输出 | 插值 HSP 反/正 | merge 后 = 第二次输入 | nlchain HSPs | 窗口 / 预筛跳过 | numPass1Dropped |
  |---|---|---|---|---|---|---|
  | pseudocat vs pseudopig（3 条 query 合计） | 17 | 22 / 0 | 39 | 39 | 11 / 0 | 0 |
  | 同上，次级链 `8800,1,0.9` | 17 | 22 / 0 | 39 | 39 | 11 / 0 | 0 |
  | aglobin human vs cow | 38 | 0 / 40 | 78 | 78 | 29 / 2 | **0** |
  | 同上，次级链 `8800,1,0.9` | 45 | 0 / 40 | 85 | 85 | 29 / 2 | 3 |

  - 计数自洽：第二次输入 = 第一次输出 + 插值 HSP（这几组 merge 都没有合并掉 HSP），也等于 `nlchainStats` 的 HSP 数。pseudopig 的比对都在反链上，所以插值 HSP 都是反向的。
  - 默认参数下 numPass1Dropped 为 0，符合预期。次级链下为 3，原因见 11.4。
- **预筛无损**：用关掉预筛的临时构建对照，3 组数据（aglobin、pseudocat、构造数据）× `--inner` 取 2000/5000/8000 × 开不开次级链，共 18 组，MAF（去掉回显命令行的注释行）、段表、stderr 全部一致。aglobin 上被跳过的窗口分别为 2、5、6 个。

### 11.3 与实施方案的差异（均已确认或有理由）

1. **numPass1Dropped 放在 `tweenerStats`**，由 `tweener_show_stats`/`tweener_generic_stats` 输出（方案原写 `nlchainStats`，但那要改 `nlchain.h/c`，与"不改 nlchain.c/h"冲突；用户确认）。新增的这几项统计只在有插值窗口时输出，非 nlchain 的 `--inner --stats` 输出不变。
2. **D7**：`--inner` 与 `--nogapped` 在原版就互斥（`lastz.c` 参数检查报 "--inner requires --gapped"），所以方案 §7 第 2 步"`--nlchain --inner --nogapped` 正常运行"不成立。保留原有报错，`nlchain_both_strands` 中按 `gappedExtend` 判断只作防御（用户确认）。
3. **W 的位置**：8 kb 空隙中 X 居中时，左右两侧各约 3.97 kb，W 不可能距 X 和两端都 ≥ 2 kb。实际放在距左端约 1.95 kb、距 X 约 2 kb 处；这个距离的用意是不让 A 和 X 的 gapped 延伸吞掉 W，默认 ydrop 下足够。W 在 X 之前，所以段表的第一段是 `S + A…W`，不是方案写的 `S + A`。
4. **自适应变体的预期改了**。方案原写"验证 W 的比对不会因为门槛被拉低而输出"，这在构造数据上做不到：自适应下主搜索本身就以很低的门槛（约 836）收集 HSP，W 不加 `--inner` 也会被找到并输出（d9db595 也一样）。另外显式给 `--gappedthresh` 时 gapped 门槛不随自适应变化，测不出 D9。改为不给 `--gappedthresh`，检查"加 `--inner` 只多出 X"；D9 本身由上面的临时调试证实。
5. **预筛的 maxSub 取整个替换矩阵的最大值**，不只看 ACGT（方案 §2.8）。这样对 N、IUPAC 等字符也是无损的；HoxD70 下两者都是 100。
6. **窗口上限不取偶**：旧 tweener 把 windowSize 向上取偶（`round_up_2`），新函数直接用 `--inner` 的窗口值。
7. 并入主表时沿用插值临时表里已分配的 hspId，不再取新号（方案只要求进临时表时传 0）。
8. `--debug=anchors` 的那一行比方案多打印 merge 后的 HSP 数，用于核对冒烟计数。

### 11.4 未解决 / 风险

- **次级链下 numPass1Dropped = 3**（aglobin，`--nlchainsecondary=8800,1,0.9`）：被丢的 3 个 HSP 正好是第一次次级链 2 的全部 HSP（cow 20988–22009 对 human 33112–35156，链分 32395，高于 T_sec）。用 `--nlchaininner=3000,3750`（第二次参数与第一次相同）时这条次级链同样消失，所以原因不是 J2 降低，而是插值 HSP 加强了这段 query 区间内的主路径（第一段 HSP 从 12 个增加到 31 个），次级链没能通过 r_sec 的区域比值检查。这属于方案 §4.1 预料到的"打开次级链时"的例外；默认参数（次级链关闭）下为 0。§8 的 G2p 要在 simTest 上再看这个数。
- **b2 下主路径被切成更多 G 段**：aglobin 上第一次 4 段，第二次 15 段（多出来的都是 G）。原因是 J2 = 1000 后，g(d) 超过约 1 kb 级跳转代价的正常边改走跳转。这与方案 §4.1"只改变段的划分"一致，不影响进入 gapped 的 HSP 集合，但 §8 用段表比较 I 段数时要只数 I，不能看总段数。
- 性能只在小数据上看过（插值计时 0.001–0.002 s）；MR_chr0 规模下窗口数约为旧版 10 倍（方案 §2.8），要在 G2p 看实际耗时。
- 回归脚本与选项清单这次也放在 scratchpad 中，会话结束后会丢失；清单已写在 11.2。构造测试的两个脚本已放进 `test_data/nlchain_inner/`。
- `README.lastz.html` 和 `lastz --help` 仍未写入任何 nlchain 选项（含 `--nlchaininner`）。

### 11.5 下一步（新会话）：远端实验（方案 §8）

1. 先提交本次实现（给用户看文件列表、diff 摘要和 commit message）。
2. 服务器构建：仿照 `scripts/bootstrap_nlr.sh`，从提交后的版本 `git archive` 构建正式版（gcc 13，`-Werror`）和 `-Dcollect_stats -DdbgTiming` 版，sha256 记入 MANIFEST；先跑 `make base_tests`、`nlchain_test`、`nlrechain_test`。远端坑见记忆文件（tar 关 xattr、ssh 后台启动包子 shell）。
3. 在 simTest 20 个作业上跑方案 §8 的各组：G0（新，`--nlchain`，应与 `nlr_default` 逐字节一致，只差命令行回显）、G1（d9db595，`--nlchain --inner=2000`；d9db595 与 8e9f350 的 `src/` 相同，可直接用服务器上的 `bin/lastz_nlr`）、G2（新，`--inner=2000`，b2）、G2b1（`--nlchaininner=2000,2500`）、G2a（`--nlchaininner=2200,2750`）、G3（`--inner=2000 --gappedthresh=2000`）、G3b（只加 `--gappedthresh=2000`）、G2p（统计版，同 G2），可选 G4（方案 B 补丁二进制，不提交）。参数组写进 `config/tune_groups.tsv`，评测和汇总同 §10.3。
4. 比较 R/P/F（np、anc）、tMul/qMul、CPU、插值与 nlchaining 耗时；插值窗口数、预筛比例、插值 HSP 数；`--nlchainreport` 中的 I 段数（G2/G2b1/G2a 对 G0）；G2p 的 numPass1Dropped。按方案 §8 的判定规则决定 D6 的默认值、是否把 gapped 门槛默认降到 innerThreshold、是否改用方案 B。

## 12. gapped 之前插值：远端实验（方案 §8，2026-11-03）

服务器与目录同 §10（`PROJECT_ROOT`）。实现已在 19fa9e5 提交。汇总：`results/nlr/nli_report.txt`（R/P/F）、`results/nlr/nli_stats.txt`（段表与 G2p 计数）；补跑的 G3a 见 12.5，默认值的决定与验证见 12.6，汇总在 `results/nlr/nli_report_g3a.txt`。

### 12.1 构建与测试

- 源码：`git archive 19fa9e5` 解到 `external/lastz-19fa9e5`（上传后核对 md5）。构建脚本 `scripts/bootstrap_nli.sh`，由 `bootstrap_nlr.sh` 改来，sha256 记在 `bin/MANIFEST_nli.tsv`：
  - `bin/lastz_nli`：正式版，`-Werror`，gcc 13.1，sha256 `522ac491…a4fd95`。
  - `bin/lastz_nli_prof`：`-Dcollect_stats -DdbgTiming`，sha256 `0862a719…f0a5d0`。5 条警告都在原有文件里（`gapped_extend.c` 2、`pos_table.c` 2、`quantum.c` 1），与本地相同；脚本检查 `nl*.c`、`tweener.c`、`lastz.c` 没有警告。
- 测试（`build/nli_tests/`）：`nlrechain_test` 2,929,359 项检查 0 失败；`nlchain_test` 518,512 项检查 0 失败；`base_tests` 19 项 SUCCESS；`lastz_D` 加 `-Werror` 通过；构建日志无 warning。
- 冒烟（pseudocat vs pseudopig，simTest 参数）：两版在 `--nlchain` 和 `--nlchain --inner=2000` 下输出与清洗报告都一致。`--inner` 下窗口 33 个，预筛跳过 3 个，插值 HSP 反/正 = 32/1，numPass1Dropped = 1（本地默认参数下为 0，参数不同）。

### 12.2 运行

- 参数组 `nli_g0`、`nli_g1`、`nli_g2`、`nli_g2b1`、`nli_g2a`、`nli_g3`、`nli_g3b`、`nli_g2p` 加在 `config/tune_groups.tsv` 末尾，原文件备份在 `logs/tune_groups.tsv.bak_nli`。G3/G3b 的 `--gappedthresh=2000` 写在公共参数的 `--gappedthresh=6000` 之后，lastz 取后一个。G2p 加 `--stats`（不带文件名时写 stderr，与计时在一起）。
- 运行脚本 `logs/nli_run.sh`、`logs/nli_run2.sh`：G0、G2、G2b1、G2a、G3、G3b 共 120 个作业并行；G2p 单独以 20 个作业并行；评测同 §10.3。`results/tune/long.tsv` 运行前备份为 `logs/long.tsv.bak_nli`。
- **G1 没有结果**：`bin/lastz_nlr`（8e9f350，与 d9db595 的 `src/` 相同）加 `--nlchain --inner=2000`，20 个作业全部以 `FAILURE: outer alignments out of order` 退出。原因是 `nl_rechain_aligns` 之后比对表不再按 beg1 排序，而 `tweener.c` 默认定义了 `verifyOrDie`，`tweener_interpolate` 入口检查到乱序就退出。重新串链之前的 `lastz_nl_step5`、`lastz_nl` 跑同一个作业正常。这是 8e9f350 引入、19fa9e5 已绕开的问题（nlchain 下不再调用 `tweener_interpolate`），§9.2 的小数据回归没有触发。用户决定不与旧版比较，G1 不再补跑；判定规则中"高于 G1"一条去掉。G4 按用户决定不做。
- **G0 = `nlr_default`**：20 个作业的 MAF（去掉 `#` 行）和段表逐字节一致。**G2p = G2**：同样逐字节一致。

### 12.3 准确度

| 组 | 参数（公共参数之外） | R_np | P_np | F_np | R_anc | P_anc | F_anc | tMul / qMul | MAF 块数 |
|---|---|---|---|---|---|---|---|---|---|
| G0 `nli_g0` | `--nlchain` | 0.7492 | 0.7022 | 0.7247 | 0.6632 | 0.7819 | 0.7167 | 0.016 / 0.009 | 54,860 |
| G2 `nli_g2` | `--inner=2000`（b2：J2=1000，T_solo2=2500） | 0.7669 | 0.7040 | 0.7340 | 0.6788 | 0.7842 | 0.7272 | 0.018 / 0.010 | 60,218 |
| G2b1 | `--inner=2000 --nlchaininner=2000,2500` | 0.7663 | 0.7042 | 0.7338 | 0.6783 | 0.7843 | 0.7269 | 0.018 / 0.010 | — |
| G2a | `--inner=2000 --nlchaininner=2200,2750` | 0.7661 | 0.7042 | 0.7338 | 0.6781 | 0.7843 | 0.7268 | 0.018 / 0.010 | — |
| **G3** | `--inner=2000 --gappedthresh=2000`（b2） | 0.7773 | 0.7110 | **0.7426** | 0.6882 | 0.7922 | **0.7360** | 0.020 / 0.012 | 78,091 |
| G3b | `--gappedthresh=2000` | 0.7536 | 0.7053 | 0.7284 | 0.6672 | 0.7854 | 0.7206 | 0.017 / 0.009 | 62,597 |
| chainNet medium/3000 | — | 0.7505 | 0.6930 | 0.7205 | 0.6696 | 0.7780 | 0.7191 | 0.000 / 0.016 | — |

- **G2 有效**：F_np +0.0093、F_anc +0.0105（对 G0），10/10 个物种对两项都不低于 G0。收益主要来自 recall（R_np +0.018），P_np 也略升（+0.0018），没有出现 §9 担心的 precision 下降。远缘物种对收益更大（Cow-Rat F_np +0.021，Mouse-Rat +0.0003）。
- **D6 三档几乎相同**：G2、G2b1、G2a 的 F_np 相差 0.0002，F_anc 相差 0.0004，都在判定线 0.001 之内。按方案 §8 的规则取 a（第二次与外层同用 J、T_solo）。b2 两项略高，但差值远小于判定线。
- **降低 gapped 门槛收益更大**：G3 对 G2 再加 F_np +0.0086、F_anc +0.0088，P 也升（P_np +0.007）。G3b 单独降门槛只有 +0.0037/+0.0039，所以 G3 的增量主要来自插值 HSP 不再被 6000 的 gapped 门槛挡掉（§2.2 的 D5a 问题）。
- G3 在 10/10 个物种对上 F_np、F_anc 都高于 chainNet medium/3000。

### 12.4 段表、计数与耗时

主路径段数（10 对合计，`nli_stats.py`）：

| 组 | D | G | I | I 段 query bp |
|---|---|---|---|---|
| G0 | 2,764 | 2,520 | 8,173 | 43.6 M |
| G2 / G3 | 8,785 | 10,360 | 12,837 | 25.3 M |
| G2b1 | 3,595 | 2,816 | 11,003 | 44.8 M |
| G2a | 3,277 | 2,579 | 10,491 | 46.2 M |

- G3 与 G2、G3b 与 G0 的段表统计相同（gapped 门槛不影响串链）。
- 插值让 I 段增加 28%（G2a）到 57%（G2）。b2 下 J2=1000 把主路径切得更碎（G 段 ×4），I 段覆盖的碱基反而减少，所以 G2 的 I 段数不能与其他组直接比较（§11.4 已预见）。

G2p（20 个作业合计，对照 `nlr_default_prof`，§10.3）：

| 项 | nlr_default_prof | G2p |
|---|---|---|
| 插值窗口 / 预筛跳过 | — | 221,015 / 53,792（19.6%） |
| 插值 HSP 反 / 正 | — | 79,229 / 117,072 |
| numPass1Dropped | — | 43 |
| gapped extensions | 100,926 | 129,178（+28%） |
| interpolation (s) | 0.0 | 3.0（占 total 0.18%） |
| nlchaining (s) | 4.3 | 10.9（两次串链，占 0.66%） |
| nlchain rechaining (s) | 3.1 | 3.5 |
| gapped extension (s) | 223.3 | 277.7 |
| seed hit search (s) | 1178.4 | 1320.6 |
| total run time (s) | 1436.9 | 1647.8 |

- 插值和两次串链本身合计约 10 s，不到 1%。主要的新增开销是 gapped 延伸：次数 +28%，约 +55 s（约 4%）。
- total 多出的 211 s 中有 142 s 在 seed hit search，这一段代码没有改，是服务器负载差异（两次运行时负载不同）。CPU 一列（G0 1775 s 对 `nlr_default` 1678 s，二者输出相同）也有同样的漂移，所以组间 CPU 不可直接比较。
- 预筛跳过比例随物种距离下降：Mouse-Rat 47%，Cow-Rat 3.6%。
- numPass1Dropped 合计 43（每对 1–9），次级链关闭。少量第一次主路径上的 HSP 在第二次串链中没有留在主路径上，与 §9"第二次串链改变主路径"一致；数量相对 120 万个 HSP 可忽略，没有逐个核实。
- G3 没有统计版，gapped 延伸增加多少没有测量（MAF 块数 78k 对 G2 60k）。

### 12.5 补跑：G3a 与降门槛的开销

G3 用的是 b2，而 D6 按规则应取 a，所以补跑三组（`logs/nli_run3.sh`）：

- `nli_g3a`：`--inner=2000 --gappedthresh=2000 --nlchaininner=2200,2750`，正式版。
- `nli_g3ap`：同 G3a，统计版；与 G3a 的 MAF（去掉 `#` 行）逐字节一致。
- `nli_g2ap`：同 G2a，统计版。与 G3ap 同批各 20 个作业并行，两组只差 gapped 门槛，负载相同，耗时可以直接比较。

| 组 | R_np | P_np | F_np | R_anc | P_anc | F_anc | tMul / qMul | MAF 块数 |
|---|---|---|---|---|---|---|---|---|
| G3（b2） | 0.7773 | 0.7110 | 0.7426 | 0.6882 | 0.7922 | 0.7360 | 0.020 / 0.012 | 78,091 |
| **G3a**（a） | 0.7761 | 0.7113 | **0.7422** | 0.6870 | 0.7925 | **0.7355** | 0.019 / 0.011 | 77,364 |

- G3a 与 G3 相差 F_np 0.0004、F_anc 0.0005，在判定线 0.001 之内。D6 取 a 不损失降门槛的收益。
- G3a 在 10/10 个物种对上 F_np、F_anc 都高于 G0、G2a 和 chainNet medium/3000。

| 项（20 个作业合计） | G2ap（门槛 6000） | G3ap（门槛 2000） |
|---|---|---|
| gapped extensions | 126,170 | 123,606 |
| gapped extension (s) | 259.0 | 256.8 |
| interpolation (s) | 2.9 | 2.9 |
| nlchaining (s) | 10.3 | 10.4 |
| nlchain rechaining (s) | 3.2 | 3.4 |
| total run time (s) | 1436.7 | 1448.8 |
| numPass1Dropped | 55 | 55 |

- **降门槛没有可见的开销**：gapped 延伸次数反而少 2%，total 相差 0.8%（seed hit search 也差 13 s，属于噪声）。gapped 门槛只决定延伸后的比对留不留，不决定做不做延伸；次数变少的原因推测是留下的比对更多，后续落在其中的 HSP 被跳过，没有核实。
- D6 取 a 时 numPass1Dropped 为 55（b2 下为 43），说明第一次主路径 HSP 在第二次被丢掉与 J2 的取值无关，来自第二次串链改走了别的路径。
- G2ap 的计数与 G2p 相同（窗口、插值 HSP 只取决于第一次串链和 `--inner`）。

### 12.6 定默认值（用户决定，2026-11-03）

- **D6 取 a**：未给 `--nlchaininner` 时，第二次串链沿用第一次的 J、T_solo（`src/lastz.c` `nlchain_both_strands`，删掉 min(J, innerThreshold/2)、min(T_solo, 1.25·innerThreshold) 的推导；`lastz.h` 注释同步）。`--nlchaininner` 保留，只给 J2 时 T_solo2 沿用外层的 T_solo。
- **gapped 门槛：代码不改（G-1）**。simTest、UCSC 式参数都显式给 `--gappedthresh`，而显式参数优先，"未给时默认降到 innerThreshold"对这些流程不起作用。改为在 `nlchain_design.md` §8.3 推荐 `--inner` 时同时给 `--gappedthresh=<innerThreshold>`。没有显式给 gapped 门槛的配置（G-2）没有实验数据，不改。
- 构造测试数据不变（X 按 b2 档校准）。`check_inner_test.py` 中依赖 X 进入主路径的各项改为显式传 `--nlchaininner=1000,2500`；"默认"一项改为检查 a 的行为（只有一段 `S + A…B`，没有 X、W），并与 `--nlchaininner=3000,3750` 逐字节一致；"只给 J2"一项改为与 `--nlchaininner=1000,3750` 逐字节一致。
- 文档：`nlchain_design.md` 第 5′ 步、§8.2 第 10 条、§8.3（表格与 gapped 门槛建议）、修订记录；`nlchain_prechain_inner_plan.md` 开头与 §4.1 加结论。

验证：

- 本地 VM（gcc 11.4）：`lastz`、`lastz_D` 加 `-Werror` 0 警告；`collect_stats` 构建新改文件无警告；`base_tests` 19 项 SUCCESS；`nlchain_test` 518,512 项、`nlrechain_test` 2,929,359 项检查，均 0 失败。
- 构造测试 20 项全部通过（含 `collect_stats` 计数 5 项）。
- 回归对照 19fa9e5（去掉 `#` 行的 stdout、段表、stderr、退出码），18 项全部通过：
  - 不经过改动的路径 5 组一致：非 nlchain 的默认 lav、`--inner=2000`、`--inner=2000 --chain`，以及 nlchain 不加 `--inner`（含次级链）。
  - 新版显式给 b2 值 = 旧版默认，5 组：aglobin `--inner=2000`/次级链/`--inner=5000`（`--nlchaininner=2500,3750`），pseudocat 默认参数与远缘参数（`--hspthresh=2200 --gappedthresh=6000 --ydrop=3400`）。
  - 新版默认 = 旧版显式给 a 值，7 组：同上 5 组，加远缘参数下 `--gappedthresh=2000`、显式 `--nlchainjump=2000 --nlchainsolo=2600`。
  - 反例：aglobin `--inner=2000` 新默认段表 4 段，旧默认 15 段，确认改动生效。
  - 回归脚本在会话 scratchpad 中，会话结束后丢失，清单如上。
- 远端：19fa9e5 源码副本只替换 `lastz.c`、`lastz.h`（上传后核对 md5），编出 `bin/lastz_nli_d6a`（`-Werror`，0 警告，sha256 `7d811a08…9327c0f`，未写入 MANIFEST）。参数组 `nli_d6a_g3`（`--nlchain --inner=2000 --gappedthresh=2000`，不给 `--nlchaininner`）20 个作业的 MAF（去掉 `#` 行）和段表与 `nli_g3a` 逐字节一致，stderr 为空。

### 12.7 未解决 / 风险

- G3b 显示不加 `--inner`、只把 gapped 门槛降到 2000 也有 +0.004 的收益。这超出本方案范围，未处理。
- 8e9f350/d9db595 的 `--nlchain --inner` 在真实数据上崩溃（12.2），19fa9e5 已不走这条路径，旧提交不修。
- numPass1Dropped 在 simTest 上为 43–55（12.4、12.5），原因推断为第二次串链改走别的路径，没有逐个核实。
- 服务器上的新脚本（`bootstrap_nli.sh`、`nli_stats.py`、`logs/nli_run*.sh`）、`tune_groups.tsv` 加的 12 行和 `bin/lastz_nli_d6a` 只在服务器上；`bin/lastz_nli*` 是 19fa9e5 的构建，不含本次默认值改动。
- `README.lastz.html` 和 `lastz --help` 仍未写入任何 nlchain 选项。

### 12.8 下一步（新会话）

1. 补写 `README.lastz.html` 和 `lastz --help` 中的 nlchain 选项（含 `--nlchaininner` 和"`--inner` 时同时降 `--gappedthresh`"的建议）。
2. 16 个最终测试物种对（§7 第 5 项）：由用户决定何时跑。服务器脚本需要先参数化；参数按 §8 加方案 A，`--inner=2000 --gappedthresh=2000`，二进制从本次提交构建。
3. 可选：单独评估不加 `--inner` 时降低 gapped 门槛的收益（G3b）。
