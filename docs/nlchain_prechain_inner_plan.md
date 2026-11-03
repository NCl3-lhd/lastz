# nlchain 模式下 gapped 之前做正反两链插值：开发文档

日期：2026-11-03。状态：方案已核实（§2）；用户已同意按本文档实施（2026-11-03），§4 的 D1–D12 按"推荐"列执行，D6 按 §4.1（暂定默认 b2）；未动代码。

**2026-11-03 更新**：§7 第 1–4 步已在本地完成并验证，记录见 `docs/nlchain_progress.md` §11；第 5 步（§8 远端实验）也已完成，见该文档 §12。按实验结果，D6 默认改为 a（见 §4.1 末尾的结论），D5 维持 a（代码不改），改为在文档中推荐 `--gappedthresh` 与 `--inner` 取同值。

本方案取代 `docs/nlchain_inner_plan.md`（"窗口内改用 nlchain 串链"，未实施）。两者都改 `--nlchain --inner` 的行为，只能选一个。

## 1. 目标与新流程

```
主种子搜索 → 正反两张 HSP 表（split_anchors：反链 anchors，正链 secondaryAnchors）
  ↓
nlchain 第一次串链（正反一起）→ 主路径；两张表压缩为主路径 + 次级链 HSP
  ↓
【新增】主路径上每段内相邻 HSP 之间的窗口：7-mer 种子，正反两个方向都搜
  ↓
【新增】新 HSP 按方向并入两张表，各自 merge，nlchain 再串一次
  ↓
反链 finish_one_strand → 正链 finish_one_strand（gapped → 过滤 → 剥链 + T_chain → 输出）  ← 不改
```

要解决的问题：现在的插值（`tweener_interpolate`）在 `nl_rechain_aligns` 之后运行（`src/lastz.c:3534–3559`），插值出的比对不经过 T_chain；而且窗口只在当前链方向上搜，找不到空隙里的反向片段。

## 2. 核实结论：原方案中需要修正的地方

按影响从大到小排列。行号以 d9db595 为准。

### 2.1 "空隙里的小倒位能通过反向 HSP 进入主路径"：在当前 DP 下多数情况实现不了

设 A、B 是同一段内相邻的正链 HSP，X 是夹在空隙里的反链 HSP。X 进入主路径要走 A→X→B，两条都是跳转边（I）。按 `nlchain_design.md` §5：

- X 单独成段，要求 s_X ≥ T_solo（E4；只有 R 能跳出），否则 X 停在 U 状态，跳不出去。
- 还要比直连 A→B 好：`s_X − P_J(A,X) − P_J(X,B) > −C_N(A,B)`。
- 平衡倒位不改变两侧的对角线，所以 d(A,B) ≈ 0，C_N(A,B) = g(d) ≈ 0–900。
- 默认参数（J = 3000，λ = 1000，D0 = 3000）下，5 kb 空隙中间夹 400 bp 倒位：两次跳转的 D 都约 2.7 kb，P_J ≈ 3000 + 1000·log2(1 + 2.7/3) ≈ 3930。两次合计约 7.9k。

所以 X 要有大约 8k 分才进得了主路径。按 HoxD70 粗算，90% 一致度约 100 bp、70% 一致度约 160 bp 就有这个分数。这种片段主种子搜索（12of19，hspthresh 3000）通常已经能找到。7-mer 能多找到的，主要是分数在 2000 到 8000 之间的反向片段，DP 会把它们丢掉。

结论：
- 好处 1 改为"可能有，需要实测"。§8 的实验要统计 `--nlchainreport` 中 I 段的数量，单独量化这一项。
- 小倒位召回不到，根子在 DP 的跳转代价，与插值无关：第一次串链里，主种子找到的同类反向 HSP 也是这样被丢掉的。如果小倒位的召回重要，应该另开方案调整 I 型跳转的代价，不靠插值。
- D6 的理由"插值 HSP 靠和邻居连成链留下，不需要单独跳转"，和好处 1 互相矛盾：反向片段只能靠跳转进入主路径。D6 选 a，就等于默认放弃好处 1 中分数低于约 8k 的那部分。D6 改为第二次单独取值，见 §4.1。

### 2.2 D5a 会让 `--inner=2000` 实际变成 3000 门槛，和旧版没法直接比

- 旧版窗口内的 gapped 门槛是 innerThreshold（`tweener.c` 中 `scoreThresh.s = _scoreThresh`，传入的是 `currParams->innerThreshold`）。
- 新版的插值 HSP 走主 gapped，门槛是 gappedThreshold = hspthresh = 3000（`lastz.c:9695`）。
- 插值 HSP 只有三种结局：
  - 被外层 HSP 的 gapped 延伸吞掉：没有新增贡献；
  - 自己延伸出来但不到 3000：被丢掉；
  - 自己延伸到 3000 以上：留下。
- 因此 D5a 下，新版只能多找回"12of19 漏掉、7-mer 找得到、gapped 分 ≥ 3000"的片段。旧版能找回到 2000 分。§8 比较新旧两组时，差异里混着这个门槛变化。

不需要改代码就能把门槛的影响单独测出来：
- 外层 HSP 都满足 s ≥ hspthresh。它的两半都是非负分，从锚点出发的 gapped 延伸的分数不低于 s。
- 所以把 gapped 门槛降到 innerThreshold，几乎只影响由插值 HSP 延伸出来的比对。
- 实验中加一组 `--gappedthresh=2000` 即可（§8 的 G3）。另加 G3b（不加 `--inner`、只降门槛）确认上面的推论。

建议：代码按 D5a 做（不改），由 G3 的数据决定是否把"降门槛"做成默认。

### 2.3 自适应 hspthresh 下，往表里加 HSP 会改变 gapped 门槛，还可能弄坏表

- `add_segment` 每加一条都更新 `st->lowScore`（`segment.c:1048`）。
- `finish_one_strand` 在自适应模式下取 `anchors->lowScore` 和 `secondaryAnchors->lowScore` 的较小者作为 gapped 门槛（`lastz.c:3341–3348`、`3468–3473`）。插值 HSP 加进去以后，gapped 门槛会被拉低到 innerThreshold 附近。
- 自适应模式下，`anchors` 带 coverageLimit（`lastz.c:1442–1451`）。`compact_tables` 不更新 `coverage`/`lowScore`（`nlchain.c:1452–1489`），`coverage` 仍是压缩前的值。如果它已达到上限，`add_segment` 会：
  - 直接丢掉分数低于 lowScore 的新 HSP（`segment.c:1015–1019`）；
  - 或者把压缩后已不是堆的数组当作最小堆来维护。
- 处理办法：加 HSP 前保存 `coverageLimit`、`coverage`、`lowScore`，把 `coverageLimit` 置 0，加完后直接给字段赋值恢复。不能调 `limit_segment_table` 恢复，它会按上限裁表。

原方案没有提到这一点。之前 8 组 `--nlchain` 回归里有自适应阈值这一组，所以必须处理（D9）。

### 2.4 边遍历边加 HSP 会读到失效指针

`add_segment` 可能 `realloc` 整张表，返回新指针（`segment.c:1023–1033`）。新函数一边按 `nlsegment.first/count` 读主路径 HSP，一边往同一张表里加，会读到已释放的内存。
- 新 HSP 先收集到两张临时表里，遍历结束后再并入。
- 接口传 `segtable**`，因为 `anchors`/`secondaryAnchors` 的指针可能变化。

### 2.5 `--nogapped`、`--self`、`--queryhspbest` 的行为要定

- **`--nogapped`**：旧版 tweener 只处理 gapped 比对，`--inner` 在这种情况下什么都不做。新流程在 gapped 之前，如果照常运行，插值 HSP 会直接作为 HSP 输出。建议跳过新流程，与原版语义一致（D7）。
- **`--self`**：新窗口跟着主路径走。如果主路径在平凡对角线附近（串联重复区），窗口会跨过平凡对角线。原版 tweener 搜索时传的也是 `selfCompare=false`（`tweener.c:800`），但平凡比对要到 gapped 阶段才由 `inhibitTrivial` 挡掉。新流程里，平凡 HSP 会先进入 nlchain DP，可能主导主路径。建议在 nlchain 下禁止 `--inner` 与 `--self` 同用（D8）。
- **`--queryhspbest`**：插值 HSP 不受数量上限约束，与旧版一致，只在文档中说明。

### 2.6 回归集需要换掉一组

`nlchain_progress.md` §9.2 的 8 组 `--nlchain` 回归中有一组是 `--inner`，这一组按设计会变。换成 query `[multiple]`，凑足 8 组；另外单独对 `--nlchain --inner` 做 §7 第 3 步的构造测试。

### 2.7 D4 的理由要改，结论不变

- HSP 表很小：MR_chr0 上 15.7 万个 HSP，每个约 40 字节，共约 6 MB（`nlchain_profiling.md` §4）。"内存翻倍"不成立。
- 选 a 的真实理由：
  - 第一次的主路径在第二次中仍然可行（表里的 HSP 都还在），第二次的主路径分不会低于第一次；
  - 第二次 DP 的规模小（主路径 HSP 约为全部 HSP 的 28%–79%）；
  - 实现简单。

### 2.8 窗口数量和耗时

- MR_chr0 上主路径 HSP 有 44,602 个，gapped 延伸 3,640 次，新窗口数约为旧版的 10 倍以上。
- 每个窗口要建一次 7-mer 位置表（4^7 个桶，`zalloc`）。
- 大部分窗口位于 gapped 延伸本来就会跨过的小空隙里。

无损预筛：一个 HSP 的长度不超过窗口较短的一边，分数不超过"长度 × 最大单列分"。所以 `min(dq, dt) × maxSub < innerThreshold` 的窗口里不可能有达标的 HSP，可以直接跳过（maxSub 取 `maskedScoring` 替换矩阵 ACGT 部分的最大值，HoxD70 为 100）。另外，target 窗口的位置表正反两次搜索共用，只建一次。

### 2.9 其他小问题

- **"query 反向互补后再搜"写反了**：调用时 query 处于反向互补状态（`lastz.c:1685`、`1727`）。反链方向直接截取即可，正链方向要把截出的片段反向互补。坐标换算见 §6。
- **`nlchainStats` 会累加两次**：第二次调用前恢复第一次调用前的快照，使统计反映最终结果（D10）。
- **T_chain 的把关比原文说的弱**：剥链只删"自成一条、链分 < 5000"的链。插值比对如果能和外层比对串成一条链（即串进去链分更高），就会留下。也就是说，T_chain 去掉的是孤立的插值比对，不是所有低分的插值比对。
- **README 目前没有任何 nlchain 选项**（`nlchain_progress.md` §10.4）。单独给 `--inner` 补一句没有上下文，建议放进 README/`--help` 统一补写那一步（rechain 方案 §4 第 10 步），本方案只更新 `docs/`。
- **原文"复用 static 工具函数"成立**：`extract_subsequence`、`collect_inner_hsps` 和 `simpleInfo` 的设置都在 `tweener.c` 里，放在同一个文件里合理。
- **窗口都在同一个 query 分区、同一个 target 分区内**：nlchain 按 query 分区分别做 DP（`nlchain.c:1563–1571`），正常边要求 p 相同。所以新窗口不需要 `try_bounded_align` 那套跨分区拆分。

## 3. 修正后的流程细节

```
nlchain_both_strands(target, query):
  1. 解析 J、T_solo（不变）
  2. merge 两张表（不变）
  3. [collect_stats] 保存 nlchainStats 快照
  4. nl_reduce_to_chains(...) → result1                       （第一次）
  5. 若 innerThreshold > 0 且 gappedExtend：
       a. tweener_collect_nl_hsps(..., &result1, &inMinus, &inPlus)
       b. 把 inMinus 并入 anchors、inPlus 并入 secondaryAnchors（§2.3 保护）
       c. 若 mergeAnchors：两张表各自 merge_segments
       d. [collect_stats] 恢复 nlchainStats 快照
       e. nl2 = nl；nl2.jumpPen = J2；nl2.soloScore = T_solo2（§4.1）
       f. nl_free_result(result1)；nl_reduce_to_chains(..., &nl2, ...) → result2   （第二次）
       g. [collect_stats] 统计第二次丢掉的第一次主路径 HSP 数（§4.1 的验证）
  6. dbgShowAnchors 信息；nl_write_report(最终结果)；nl_free_result
finish_one_strand:
  - 插值调用加条件 (!currParams->nlchain)，nlchain 下不再走旧 tweener（D1a）
```

`nl2` 从第一次解析好的 `nl` 复制，只改 J、T_solo，不重新按 hspthresh 解析。如果重新解析，自适应模式下 J 会被插值 HSP 的低分拉低。

## 4. 设计决策

D1–D6 是原方案中的编号。

| # | 问题 | 推荐 | 说明 |
|---|---|---|---|
| D1 | nlchain 下旧的 gapped 后插值 | **a 关掉** | 同意原方案 |
| D2 | 窗口怎么划 | **a 段内相邻 HSP** | 同意。主路径段内相邻的两个 HSP 就是"同方向、走正常边"的一对。只用主路径（chainId = 0），次级链段不划窗口 |
| D3 | 主路径两端外延 | **a 不搜** | 同意 |
| D4 | 第二次串链的输入 | **a 压缩后的表 + 新 HSP** | 结论同意，理由改为 §2.7 |
| D5 | gapped 阈值 | **a 代码不改**，实验加 G3/G3b | 见 §2.2。由数据决定是否默认把 gapped 门槛降到 innerThreshold |
| D6 | 第二次串链的 J、T_solo | **第二次单独取值，全局生效**：J2 = min(J, ⌊innerThreshold/2⌋)，T_solo2 = min(T_solo, ⌊5·innerThreshold/4⌋)；新选项 `--nlchaininner` 可覆盖 | 用户意见：第二次的值可以调整。取值与理由见 §4.1；默认值是暂定的，由 §8 的扫描决定 |
| D7（新） | `--nogapped` 下的 `--inner` | 跳过新流程，不报错 | 与原版"`--inner` 只作用于 gapped"一致 |
| D8（新） | `--nlchain --inner --self` | 报错 | §2.5；备选：跳过与平凡对角线相交的窗口 |
| D9（新） | 自适应 hspthresh | 保存并恢复 `coverageLimit`/`coverage`/`lowScore` | §2.3；备选：nlchain 下禁止 `--inner` 与自适应阈值同用 |
| D10（新） | `nlchainStats` | 统计只反映第二次串链 | 另加插值计数（§5.4） |
| D11（新） | 无损预筛 | 做 | §2.8；只跳过不可能产出达标 HSP 的窗口，不改变结果 |
| D12（新） | 对照方案 B | 做成实验用的补丁二进制，不提交 | 见 §8 的 G4 |

### 4.1 D6：第二次串链的 J2、T_solo2

> **结论（2026-11-03，§8 实验后）**：默认改为 a（J2 = J、T_solo2 = T_solo）。simTest 上 b2、b1、a 三档 F 相差 ≤ 0.0004，在 §8 判定线 0.001 之内；降 gapped 门槛后（G3 对 G3a）同样只差 0.0005。下文是实验前的分析，保留作记录。

**全局生效就够了，不用只对插值 HSP 生效。** 按 D4a，第二次的输入只有第一次的主路径（和次级链）HSP 加上插值 HSP。
- 外层 HSP 本来就都在路径上，J 降低后，第二次对它们能做的只有一件事：把 g(d) 大于 P_J 的正常边改切成跳转（G 段）。这只改变段的划分，不改变进入 gapped 的 HSP 集合。
- 所以全局降低 J2，实际效果几乎只作用在插值 HSP 上。
- "只对插值 HSP 生效"要做的事情多得多：
  - 给 `nl_reduce_to_chains` 加逐 HSP 标记（公共接口）；
  - 改 P_J 和 T_solo 的取法；
  - 跳转池的剪枝界要改用 min(J, J2)；
  - 扩充 `nlchain_test` 的对拍。
  - 不值得。
- 例外：打开次级链（T_sec > 0）时，主路径可能多吸收一部分次级链 HSP，剩下的部分可能低于 T_sec 而被丢掉。
- 验证：`collect_stats` 构建中统计"第一次主路径上、第二次被丢掉的 HSP 数"，预期为 0 或接近 0（§3 第 5g 步）。

**取值的目标：让分数够过 T_chain 的倒位片段进得来。**
- 孤立的反向片段 X 在剥链里自成一条链（剥链按链方向分别做），链分约等于它的 gapped 分，低于 T_chain（5000）就会被删掉。
- 所以 X 进入主路径的门槛压到 T_chain 附近就够了，再低只会多做 gapped 再删掉。
- 门槛 ≈ 2·(J2 + λ·log2(1 + D/D0))，D 为两次跳转的跨度。λ = 1000、D0 = 3000 时：

| 空隙（X 居中） | 跨度项（每次） | J2 = J（D6a，h=2200 / 3000） | J2 = innerThreshold = 2000 | J2 = innerThreshold/2 = 1000 |
|---|---|---|---|---|
| 5 kb | ≈ 910 | 6.2k / 7.8k | 5.8k | 3.8k |
| 8 kb | ≈ 1250 | 6.9k / 8.5k | 6.5k | 4.5k |
| 10 kb | ≈ 1440 | 7.3k / 8.9k | 6.9k | 4.9k |
| 20 kb | ≈ 2130 | 8.7k / 10.3k | 8.3k | 6.3k |

**三档取值**：

| 选项 | J2 | T_solo2 | 说明 |
|---|---|---|---|
| a（原 D6a） | J | T_solo | 放弃 2k–8k 的倒位 |
| b1 | innerThreshold | 1.25 × innerThreshold | 与外层"J = hspthresh、T_solo = 1.25h"同形。但 simTest 用的是 UCSC 远缘参数（h = 2200），J 只从 2200 降到 2000，门槛几乎不变 |
| **b2（推荐，暂定默认）** | **innerThreshold / 2** | **1.25 × innerThreshold** | 10 kb 以内的空隙，门槛落在 T_chain 以下（3.8k–4.9k），所以够过 T_chain 的倒位片段都能进来 |

- **"第二次不比第一次严"**：两者都与外层取较小值（min）。用户显式给的 `--nlchainjump`/`--nlchainsolo` 比推导值还小时，第二次跟着用小的。
- **λ、D0、L 不改**：它们是绝对值，已在 simTest 上调定；窗口不超过 20 kb，跨度项有上界（约 2.1k）。如果 §8 显示 I 段召回受跨度项限制，再考虑单独给第二次设 λ。
- **T_solo2 不跟着 J2 减半**：T_solo2 要明显高于 innerThreshold，U 状态才有意义（`nlchain_design.md` §5.2）。它挡的是偏离对角线的随机插值 HSP：
  - J2 = 1000 时，偏离对角线几百 bp 的同向随机 HSP W 走两次跳转的代价约 2.3k，比走两条正常边（2·g(d)）便宜；
  - 但 W 要先满足 s_W ≥ T_solo2 = 2500 才能跳出。
  - 这是 b2 主要的 precision 风险，§8 要看 P_np 的变化。
- **覆盖选项** `--nlchaininner=<J2>[,<T_solo2>]`：
  - 只能与 `--nlchain`、`--inner` 一起用，不能为负；
  - 只给 J2 时，T_solo2 仍按推导取值；
  - 显式给的值不再和外层取 min；
  - 用于 §8 的扫描；扫描定下默认后，这个选项仍然保留。

**对照方案 B（只用于实验）**：保留旧 tweener，只把 `nl_rechain_aligns` 挪到 `tweener_interpolate` 之后。这样插值比对也要过 T_chain，改动约 10 行。它与新方案的差别正好是"gapped 之前插值 + 反向搜索 + 经过 nlchain DP"这部分，用来判断这部分有没有带来收益。

## 5. 代码改动

共改 4 个源文件：`src/tweener.h`、`src/tweener.c`、`src/lastz.h`、`src/lastz.c`。不新增源文件；`tweener.o` 已经在 Makefile 里。`src/nlchain.h`/`nlchain.c` 不改（D6 全局生效，§4.1）。

### 5.1 `src/tweener.h`

```c
#include "nlchain.h"			// non-linear chaining stuff

void tweener_collect_nl_hsps (seq* target, seq* query, const u8* qToComplement,
                              const s8 upperCharToBits[], seed* innerSeed,
                              scoreset* maskedScoring, score xDrop,
                              score innerThreshold, u32 windowSize,
                              const nlchainresult* r,
                              const segtable* stMinus, const segtable* stPlus,
                              segtable** innerMinus, segtable** innerPlus);
```

- 只读两张主表，结果放进调用者给的两张临时表，不直接改主表（§2.4）。
- `tweenerStats` 增加 `numNlWindows`、`numNlWindowsSkipped`、`nlWindowArea`、`numNlHsps[2]`（下标 0 为反链，1 为正链），并在 `show/generic_stats` 中输出。

### 5.2 `src/tweener.c`

1. **文件头注释**：加一段说明 `tweener_collect_nl_hsps` 的用途，并指向本文档。
2. **新增函数** `tweener_collect_nl_hsps`：
   - **setup**：仿照 `tweener_interpolate` 的 setup 段，给这个函数用到的私有全局变量赋值：`seq1`、`seq2`、`upperCharToBits`、`innerSeed`、`maskedScoring`、`xDrop`、`scoreThresh`、`windowSize`。
     - 新建 `tweenSeq1`（target 窗口）、`tweenSeq2`（query 反链窗口）、`tweenRev2`（query 正链窗口；这里借用变量名，存的是反向互补而不是反向）。
     - 设置 `simpleInfo.hp`，与 `tweener_interpolate` 的第 314–328 行相同；`innerAnchors` 用 `new_segment_table (numDefaultAnchors, 0)` 新建。
     - 算出 `maxSub`（§2.8）。
   - **遍历**：对 `r->seg[s]`，只处理 `isSecondary == 0` 且 `count ≥ 2` 的段。
     - 段所在的表 `st = (strand < 0) ? stMinus : stPlus`。
     - 对段内每一对相邻 HSP（k, k+1）按 §6 计算窗口；任一边 ≤ 0 或 > windowSize 时跳过；满足 §2.8 预筛条件时，计入 `numNlWindowsSkipped` 后跳过。
   - **每个窗口**：
     ```c
     extract_subsequence (target, t0, t1, tweenSeq1);
     seq1Positions = build_seed_position_table (tweenSeq1, 0, tweenSeq1->len,
                                                upperCharToBits, innerSeed, 1);

     extract_subsequence (query, r0, r1, tweenSeq2);		// 反链方向（query 当前状态）
     search_window (seq1Positions, tweenSeq2);				// HSP 进 innerAnchors
     append_shifted (innerMinus, t0, r0, query->revCompFlags);

     extract_subsequence (query, r0, r1, tweenRev2);
     rev_comp_sequence   (tweenRev2, qToComplement);			// 正链方向
     search_window (seq1Positions, tweenRev2);
     append_shifted (innerPlus, t0, f0, rcf_forward);

     free_position_table (seq1Positions);
     ```
     - `search_window`：先 `empty_segment_table (innerAnchors)`，再设 `simpleInfo.hp.seq2` 和 `tweenSeq2` 指向本次的 query 片段（`collect_inner_hsps` 的调试输出用），然后调用 `seed_hit_search`，参数与 `bounded_align` 第 799–809 行相同。
     - `append_shifted`：把 `innerAnchors` 中每条 HSP 的 pos1 加 t0、pos2 加 r0（或 f0），`add_segment` 进临时表，id 取所在链方向的 rcf 值（反链为 `query->revCompFlags`，正链为 `rcf_forward`，与 `collect_hsps` 一致），hspId 传 0。
     - 统计：`numNlWindows`、`nlWindowArea`、`numNlHsps[]`。
   - **cleanup**：释放三个序列和 `innerAnchors`。
3. **不改动**：`tweener_interpolate`、`bounded_align`、`collect_inner_hsps` 一律不改（非 nlchain 的路径要逐字节不变）。

### 5.3 `src/lastz.h` / `src/lastz.c`

1. **`nlchain_both_strands`**（第 3735 行）：按 §3 的第 3–5 步修改。
   - 并入 HSP 的保护（D9）写成一个私有函数 `append_inner_hsps (segtable** dst, const segtable* src)`：
     - 保存 `coverageLimit`、`coverage`、`lowScore`；
     - `coverageLimit = 0`，逐条 `add_segment`；
     - 再把三个字段赋值恢复。
   - 两张临时表在函数内新建、用完释放。
   - 插值计时：搜索部分包在 `debugClockInterpolation` 里（nlchain 下旧 tweener 不再运行，这个时钟不会重复计时）；两次串链都计入 `debugClockNlChaining`。
   - `dbgShowAnchors` 时，多打印一行：第一次后的 HSP 数、加入的插值 HSP 数（反链/正链）、第二次后的 HSP 数、J2 和 T_solo2。
   - J2、T_solo2 的取法（§4.1）：
     ```c
     // 第二次串链单独取 J、T_solo，且不比第一次严（docs/nlchain_prechain_inner_plan.md §4.1）
     nl2 = nl;
     if (currParams->nlInnerJump >= 0) nl2.jumpPen = currParams->nlInnerJump;
     else nl2.jumpPen = min (nl.jumpPen, currParams->innerThreshold / 2);
     if (currParams->nlInnerSolo >= 0) nl2.soloScore = currParams->nlInnerSolo;
     else nl2.soloScore = min (nl.soloScore, 5 * currParams->innerThreshold / 4);
     ```
   - `collect_stats` 下的验证计数 `numPass1Dropped`：
     - 第一次压缩后，把两张表复制一份，按 (pos1, pos2, length) 排序；
     - 第二次之后，逐条检查这些 HSP 是否还在表里；被 merge 合并改变了坐标的，按丢掉计（偏保守）。
     - 计数加进 `nlchainStats`，并在 `nlchain_show_stats` 中输出。它不在 D10 的快照恢复范围内。
     - 只在统计构建下复制，默认构建没有开销。
2. **`src/lastz.h`**：参数结构在 `nlr` 之后加两项，`lastz.c` 的 `defaultParams` 在对应位置补初值 `-1, -1`：
   ```c
   score		nlInnerJump;	// (I) --nlchaininner 的 J2；-1 => min(J, innerThreshold/2)
   score		nlInnerSolo;	// (I) --nlchaininner 的 T_solo2；-1 => min(T_solo, 1.25*innerThreshold)
   ```
3. **`finish_one_strand` 的插值调用**（第 3546 行）：条件改为 `(alignList != NULL) && (currParams->innerThreshold > 0) && (!currParams->nlchain)`，并加注释指向本文档。
4. **命令行**：仿照 `--nlchainsolo`（第 6168 行）解析 `--nlchaininner=<J2>[,<T_solo2>]`。
   - 值不能为负；设 `haveNlchainOption = true`，再加一个局部标记 `haveNlchainInner`。
   - 解析两个逗号分隔值的写法照抄 `--nlchainjump`。
5. **参数检查**（第 9366 行起的 nlchain 块）：
   - D8：`if ((lzParams->innerThreshold > 0) && (lzParams->selfCompare)) chastise ("--nlchain can't be used with both --inner and --self\n");`
   - `if ((haveNlchainInner) && (lzParams->innerThreshold <= 0)) chastise ("--nlchaininner requires --inner\n");`
   - D7 不报错，在 `nlchain_both_strands` 中按 `currParams->gappedExtend` 判断是否运行。
6. **`--help`**：本方案不改（见 §2.9），`--nlchaininner` 与其他 nlchain 选项一起补写。

### 5.4 文档

- `docs/nlchain_design.md`：
  - §7 回溯流程的第 5、6 步之间插入"5′. 插值（`--inner`）"，说明窗口、正反搜索、第二次串链；
  - §8.2 加一条，记录调用点与 D7/D8/D9；§8.3 默认值表加 J2、T_solo2 两行（标注"暂定，待 §8 扫描"）；
  - §11 修订记录加一行。
- `docs/nlchain_progress.md`：新增一节，写实现内容、验证结果、实验数字、未解决问题。
- `docs/nlchain_inner_plan.md`：在顶部标注"已被 `nlchain_prechain_inner_plan.md` 取代，未实施"（或删除，由用户定）。

## 6. 坐标换算

记号：
- query 分区的范围为 `[qF, qT)`，即 `qF = sepBefore + 1`，`qT = sepAfter`；无分区时为 `[0, Q)`。
- 正链坐标 f 与反链坐标 r 的区间换算：`[r0, r1) = [qF + qT − f1, qF + qT − f0)`，与 `convert_hsp` 的口径一致（`nlchain.c:1434–1435`）。
- 表中 HSP 的 pos2 都按本链方向的坐标存储：正链表存正链坐标，反链表存反链坐标。

**窗口**：段内相邻两个 HSP，记 u 为在本链方向坐标上靠前的一个，v 为另一个。
- 正链段：u = k，v = k+1。
- 反链段：沿正链 query 前进时，反链坐标是递减的，所以 u = k+1，v = k。

```
t0 = u.pos1 + u.length,   t1 = v.pos1                 （target）
本链方向 query 区间 [a0, a1) = [u.pos2 + u.length, v.pos2)
若本链为正链：f0 = a0, f1 = a1, 由此换算出 r0, r1
若本链为反链：r0 = a0, r1 = a1, 由此换算出 f0, f1
要求 t1 > t0、a1 > a0，且 t1 − t0、a1 − a0 都 ≤ windowSize
```

**HSP 回到全局坐标**：窗口片段内的 HSP 记为 (p1′, p2′, len)。
- 反链方向搜索：`pos1 = t0 + p1′`，`pos2 = r0 + p2′`，进 innerMinus。
- 正链方向搜索：`pos1 = t0 + p1′`，`pos2 = f0 + p2′`，进 innerPlus。

正链片段是反链片段的反向互补：`tweenRev2[i] = comp(query->v[r1 − 1 − i])`，对应正链坐标 `f0 + i`，所以 `pos2 = f0 + p2′` 成立。

**不需要处理分隔符**：窗口在同一个 query 分区、同一个 target 分区内（§2.9），片段中没有分隔符。

## 7. 实施步骤与验证

1. **实现**（§5.1–5.3）。
   - 验证：本地 OrbStack VM 上用 gcc `-Werror` 构建 `lastz`、`lastz_D`、`nlchain_test`、`nlrechain_test`，无警告；`-Dcollect_stats -DdbgTiming` 构建通过（不加 `-Werror`，原因见 `nlchain_progress.md` §10.1）。
2. **回归**，基准为 d9db595 构建的二进制。
   - `make base_tests`、`make nlchain_test`、`make nlrechain_test` 全过。
   - 不加 `--nlchain` 的 22 组选项（含 `--inner`）：stdout、stderr、退出码都逐字节一致。
   - `--nlchain` 不加 `--inner` 的 8 组：把原来的 `--inner` 组换成 query `[multiple]`（§2.6），逐字节一致。
   - 互斥：`--nlchain --inner --self` 报错；`--nlchaininner` 不带 `--inner`、取负值都报错；`--nlchain --inner --nogapped` 正常运行，输出与 `--nlchain --nogapped` 一致（D7）。
3. **构造测试**：生成脚本放进 `test_data/nlchain_inner/`（沿用上一版方案的 D5，避免丢失）。用显式参数把各个门槛拆开，使每个预期都可以判定：
   ```
   --nlchain --hspthresh=50000 --gappedthresh=3000 --nlchainjump=3000 --nlchainsolo=3750 --inner=2000
   ```
   - 外层门槛调高，使下面的 X 主搜索找不到。
   - 第一次串链显式取 h = 3000 时的默认 J、T_solo；λ、D0 用默认值。
   - 按推导：J2 = min(3000, 1000) = 1000，T_solo2 = min(3750, 2500) = 2500。
   - 序列：target = 随机 5 kb + A（3 kb）+ 空隙 G（8 kb）+ B（3 kb）+ 随机 5 kb；query = A′ + G′ + B′，A′/B′ 为 90% 一致的拷贝，G 与 G′ 等长（平衡）。G′ 中放两段：
     - X：G 正中一段的反向互补，长度按分数校准：s_X 约 6500，并且 gapped 分 ≥ 5000（T_chain）。
       - 按 §4.1 的表，8 kb 空隙下两次跳转的代价：J2 = 1000 时约 4.5k，J = 3000 时约 8.5k；
       - 6500 落在两者之间，所以 X 进不进主路径取决于 J2；
       - 脚本按实际坐标用 §2.1 的公式算出两档代价，确认 s_X 落在中间。
     - W：在 G′ 中距离 X 和两端都 ≥ 2 kb 处，放一段约 25 bp 的精确同向拷贝，与 A/B 在同一条对角线上，HSP 分数在 2000–2500 之间。
   - 预期（新版，推导值）：
     - `--nlchainreport` 依次为 `S + A`、`I − X`、`I + B`，X 的坐标与构造位置一致；
     - 输出 MAF 中有 X 的反链比对；
     - W 走正常边留在主路径上（`dbgShowAnchors` 可以看到），但 gapped 分 < 3000，最终输出里没有；
     - 加 `--gappedthresh=2000` 后，W 出现在输出中。
   - 预期（新版，`--nlchaininner=3000,3750`，即 D6a）：
     - 段表只有一段 `S + A…B`，没有 X；
     - W 的情况同上。这说明 X 的进出只由 J2 决定。
   - 预期（新版，只给 `--nlchaininner=1000`）：与推导值的结果逐字节一致（T_solo2 仍按推导取 2500）。
   - 预期（d9db595，去掉 `--nlchaininner`，其余参数相同）：没有 X；W 由旧 tweener 输出。
   - 变体（每个都检查坐标与构造一致）：
     - query 整体反向互补后再跑一次：段表的 σ 对调，比对坐标与正向运行一致；
     - target 用 `[multiple]`（两条染色体，A/G/B 都在第二条上）；
     - query 用 `[multiple]`（两条 query，构造只在第二条上），段表与单独跑这条 query 时一致；
     - 自适应 `--hspthresh=top50%`：在构造中让 W 能被找到；验证 W 的比对不会因为门槛被拉低而输出（D9）。用一个临时调试输出确认 `lowAnchorScore` 加 HSP 前后相同，确认后去掉。
   - `-Dcollect_stats` 构建：`numNlWindows`、`numNlWindowsSkipped`、`numNlHsps[]` 与构造的预期一致，`numPass1Dropped` = 0。
4. **pseudocat/aglobin 冒烟**：`--nlchain --inner=2000` 正常运行，统计自洽：第二次的输入 HSP 数 = 第一次的输出 + 插值 HSP 数 − merge 合并掉的。
   - 默认参数和次级链 `--nlchainsecondary=8800,1,0.9` 各跑一次，记录 `numPass1Dropped`。
   - 默认参数下预期为 0。如果不为 0，先停下分析原因，再决定是否需要"只对插值 HSP 生效"的做法（§4.1）。
5. **远端实验**（§8）。
6. **记录**：更新 `docs/nlchain_progress.md`、`docs/nlchain_design.md`；提交前先给用户看文件列表、diff 摘要和 commit message。

## 8. 远端实验设计

数据与流程同 `nlchain_progress.md` §10.3：simTest 20 个作业，mafComparator `--samples=10000000 --seed=1`，两份真值，`tune_report.py` 汇总。二进制按 `bootstrap_nlr.sh` 的做法构建，sha256 记入 MANIFEST。

| 组 | 二进制 | 参数 | 作用 |
|---|---|---|---|
| G0 | 新 | `--nlchain` | 应与 `nlr_default` 逐字节一致（只差命令行回显） |
| G1 | d9db595 | `--nlchain --inner=2000` | 旧插值 |
| G2 | 新 | `--nlchain --inner=2000` | 本方案，D6 取推导值 b2（J2 = 1000，T_solo2 = 2500） |
| G2b1 | 新 | `--nlchain --inner=2000 --nlchaininner=2000,2500` | D6 取 b1 |
| G2a | 新 | `--nlchain --inner=2000 --nlchaininner=2200,2750` | D6 取 a（simTest 的 h = 2200 下与外层相同） |
| G3 | 新 | `--nlchain --inner=2000 --gappedthresh=2000` | D5 的门槛影响 |
| G3b | 新 | `--nlchain --gappedthresh=2000` | 确认降门槛几乎不影响外层比对（§2.2） |
| G4（可选） | 方案 B 补丁 | `--nlchain --inner=2000` | 对照：旧插值 + 过 T_chain（D12） |
| G2p | 新，`-Dcollect_stats -DdbgTiming` | 同 G2 | 计数与各阶段耗时 |

**比较**：
- R/P/F（np、anc）、tMul/qMul、CPU、interpolation 与 nlchaining 的耗时；
- 插值窗口数、预筛跳过的比例、插值 HSP 数（正/反）；
- `--nlchainreport` 中 I 段数：G2、G2b1、G2a 各自对 G0（量化 §2.1 的好处 1，以及 J2 的作用）；
- G2p 的 `numPass1Dropped`（§4.1 的"全局 ≈ 只作用于插值 HSP"是否成立）。

**判定**：
- G2（或 G3）的 F_np、F_anc 都高于 G0 和 G1，才算有效。
- D6 的默认值按 G2 / G2b1 / G2a 中 F_np 与 F_anc 的平均值来定；三者相差不超过 0.001 时取 a，参数最少、行为与外层一致。
- 如果 G4 与 G2 相当，说明收益来自"过 T_chain"，而不是 gapped 之前插值，应改用更简单的方案 B。
- 新增耗时以 G2p 计，记录占总 CPU 的比例；先不设门槛，看数据再定。

## 9. 风险

- **收益可能很小**：
  - 倒位一项受 §2.1 的 DP 代价限制；
  - D5a 下，同向弱片段要 gapped 分 ≥ 3000 才能留下；
  - 被外层 gapped 吞掉的插值 HSP 没有贡献。
  - 三者叠加，G2 可能与 G0 接近。G3、G4 就是为了区分原因而设的。
- **precision**：
  - 7-mer 门槛 2000 在 20 kb 窗口里会产生随机 HSP。DP 只在偏离对角线约 100 以内时才会把它们串进来（g(d) 的两次代价低于 s_W），之后还要经过 gapped 门槛和剥链。
  - 旧 `--inner` 在 baseline 上曾使 precision 下降（`lastz_profiling.md` §6）。
- **J2 降低带来的 precision 风险**（§4.1）：
  - 偏离对角线、分数 ≥ T_solo2 的同向随机插值 HSP，可以通过两次便宜的跳转进入主路径；
  - 分数 ≥ 约 4k 的随机反向 HSP，也可以作为倒位进入主路径。
  - 它们之后还要过 gapped 门槛和 T_chain。孤立的反向片段，gapped 分不到 5000 会被删掉。
- **第二次串链改变主路径**：插值 HSP 可能让 DP 改走另一条路径。第二次的主路径分不低于第一次，但段的划分可能变化。用段表对比 G0 和 G2 来观察。
- **merge 后 hspId 不可靠**：`merge_segments` 只写 pos/length/s，id/hspId 沿用目标槽位的旧值（`segment.c:1572–1576`、`1587–1591`）。统计不要依赖 hspId 来标记插值 HSP。
- **`--masking`（动态掩蔽）**：旧版正链的插值在反链输出并掩蔽之后才运行；新版两条链都在掩蔽之前搜索。只在 `--masking` 下有差别，记入文档，不处理。

## 10. 待用户确认

1. D5、D7–D12 的推荐项（尤其是 D8 报错、D12 是否做 G4）；D6 的 §4.1：暂定默认 b2、全局生效、新选项 `--nlchaininner`。
2. `docs/nlchain_inner_plan.md` 标注为"已取代"还是删除。
3. 是否需要另开方案处理 §2.1 的根因（I 型跳转代价）。
