# nlchain 方案 A 实施方案：gapped 之后重新串链，加链分下限

日期：2026-11-02。

**状态**
- 已定方向（用户确认）：方案 A 在 `--nlchain` 下默认开启；次级链按用途选择。
  - 只要 1:1 直系同源：用默认参数。
  - 需要旁系同源：加 `--nlchainsecondary=<4h>,1,0.9`。
- 第 5 节的实现细节，要由新会话按全局规范先列出、经用户确认后再编码。
- **2026-11-02**：D1–D6 全部按推荐确认，另确认补充细节 X1–X6 和对拍门槛 0.995。§4 第 1–6 步已在本地完成，结果、补充细节和与本方案的差异见 `docs/nlchain_progress.md` §9。第 7 步起尚未做。

**依据**：`docs/nlchain_precision_gap.md`，包括根因、模拟数据和方案比较。

## 1. 目标

1. **补上 precision 差距**。nlchain 的 gapped 比对里夹着 50–1000 bp 的错位碎片：比对总分达标，但其中一段与别的比对重复覆盖同一段碱基。chainNet 流程用 axtChain 重新串链、再用链分下限把这些碎片删掉，nlchain 目前缺这一步。
2. **在 lastz 内部做这一步**，不依赖 kent 工具。算法以 axtChain 为准，原因有两个：
   - 模拟用的就是 axtChain，数字有据；
   - 偏离 axtChain 就要重新验证。
3. **效果目标**：与模拟值相符。模拟方式为 axtChain `-minScore=0 -linearGap=medium` 加 HoxD55，然后保留链分 ≥ 5000 的链，不跑 chainAntiRepeat。

   simTest 10 个物种对的平均 F：

   | 配置 | F_np | F_anc |
   |---|---|---|
   | 默认参数 + A | 0.7247 | 0.7167 |
   | `--nlchainsecondary=8800,1,0.9` + A | 0.7208 | 0.7194 |
   | 对照：原版 + chainNet medium/3000 | 0.7205 | 0.7191 |
   | 对照：现状 nlchain | 0.7085 | 0.7054 |

## 2. 算法（逐步对齐 axtChain）

下面写的是 kent 源码的行为，2026-11-02 已从 `src/lib/chainBlock.c`、`chainConnect.c`、`gapCalc.c` 和 `src/hg/mouseStuff/axtChain/axtChain.c` 核对过。

kent 源码有自己的许可条款，**不要拷贝代码**，按本节的描述独立实现。lastz 本身是 MIT 许可。

**输入**：一次 `finish_one_strand` 中 `gapped_extend` 之后、通过各项过滤的 `alignList`。这一次调用只涉及一条链方向；反链调用时 query 已经反向互补，坐标直接用 alignel 里的坐标即可。

**分组**：按 (target 分区, query 分区) 分组，链不跨序列。axtChain 按 (qName, strand, tName) 分组，与此一致。

**步骤**

1. **拆块**。比对 edit script 里每一段连续的 `editopSub` 是一个无 gap 块，记录：
   - `tStart/tEnd`、`qStart/qEnd`（半开区间）；
   - 来源比对；
   - 在来源比对中的次序。
2. **块分**：块内替换分之和，用 `scoring->sub`，不含 gap 罚分。axtChain 用 `axtScoreUngapped`，矩阵相同，都是 HoxD55。
3. **合并完全重复的块**：按 (qStart, tStart) 排序，起点完全相同的块并成一个，终点取较大者（对应 axtChain 的 `removeExactOverlaps`）。
4. **DP**：块按 tStart 升序处理。每个块的初始分 `total = score`。
   - 合法前驱 j：`qStart_j < qStart_i` 且 `tStart_j < tStart_i`，两个都是严格小于。
   - 候选分 = `total_j + score_i − connectCost(j, i)`，比当前 `total_i` 大才更新。
   - `connectCost(a, b)`：
     - 令 `dq = b.qStart − a.qEnd`，`dt = b.tStart − a.tEnd`。
     - 若 dq < 0 或 dt < 0：重叠量 `ov = −min(dq, dt)`。
       - ov ≥ 任一块的长度（完全包含）时，代价取 1e8，相当于禁止。
       - 否则加上交叉点调整（见下），再令 `dq += ov`、`dt += ov`。
     - 再加 `gapCalcCost(dq, dt)`。
   - 交叉点调整（`cBlockFindCrossover`）：
     - 在重叠的 ov 列上，先假设全部用右块（分数为 rScore），再从左往右逐列"改用左块"，记录总分最高的位置 bestPos，平局取靠前的位置。
     - 调整量 = `rScore + lScore − bestScore`，截断为整数。
   - `gapCalcCost(dq, dt)`：负值先归零。
     - dt = 0：用 qGap 行，按 dq 查表；
     - dq = 0：用 tGap 行，按 dt 查表；
     - 两者都 > 0：用 **bothGap 行，按 dq + dt 查表**。
     - 表内线性插值，结果截断为整数；1..110 先预计算；超出表尾按最后一段的斜率外推（medium 为每 bp 1.0）。
     - medium 表：position 1 2 3 11 111 2111 12111 32111 72111 152111 252111；qGap = tGap = 350 425 450 600 900 2900 22900 57900 117900 217900 317900；bothGap = 750 825 850 1000 1300 3300 23300 58300 118300 218300 318300。loose 表见 kent `gapCalc.c`。
   - 前驱搜索：axtChain 用 2-d 树，结点记子树内的 `maxQ`、`maxT`、`maxScore`，用"子树最大分 + score_i − gapCost(下界)"剪枝。实现可以照这个思路做，也可以先用暴力 O(n²) 做参考解。
5. **剥链**：
   - 块按 `total` 降序排列。
   - 从每个还没用过的块出发，沿前驱回溯并标记已用；遇到没有前驱、或前驱已被更高的链用掉时停止。
   - 不重跑 DP。
6. **链内整理**，按以下顺序：
   1. 部分重叠裁剪（`chainRemovePartialOverlaps`）：
      - 链内相邻块 a、b 重叠时，在交叉点 c 切开：`b.start += c`，`a.end −= (ov − c)`。
      - 若一块完全落在另一块里，删掉 b。
      - 裁剪后若 a 为空，删掉 a 并从头再扫。
   2. 合并首尾正好相接的相邻块（`chainMergeAbutting`）。
   3. 重算链分（`chainCalcScore`）：块分之和减去相邻块之间的 `gapCalcCost`。
7. **下限**：链分 < T_chain 的链整条删除（`axtChain -minScore` 也是用重算后的分数做判断）。
8. **输出**（推荐做法 (a)，见 §5 D2）：
   - 把保留的块放回各自的来源比对。同一条链中、在来源比对里又是连续的块（可以包括被裁短的块），合成一段，块与块之间沿用原比对的 gap 操作。
   - 每段构造成一个新的 `alignel`：beg/end 从 1 起、闭区间；edit script 用 `edit_script_new` 和 `edit_script_add` 构造；分数用 `score_alignment` 重算。
   - 删掉原 `alignList`，换成新列表。
   - 新列表的碱基对集合，与保留链的块的碱基对集合完全相同。

**与模拟的已知差别**：axtChain 的 HoxD55 矩阵对 N 等非 ACGT 碱基的计分，可能与 lastz 的 `--ambiguous=iupac` 不同，影响很小。验收时用逐碱基对比对来核对（§4 第 7 步）。

## 3. 代码落点

- **新文件** `src/nlrechain.h`、`src/nlrechain.c`：第 2 节的步骤 1–8，对外只暴露一个入口，例如：

  ```c
  alignel* nl_rechain_aligns (seq* target, seq* query, scoreset* scoring,
                              alignel* alignList, const nlrechainparams* params);
  ```

  理由：这是一个独立的后处理，与 `nlchain.c` 的 HSP 级 DP 是两套模型（代价和重叠处理都不同，见 §5 D1），放在一起会让 §8.4 要求的"转移方程一眼能对上公式"变得难读。
- **`lastz.c`**
  - 调用位置：`finish_one_strand` 中 identity/coverage 等过滤之后、`tweener_interpolate` 和 `print_align_list` 之前，条件为 `currParams->nlchain && T_chain > 0`。
  - 这样之后的动态掩蔽（census）也用删过碎片的比对，与 chainNet 输出的语义一致。是否接受这一点，见 §5 D6。
  - 还要加：参数、默认值、CLI 解析、互斥检查（只在 `--nlchain` 下允许）、`dbg_timing` 新增 `debugClockNlRechain`，以及 `-Dcollect_stats` 下的计数（块数、链数、删除的链数和碱基对数）。
- **`src/lastz.h`**：`lastzParams` 增加 rechain 参数。
- **`src/Makefile`**：加入新文件；单测可以并入 `make nlchain_test`，也可以新开一个目标。
- **CLI**（名字待确认，见 §5 D3）：`--nlchainclean=<T_chain>`，默认 5000，取 0 关闭。gap 表沿用 `--nlchaingapcost=medium|loose`。
- **调试输出**（可选，建议做）：`--nlchaincleanreport=<file>` 以 kent chain 格式写出剥链结果（删除前，并标出哪些链被删）。用于与 axtChain 逐链对拍。

## 4. 实施步骤（每步附验证方式）

1. **gap 表与代价函数**：medium/loose 两张表的三行，`gapCalcCost`，交叉点函数。
   - 验证：单测核对 11 个表结点、插值截断、外推，以及按 dq + dt 查 bothGap；交叉点用手算样例，覆盖平局取靠前位置。
2. **拆块与块分**。
   - 验证：单测检查所有块的碱基对并集等于原比对的碱基对；块分之和减去原比对的 gap 罚分等于 `score_alignment`。
3. **DP、剥链、链内整理、下限**。
   - 验证：随机小块集（n ≤ 12），每个块的最优前驱和 total 与暴力枚举一致；剥链结果与按 §2 第 5 步手写的参考实现一致；部分重叠裁剪覆盖"完全包含"和"裁空"两种情况。
4. **输出重建**。
   - 验证：单测检查输出比对的碱基对等于保留链块的碱基对，每条输出比对的分数等于 `score_alignment`。
5. **接入 lastz.c**。
   - 验证：本地 `make`，并用 gcc 加 `-Werror` 构建（注意 7939842 那次 `-Wempty-body` 的教训）；`make base_tests` 全过；不加 `--nlchain` 时 22 组选项的输出与原版逐字节一致；`--nlchain --nlchainclean=0` 与当前 nlchain（06533fd）逐字节一致。
6. **本地冒烟**：pseudocat vs pseudopig、aglobin。
   - 验证：输出合法，`-Dcollect_stats` 下的计数自洽（删除的碱基对 + 保留的碱基对 = 原有碱基对，扣除交叉点裁掉的部分）。
7. **远端与 axtChain 对拍**（simTest）。
   - 输入：`runs/tune/nl_step5default/` 的 MAF 已经有对应的 axtChain 结果 `runs/prec_sim/nl_step5default.medium0/<物种对>/chain_raw/`。
   - 新二进制加 `--nlchainreport`、`--nlchaincleanreport` 跑 20 个作业。
   - 验证：保留链的碱基对集合与"axtChain + ≥5000"（`runs/prec_sim/xover/<物种对>/raw5000.chain`）的 Jaccard ≥ 0.995（门槛待确认）。不达标时先查差异链，不急着调参。
8. **远端评测**：两套配置各跑 20 个作业，都在两份真值上评；还要记录 CPU。
   - 两套配置：默认参数；`--nlchainsecondary=8800,1,0.9`。
   - 验收：
     - 各项 F 与 §1 的模拟值相差不超过 0.002；
     - 默认配置 F_np > 0.7205；
     - 新增阶段的 CPU 不超过 lastz 总 CPU 的 1%。
9. **（可选）在 A 之上重扫次级链参数**：T_sec ∈ {3h, 4h, 6h}，N_sec = 1，r_sec = 0.9。A 删掉碎片后，开次级链的代价可能变小，推荐值可能会变。
10. **文档与提交**：更新 `docs/nlchain_progress.md`、`docs/nlchain_design.md`（§7 回溯流程加一步、§8 接口），补写 `README.lastz.html` 和 `lastz --help`。commit 前先展示内容，经用户确认。

16 个最终测试物种对，等第 8 步通过后由用户决定何时跑。

## 5. 编码前要请用户确认的设计选择

| 编号 | 问题 | 选项 | 推荐 |
|---|---|---|---|
| D1 | 链化模型 | (a) 照搬 axtChain：bothGap、交叉点调整、不限 gap 长度、独立 2-d 树；(b) 复用 `nlchain.c` 的正常边和 K-d 树：g(d) 只有 qGap 行，没有 bothGap；重叠用 ρ 撤回而不是交叉点；有 maxGap/maxDiag 限制 | **(a)**。模拟数字来自 axtChain；(b) 的代价模型不同，要重新验证，省下的代码量也有限 |
| D2 | 输出形式 | (a) 把原比对在被删块处切开，分段输出；(b) 按链重建比对，块之间用 gap 连接 | **(a)**。改动小，碱基对与 (b) 相同；(b) 的输出更像 chainNet，但要处理跨比对的两侧 gap |
| D3 | 参数名与默认值 | `--nlchainclean=<T_chain>`，默认 5000，0 关闭 | 名字可改。5000 已验证：3000–12000 之间 F 只变化 0.0007 |
| D4 | T_chain 用绝对值还是随 hspthresh 缩放 | 绝对值 / k×hspthresh | **绝对值**（同 chainAntiRepeat）。目前只在 h=2200 上验证过 |
| D5 | 是否也跑 chainAntiRepeat 的重复序列和低复杂度调整 | 不跑 / 跑 | **不跑**。模拟中跑与不跑只差 0.0005（0.7242 对 0.7247） |
| D6 | 动态掩蔽（`--masking`）用删碎片前还是删碎片后的比对 | 删碎片后 / 删碎片前 | **删碎片后**。但这会改变后续 query 的掩蔽，模拟时没有体现；第 8 步的结果会包含这一效应 |

## 6. 风险

- **细节差异**：自己实现的链化与 axtChain 在细节上可能不一致，比如整数截断、平局规则、N 碱基计分。第 7 步的对拍就是为此设的。
- **动态掩蔽**：D6 会让结果与模拟不完全一致，因为模拟是在固定的 nlchain 输出上做后处理的。
- **输出变化**：比对被切成更多、更短的段，块数会增加。下游如果按块数或块长统计，需要知道这一点。
- **验证范围**：近缘物种对和全基因组规模都没有验证过。性能预期仍是线性，但要在第 8 步之后用全量数据确认。
