# lastz 非线性串联（nlchain）可行性与设计

日期：2026-10-31。状态：设计稿，未实现；同日按源码评审修订一次（要点见 §11）。已确定：正常边用 K-d 树（§6.1 方案 A），旁系同源用 L 状态加次级回溯（§5.4），正常边 gap 代价用 axtChain 式查表（§5.1）。参照对象：DuetMap 两状态 chainer（`/Users/ncl3/Projects/map/src/mapping/chainer.cpp`，设计文档 `docs/chainer-two-state.md`）。

## 0. 结论先行

**可行**，lastz 里已经有的东西能支撑大部分改动；但有两点要先说清楚，它们决定了这件事值不值得做：

1. **提速上限约 1.5×。** nlchain 和现有 `--chain` 一样，作用是在 gapped 延伸**之前**筛掉 HSP。`lastz_profiling.md` §5 显示 gapped 延伸只占 32%，种子查询加无 gap 延伸占 67%，这一步不受影响。所以最多快 1/(1−0.32) ≈ 1.47×（实测 `--chain` 为 1.46×）。nlchain 的价值在**准确度**：现有 `--chain` 的 recall 只有 0.236（baseline 0.647；两者都是在含旁系同源的全真值上评的，`--chain` 的 1:1 输出还没在无旁系同源版真值上评过），nlchain 的目标是在速度相近的前提下把 recall 拿回来，同时内置一份 query 端 1:1 的输出，替代 axtChain/chainNet。
2. **主路径是 query 端 1:1。** 一条 query 单调路径天然不含 query 端的旁系同源，只能和"无旁系同源版真值"公平比较；但它可以经 D 回跳重复使用同一个 target 位点，所以不是 target 端 1:1（评测侧别见 §9 第 0、4 步）。旁系同源由 §5.4 的次级链覆盖：在同一次 DP 里多维护一个只走共线边的局部状态 L，再按 minimap2 的方式从未访问的锚点多次回溯。**不能直接在非线性状态上多次回溯**，原因见 §5.4。

建议先做 §9 第 0 步（不写代码，估算 1:1 选择的 recall 上限）。如果上限和 `--chain` 在同一真值上的 recall 差不多，就没有必要继续做。

## 1. 现状：lastz 的 `--chain` 做了什么

- 入口在 `finish_one_strand()`（定义在 `src/lastz.c:3262`）：先 `merge_segments`（3299 行）、打分，再调用 `try_reduce_to_chain()`（3349 行），然后才做 `gapped_extend()`。
- `reduce_to_chain()`（`src/chain.c:497`）的做法是：按 pos1 排序，用 K-d 树（轴为 pos2 和对角线）找最佳前驱，回溯出**一条** pos1、pos2 都严格递增的链，链外的 HSP 全部 `filter`。
- **已知缺陷**：`best_predecessor()` 在 pos2 切分节点上的两处递归调用（`chain.c:963-964`）把 `axis` 和 `lowerBound` 传反了（对照 976 行和 168 行的声明）。G=R=0 时 lowerBound 恒为 0，碰巧无害；G 或 R 非 0 时会把对角线切分节点当成 pos2 切分来剪枝，漏掉前驱。所以现有 0.236 不受影响，但用 `--chain=G,R` 做对照的结果不可信，复用 `chain.c` 时这段也不能照抄。
- 粒度：**每条链方向 × 每个 (target 分区, query 分区) 批次一条链**（`try_reduce_to_chain` 按分区切批）。正反链在不同的 `finish_one_strand` 调用里各算各的。
- 连接罚分 `chain_connect_penalty()`（`src/lastz.c:3687`）：`diagDiff*G + numSubs*R`，重叠时按匹配分撤回。**默认 G=R=0**，也就是任意单调子集都不扣分。

recall 掉到 0.236 的原因因此很明确：人–鼠一对染色体之间的同线性被打碎成大量块，块之间有顺序和方向的变化，单调链只能保留其中顺序一致的一部分；反链块和正链块也不能连在一起。

## 2. 可以直接复用的 lastz 机制

| 需求 | 现有机制 | 位置 |
|---|---|---|
| 两条链的 HSP 放进同一张表 | `collectHspsFromBoth`：adaptive HSP 阈值时已经这样做，seg->id 记录链方向 | `lastz.c:1426`、`split_anchors()` `lastz.c:3618` |
| 合并重复/重叠的 HSP | `merge_segments`：**只能用于单链表**。它只按对角线号合并、不看 `id`，写回时也不设置 `dstSeg->id`，用在双链合表上会把正反链 HSP 合并并打乱链方向标记 | `segment.c:1527` |
| 判断 HSP 在哪条链 | `seg->id & rcf_rev` | `segment.c:2056` |
| 两条链收集完后拆成两张单链表，DP 后各自做 gapped | `split_anchors()` + `swap_anchor_sets()` | `lastz.c:1670–1690` |
| target 分区（`[multiple]`，跨染色体） | `lookup_partition(seq, pos)`、`partition.sepBefore/sepAfter` | `sequences.h:240,697` |
| 共线前驱的精确搜索 | K-d 树 + 子树最大分剪枝 | `chain.c` `build_kd_tree/best_predecessor` |
| 重叠只撤回一次（沿一个轴） | `chain_connect_penalty` 的 numSubs < 0 分支 | `lastz.c:3687` |

**需要实测核对的一点**：反链 HSP 的 pos2 应该是反向互补 query 上的坐标（种子搜索在 `rev_comp_sequence` 之后进行）。§9 第 2 步用一个人造倒位来验证。

## 3. 与 DuetMap 的差异，以及因此要改的地方

| 维度 | DuetMap（read mapping） | lastz WGA | 设计取舍 |
|---|---|---|---|
| 锚点 | k-mer，span 15，分数 = 覆盖碱基数 | HSP，长几十到几百 bp，分数为替换矩阵分（HOXD70 每个匹配约 100），≥ hspthresh | 用 HSP 分数 s 作为收益，所有罚分都用**同一套矩阵分单位** |
| query 长度 | 约 20 kb，单一来源 | 一条染色体，数十到数百 Mb | 路径起点和终点都自由（局部路径）；query 上的长空白靠跳转跨过 |
| 共线判定 | dq ≤ 2000，\|dq−dt\| ≤ 50 | HSP 之间间隔可达数 kb，indel 也更大 | 上限改为参数，默认放宽到 kb 级 |
| 每段最少证据 | 段内至少 2 个锚点（U→R） | 一个高分 HSP 本身就是可信证据 | U/R 改为"段内 ≥ 2 个 HSP **或** 单个 HSP 分数 ≥ T_solo" |
| 局部候选 | 下标窗口 W=1000 | HSP 密度随区域变化很大，下标窗口对应的 query 距离不稳定 | 正常边改用 K-d 树精确搜索（§6.1，已选定方案 A） |
| 备选/旁系 | `score_segments` 按段起点从 DP 状态收集备选 | 旁系同源是 recall 的重要来源 | 非线性状态的前驱会被主路径"截胡"（§5.4），改为独立的局部共线状态 L 加次级回溯 |
| 跳转候选 | Top-H 池（H=500），有损 | 同上 | 利用跳转罚分有上界，做**精确**的界剪枝池（§6.2） |
| 处理单位 | 一条 read 两条链一起 | 两条链分开处理 | 强制 `collectHspsFromBoth`，在两条链都收集完之后做 DP |
| 跨染色体 | target_id | 必须用 `target[multiple]` 把全部 target 染色体放进同一次运行 | 否则 1:1 路径会把 query 碱基强行分到当前这一条 target 染色体上，造成假阳性。现有实验是 36 个染色体对作业，**nlchain 的运行方式必须改成每条 query 染色体 × 全基因组 target** |

## 4. 统一坐标与记号

对每个 HSP i（按链方向分别 `merge_segments` 之后，见 §8.2 第 5 步；HSP 在收集时已带分数，seed_search.c 设置 `haveScores`）定义：

```text
len_i, s_i                       长度、HSP 分数；ρ_i = s_i / len_i 为平均每碱基分数
σ_i ∈ {+1,-1}                    +：正链；-：seg->id & rcf_rev
p_i                              target 分区下标（无分区时恒为 0）
qb_i, qe_i                       正链 query 上的半开区间
    σ=+ : qb = pos2,            qe = pos2 + len
    σ=- : qb = Q - pos2 - len,  qe = Q - pos2       （Q 为 query 长度，pos2 为反向互补坐标）
tb_i, te_i                       target 半开区间 [pos1, pos1+len)（全局坐标，跨分区）
tin_i, tout_i                    沿 query 正向行走时进入/离开 target 的位置
    σ=+ : tin = tb, tout = te
    σ=- : tin = te, tout = tb
```

DP 按 `(qb, p, tb, σ, len)` 升序处理（与 DuetMap 一致，`qb` 相同的一批一起处理）。对 j→i 定义：

```text
dq  = qb_i - qe_j                         query 上的间隔（负数表示重叠）
dt  = σ·(tin_i - tout_j)                  沿链方向 target 上的间隔（仅同 p、同 σ 时有意义）
ov  = max(0, -dq)                         query 重叠
ov_t = max(0, -dt)                        target 重叠（仅正常边用）
gain(j,i) = s_i - ρ_i · ov                重叠部分不重复计分（对应 DuetMap 的 b(j,i)）
ov_max_i = f_ov · len_i                   i 的任一前驱的重叠上界（搜索层用）
```

所有边都必须满足 `qb_j < qb_i`、`qe_j < qe_i`（query 上有推进），且 `ov ≤ f_ov · min(len_i, len_j)`（因此 `ov ≤ ov_max_i`）。跳转边另外要求 `ov = 0`（`qe_j ≤ qb_i`）：跳转两侧本来就不共线，没有理由重叠，这也让 §6.2 的待释放条件变成精确的。

`merge_segments` 合并后的分数取各片段的最大值，不重新打分，所以合并后的 ρ_i 偏低，重叠撤回会略少（重叠部分有少量重复计分）。现有 `--chain` 也是这样，影响很小，不单独处理。

## 5. DP 模型

### 5.1 两类边

**正常边 N（共线）**：同时满足

```text
p_j = p_i,  σ_j = σ_i,  σ·(tin_i - tin_j) > 0          target 起点沿链方向严格推进
-ov_max_i ≤ dq ≤ maxGap,  -ov_max_i ≤ dt ≤ maxGap
d = |dq - dt| ≤ maxDiag                                 对角线偏移（≈ indel 长度）
```

代价：

```text
C_N(j,i) = g(d) + R·max(0, min(dq, dt)) + ρ_i·max(0, ov_t - ov)
g(d)     = 0                                  d = 0
         = 按 gap 代价表在 d 上分段线性插值      d > 0（表外按最后一段斜率外推）
```

- `g(d)`：对角线偏移 d 约等于一个长度为 d 的 indel。形式和默认表照搬 axtChain 的 linearGap（`kent/src/lib/gapCalc.c`），提供两张内置表：`medium`（axtChain 用于人–鼠）和 `loose`（用于人–鸡），与 `wga_eval_plan.md` 近缘/远缘两套参数对应。`medium` 的 qGap 行为：

  | d | 1 | 2 | 3 | 11 | 111 | 2111 | 12111 | 32111 | 72111 | 152111 | 252111 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | g(d) | 350 | 425 | 450 | 600 | 900 | 2900 | 22900 | 57900 | 117900 | 217900 | 317900 |

  分数单位与 lastz 的替换矩阵分一致（axtChain 本来就用 lastz 的打分）。
- **不用碱基级的 affine（O_c + 30·d）**：gap extend 30 是碱基级罚分，搬到 HSP 层太陡。按旧初值，d≈88 bp 时正常边代价就超过跳转代价，超过约 90 bp 的 indel 全被切成跳转，maxDiag 形同虚设。
- `R·min(dq,dt)`：两个 HSP 之间未被 HSP 覆盖、但仍在同一对角带上的区域（与 lastz 的 numSubs 项含义相同）。axtChain 对两侧都有间隔的情形另收 bothGap，这里先用 R 近似，默认 0，§9 第 3 步再看是否需要。
- 最后一项：重叠只撤回一次，总撤回量为 `ρ_i·max(ov, ov_t)`。`gain` 已经撤回了 query 端的 `ρ_i·ov`，这里只补 target 端超出的部分。旧写法 `ρ_i·max(0,-dt)` 在两端同时重叠时会撤回两次（dq = dt = −k 时扣 2ρk），与 `chain_connect_penalty` 只沿一个轴撤回的做法不一致。
- R 与现有 `--chain=<diag,anti>` 的 anti 含义一致。g(d) 不是线性的，与 `--chain` 的 diag 不同义，所以 nlchain 不复用 `--chain` 的参数。
- 两类边的代价交点：用 `medium` 表和 §8.3 的跳转初值，d ≈ 2.6 kb 时 g(d) 与 P_J 相等，更大的 indel 走跳转（记为 `G`）更便宜。这只对 R 状态的前驱成立；U 状态不能跳出，只能走正常边，所以 maxDiag 实际约束的是 U 前驱和候选搜索范围。取 10 kb（此时 g ≈ 18.7k，已远大于跳转代价上限）。

**异常边 J（跳转）**：任意 j→i（只需满足 §4 的推进约束，且 query 上不重叠：`ov = 0`）。代价：

```text
P_J(j,i) = J + λ · min(L, log2(1 + D(j,i)/D0))

D(j,i) = | dt - dq | = d                    同 p、同 σ（长 indel 或回跳，回跳时 dt 为负）
       = | tin_i - tout_j |                 同 p、换向（倒位）
       = +∞（取 L）                          不同 p（易位、跨染色体）
```

这与 DuetMap 的 `abnormal_penalty` 同形：J 为基础代价，跨度项随几何距离以对数增长，并且有上限 L。同一对 (j,i) 两种边都合法时，DP 自然取较优的一种。

### 5.2 状态

| 状态 | 含义 | 可以做什么 |
|---|---|---|
| `U[i]` | 以 i 结尾的最优非线性路径；i 所在段**只有 i 一个 HSP**，且 s_i < T_solo | 只能沿正常边继续；不能跳出，不能作为终点 |
| `R[i]` | 以 i 结尾的最优非线性路径；i 所在段至少有 2 个 HSP，或 s_i ≥ T_solo | 可以沿正常边继续，可以跳出，可以作为主路径终点 |
| `L[i]` | 以 i 结尾的最优**局部共线链**分数：只走正常边，起点自由，不继承任何跳转前的分数 | 次级链终点（§5.4）；跳出门槛 T_src 的判据（§6.2） |

U/R 的设计让每一段要么至少有两个共线 HSP，要么有一个足够强的 HSP。目的是防止单个弱 HSP（通常是重复序列的随机命中）被当成"跳转孤岛"插进路径。这要求 T_solo 明显高于 HSP 的入选门槛：每个 HSP 都满足 s ≥ hspthresh，T_solo 如果不高于 hspthresh，U 状态就永远用不上（§8.3）。

L 与 U/R 相互独立：L 不读 U/R，U/R 也不读 L 的分数（只在 T_src 判断时读 L，见 E5）。L 就是 minimap2 式的线性局部 chain，和主路径共用同一次候选搜索。

### 5.3 转移方程

方程编号 E1–E6 在代码注释中引用（§8.4）。

```text
(E1) 起点      U0[i] = s_i
(E2) 跳入      UJ[i] = max_{j ∈ Jump(i)}   [ R[j] + gain(j,i) - P_J(j,i) ]
(E3) 共线延伸  RN[i] = max_{j ∈ Normal(i)} [ max(U[j], R[j]) + gain(j,i) - C_N(j,i) ]
(E4) 定状态    若 s_i ≥ T_solo： R[i] = max(RN[i], U0[i], UJ[i])，U[i] = -∞
               否则：           R[i] = RN[i]，              U[i] = max(U0[i], UJ[i])
(E5) 跳出资格  Jump(i) 只含 L[j] ≥ T_src 的 j（T_src = 0 时不限制）
(E6) 局部共线  L[i] = max( s_i,  max_{j ∈ Normal(i)} [ L[j] + gain(j,i) - C_N(j,i) ] )
```

- `Normal(i)`、`Jump(i)` 是 §5.1 定义的合法前驱集合；`gain`、`C_N`、`P_J` 见 §4、§5.1。
- 只有 R 能跳出，所以不会出现连续的异常边。
- E5 用 L 判断，而不是用"保留赢家路径的段内分"：L[j] 是 j 结尾的最优共线链分数，与 U/R 保留了哪条路径无关，所以这个门槛是精确的。DuetMap 用 `segment_support` 做同类门槛，`chainer-two-state.md` 记录过它因此丢掉合格候选的反例，这里不会出现。
- 输出用的段内分不另存：主路径上每段的分数在回溯时按 E3 的边代价重新累加（§7）。
- 平局规则（保证输出确定），写在**一个**比较函数里，所有状态共用：分数高者优先，其次段内 HSP 数多者，再次前驱下标小者，最后前驱状态按 U < R < 起点的固定次序。分数用 `double`，不逐边取整（与 DuetMap 和 lastz 的 `bigscore` 一致）。
- 主路径终点：所有 `R[i] > 0` 中按同一平局规则取最大者。

### 5.4 次级链（覆盖旁系同源）

**为什么不能在 U/R 上多次回溯。** 设主路径在 query 区间 X 上，P′ 是 X 的另一个拷贝。P′ 中第 k 个 HSP 有两种前驱：

```text
沿 P′ 延续：   R_prim(P′起点) − P_J + (P′ 到 k−1 的段内分) + gain − C_N
从主路径重跳：R_prim(k−1)    − P_J + s_{k−1}             + gain − C_N
```

延续胜出要求 P′ 自身的累积分超过主路径在同一段 query 上的涨分。主路径之所以是主路径，就是因为它在这里更好，所以这个条件通常不成立。于是 P′ 的前驱指针每隔一两个 HSP 就跳回主路径，回溯出来全是碎片。根本原因是 U/R 存的是"全局最优路径"，在路径外的区域它等于"主路径加一小段绕行"。

**做法。** 用 E6 的 L 状态做次级链，回溯规则照搬 minimap2：

1. 先回溯主路径（§7 第 1–3 步），路径上的 HSP 标为已访问。
2. 候选终点：所有 `L[i] ≥ T_sec` 且未访问的 i，按 L 降序（同一平局规则）排序。
3. 依次从每个候选终点沿 L 的前驱回溯，遇到已访问的 HSP 就停；沿途 HSP 标为已访问。
4. 次级链分数 = `L[end] − L[stop]`（stop 为第一个已访问的前驱；走到起点则为 0）。低于 T_sec 丢弃，其 HSP 恢复为未访问。
5. 每段 query 区间限额：统计已保留的次级链中，与它在 query 上的重叠超过它自身 query 长度一半的有几条；达到 N_sec 条就丢弃。
6. 相对主路径的比例（类似 minimap2 的 `pri_ratio`）：要求次级链分数 ≥ r_sec × 主路径**在同一段 query 区间内**的分数。后者用主路径沿 query 的累计分之差求出：回溯主路径时记下每个 HSP 处的累计分，按 query 坐标二分查找区间两端即可。这段 query 区间内没有主路径时，这一条不限制。不和整段主路径的分数比，因为一段主路径可以长达 Mb、分数上百万，段内一个几 kb 的片段重复拷贝几乎不可能达到它的 r_sec 倍。

**仍会丢的。** 分数低于 T_sec 或 r_sec 的拷贝（由门槛决定，可调）；种子阶段就没有产生 HSP 的拷贝（软屏蔽、动态屏蔽、`--queryhspbest`），这与 chaining 无关。query 端多拷贝（多段 query 对应同一 target 位点）不需要次级链，主路径通过回跳 D 本来就能重复访问同一个 target 位点。

**代价。** 次级链上的 HSP 都要做 gapped 延伸。留下的越多，速度越接近 baseline，§0 说的 1.5× 上限还会缩水。T_sec、N_sec、r_sec 就是用来控制这一点的。

## 6. 候选搜索启发式

### 6.1 正常边：K-d 树（已选定方案 A）

- 每个 (p, σ) 一棵 2-d 树，坐标为 j 的**出口** (qe, σ·tout)，两者沿行走方向都递增。
- 对 i 做 box 查询 `[qb_i - maxGap, qb_i + ov_max_i] × [σtin_i - maxGap, σtin_i + ov_max_i]`，它恰好对应 `-ov_max_i ≤ dq, dt ≤ maxGap`。取出候选 j 后再按 §5.1 的完整条件判断是否属于 `Normal(i)`。
- 不用入口坐标 (qb, σ·tin)：由 `dq = qb_i − qe_j ≤ maxGap` 只能推出 `qb_j ≥ qb_i − maxGap − len_j`，按 `qb_i − maxGap` 截断会漏掉长 HSP 这类前驱，与"精确"的说法矛盾。
- box 里可能有还没处理的 j（qb_j ≥ qb_i），它们的值是 -∞，不会被选中；同一 qb 批次内的 j 由 `qb_j < qb_i` 条件排除。
- 节点保存子树内两个最大值：`maxN = max(U,R)`（供 E3）和 `maxL`（供 E6）。因为 `gain ≤ s_i` 且 `C_N ≥ 0`，子树满足 `maxN + s_i ≤ 当前 RN 最优` **且** `maxL + s_i ≤ 当前 L 最优`时整棵跳过。E3 和 E6 共用这一次遍历。
- 未计算的节点值为 -∞；节点 i 算完后向上更新两个最大值（对应 `chain.c` 的 `propagate_max_score`）。
- 复用 `chain.c` 的建树和划分代码，改动是：按 (p,σ) 建多棵、坐标轴改为 (qe, σ·tout)、节点值改为两个 double、查询改为 box 查询。`best_predecessor` 的遍历不照抄（参数传反的缺陷见 §1，而且它的剪枝依赖线性 diag/anti 罚分，不适用于 g(d)）。
- 未采用的方案 B（DuetMap 式下标窗口 W）：会漏候选，HSP 稀疏的区域 W 个下标对应的 query 距离过长，重复区密集时又过短，在 WGA 里没有合理的默认值。

### 6.2 跳转边：界剪枝池（精确）

跳转边不允许重叠（§4），所以 `gain = s_i`；又因为 `P_J ≥ J`，候选 j 的值有上界：

```text
cand(j,i) = R[j] + s_i - P_J(j,i) ≤ R[j] + s_i - J
```

- 池用按 R 降序的有序结构保存所有已释放且满足 `L[j] ≥ T_src`（E5）的 R 状态。从顶部开始逐个算出真实的 `cand(j,i)`，同时记下目前最好的值 best；一旦 `R[j] + s_i - J ≤ best` 就停止，后面的 j 的上界都不超过 best。这在模型内是精确的，DuetMap 的固定 H=500 是有损的。
- 停止条件和"目前最好的真实候选"比，不和"池顶 R\*"比。因为 `ov = 0`，入池的 j 对当前批次的每个 i 都合法，所以按 R\* 推界（`R[j] < R* − λL` 即停）也是精确的；但 best 至少是池顶给出的值，现在的规则不会比按 R\* 推界停得晚，代码也不需要单独维护 slack。
- 推进约束：j 在 `qe_j ≤ qb_i` 时才进池（这同时保证了 `ov = 0` 和 `qe_j < qe_i`）。用一个按 `qe_j` 排序的待释放小顶堆实现，每处理一个新的 qb 批次前，把 `qe_j ≤ qb` 的 j 全部释放入池。
- 安全阀：单个 i 的扫描数超过 `H_max`（默认 1000）就截断，并在 `--stats` 中计数。只有截断发生时结果才是有损的，这一点要在统计里可见。
- 预期：R 是路径的累计分，越靠近当前 query 位置越大，所以 R 落在 best 以下 λL 范围内的通常很少。但任何 HSP 都能以不超过 J+λL 的代价从池顶跳入，所以一段 query 在 target 上命中很多位点时（串联重复、未屏蔽的散在重复），这些 HSP 的 R 都很接近，是退化情形，由 H_max 兜底。
- 同批冻结：与 DuetMap 相同，同一 qb 批次的 HSP 都读取批次开始前的池快照，整批算完再入池，保证输入顺序不影响结果。

### 6.3 其余启发式

| 名称 | 作用 | 是否有损 |
|---|---|---|
| T_solo | 单 HSP 段直接成为 R | 否（属于模型定义） |
| T_src | 跳出资格：`L[j] ≥ T_src`（E5） | 否（L 是独立最优值，不依赖保留了哪条路径） |
| f_ov | 相邻 HSP 最大重叠比例 | 否（属于模型定义） |
| H_max | 跳转扫描上限 | 仅在触发时有损 |
| T_sec / N_sec / r_sec | 次级链的分数门槛、每段 query 区间的条数上限、相对同一 query 区间内主路径分数的比例 | 按定义丢弃，不影响主路径 |
| 次级回溯遇已访问即停 | 一个 HSP 只属于一条链 | 是：后处理的次级链只能拿到被先处理的链剩下的部分（与 minimap2 相同） |

## 7. 回溯流程

1. **选主路径终点**：`e = argmax R[i]`（§5.3 平局规则）。没有正分 R 时，主路径为空，直接进入第 4 步。
2. **回溯主路径**：U、R 状态各存 `(pred, predState, edgeType)`；从 `(e, R)` 出发沿 pred 走到 `no_pred`，依次压入 `(i, state, edgeType)`，最后反转。路径上的 HSP 标为已访问。
3. **切段**：在 edgeType = J 的位置切开，得到若干共线段。每段记录：
   - HSP 下标区间、(p, σ)、query 区间、target 区间。
   - 段分数：段首 `s` 加段内每条正常边的 `gain − C_N`，调用与 E3 相同的代价函数重新累加。
   - 进入方式 `jumpIn`：`S` 起点；`I` 换向（倒位）；`T` 换分区（易位/跨染色体）；`D` 同向回跳（dt < 0，串联重复或重复）；`G` 同向长间隔（大 indel 或 query/target 特有插入）。
4. **次级链**：按 §5.4 的 1–6 步，沿 L 状态的前驱 `predL` 回溯。每条次级链只有一段，`jumpIn = S`，标记为 secondary。
5. **筛 HSP**：主路径和保留的次级链上的 HSP `filter=false`，其余 `filter=true`，调用 `filter_marked_segments`。
6. **两条链各自 gapped**：两张表在 DP 之前已经拆开（§8.2 第 5 步），DP 后直接沿用现有流程：反链表在 `anchors`、正链表在 `secondaryAnchors`，先对反链调用 `finish_one_strand`，`swap_anchor_sets` 后再对正链调用（`reduce_to_points` → `gapped_extend` → 过滤 → 输出）。`finish_one_strand` 里会再做一次 `merge_segments`，对已合并的单链表无副作用。nlchain 打开时 `finish_one_strand` 不进入自身的 `try_reduce_to_chain`（两者互斥）。
7. **（可选）段表输出**：`--nlchainreport=<file>` 写出 TSV（primary/secondary、链号、段号、jumpIn、p、σ、q/t 区间、段分数、HSP 数），用于诊断、断点分析，以及按主链/全部链分别评测。
8. **gapped 之后重新串链（方案 A，2026-11-02）**：每次 `finish_one_strand` 中 gapped 比对过滤之后、插值和输出之前，把该链方向的比对拆成无 gap 块，按 axtChain 的模型重新串链，删掉链分低于 T_chain（`--nlchainclean`，默认 5000）的链，再在被删处切开原比对输出。实现在 `src/nlrechain.c`，算法见 `docs/nlchain_rechain_plan.md` §2，实现记录见 `docs/nlchain_progress.md` §9。

**已知问题，待定**：gapped 延伸从断点两侧的段各自出发，y-drop 可能越过断点，相邻块在 query 上出现重叠。候选处理方式：按段的 query 区间（加一个边距）裁剪输出，或者在输出前后处理。是否需要处理，等 §9 第 3 步看到实际重叠量再定。

## 8. 接口

### 8.1 新文件 `src/nlchain.h` / `src/nlchain.c`（C，沿用 lastz 风格）

```c
typedef struct nlgaptable		// g(d) 的分段线性表（axtChain linearGap 的 qGap 行，§5.1）
	{
	u32			len;
	const unspos* position;		// 严格递增的 d
	const double* cost;			// 对应的 g(d)
	} nlgaptable;

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

typedef struct nlsegment
	{
	u32			first, count;		// 在该链方向压缩后的 segtable 中的下标区间
	u32			part1;				// target 分区下标
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

// stMinus / stPlus：反链 / 正链的单链 HSP 表，各自已 merge、已带分数；反链 pos2 为反向互补坐标。
// 返回主路径分（未缩放）；两张表各自原地压缩为主路径和次级链的 HSP，out 按 (chainId, qBeg) 列出各段。
double nl_reduce_to_chains (seq* target, seq* query,
                            segtable* stMinus, segtable* stPlus,
                            const nlchainparams* params, nlchainresult* out);
void   nl_free_result      (nlchainresult* out);
```

内部状态（字段名与 §5 的记号一一对应，见 §8.4）：

```c
typedef struct nlpath			// U 或 R 的一个状态
	{
	double		score;			// U[i] / R[i]；未达为 -inf
	u32			pred;			// 前驱 HSP；noPred 表示起点
	u32			segHsps;		// 当前段内 HSP 数（平局规则用）
	u8			predState;		// 前驱状态：stateU / stateR
	u8			edge;			// 进入 i 的边：edgeStart / edgeNormal / edgeJump
	} nlpath;

typedef struct nlstate			// 一个 HSP 的全部 DP 状态
	{
	nlpath		U, R;			// E1–E4
	double		L;				// E6
	u32			predL;			// L 的前驱；noPred 表示起点
	} nlstate;
```

`nlstate` 每个 HSP 约 64 字节。另外还有：`nlhsp`（qb、qe、tin、tout、σ、p、s、ρ，约 40 字节）、K-d 树节点（含两个 double 最大值）、跳转池和待释放堆。合计约每个 HSP 150 字节，一条 query 染色体 2M 个 HSP 约 300 MB。若 `unspos` 为 64 位（大基因组编译选项），还要再多一些。

### 8.2 `lastz.c` 改动点

1. `lastzParams` 增加 `int nlchain; nlchainparams nl;`，并设置默认值。
2. CLI：
   - `--nlchain`：开启。
   - `--nlchainjump=<J>[,<λ>,<D0>,<L>]`
   - `--nlchaingap=<maxGap>,<maxDiag>[,<R>]`
   - `--nlchaingapcost=medium|loose`：选 g(d) 的内置表（§5.1）。先不支持从文件读表，需要时再加。
   - `--nlchainsolo=<T_solo>`、`--nlchainsource=<T_src>`、`--nlchainoverlap=<f_ov>`、`--nlchainscan=<H_max>`
   - `--nlchainsecondary=<T_sec>[,<N_sec>,<r_sec>]`（`T_sec=0` 关闭次级链）、`--nlchainreport=<file>`
3. 互斥检查：
   - 与 `--chain`、`--segments` 互斥。
   - 与所有打开 `hspImmediate` 的选项互斥（`--anyornone`/`--stopafterone`、`--limitperquery=`/`--stopafter=`）：这些选项下 HSP 不收集进表，而是立即做 gapped 延伸。
   - 要求两条链都搜索（`whichStrand > 0`），且不带 `[chores=...]`。`--strand=plus/minus` 和 chore 都会跳过反链代码块（`lastz.c:1623`、`1632`），调用点不会执行，会静默退化成不做 chaining，所以直接报错。
   - query 带 `[multiple]` 分区时按 query 分区分别运行（lastz 是否允许 query 用 `[multiple]` 还要核实）。target 不是 `[multiple]` 且有多条 target 序列时给出警告（§3 最后一行的原因）。
4. `collectHspsFromBoth` 的条件加上 `|| currParams->nlchain`，并让 `collectHspsSeparately = false`。
5. 调用点在反链 `start_one_strand` 之后、反链 `finish_one_strand` 之前（`lastz.c` 1677–1681 行之间）：
   1. 现有的 `split_anchors(query->revCompFlags)` 照常执行：反链留在 `anchors`，正链移到 `secondaryAnchors`。
   2. `mergeAnchors` 为真时，对两张表**分别**调用 `merge_segments`。不能在拆表前对合表调用，原因见 §2。
   3. 调用 `nl_reduce_to_chains(target, query, anchors, secondaryAnchors, ...)`。

   此时 query 处于反向互补状态，但 nlchain 只用坐标和已有分数，不读序列，所以不受影响。之后的反链 `finish_one_strand` → `swap_anchor_sets` → `rev_comp_sequence` → 正链 `finish_one_strand` 不变。
6. `finish_one_strand` 中 `if (currParams->chain)` 不变，nlchain 时不进入。
7. `--stats`：增加 HSP 数、主路径段数、各类 jumpIn 计数、次级链数及其 HSP 数、H_max 截断次数、DP 用时（`dbg_timing` 新增 `debugClockNlChaining`）。
8. 不加 `--nlchain` 时，代码路径与现在完全一致，输出逐字节不变。
9. （方案 A，2026-11-02）`lastzParams` 增加 `nlrechainparams nlr` 和 `nlrReportFilename`；CLI 增加 `--nlchainclean=<T_chain>`、`--nlchaincleanreport=<file>`；`--nlchaingapcost` 同时选重新串链的三行代价表；`finish_one_strand` 在过滤之后调用 `nl_rechain_aligns`（`T_chain > 0` 时）。`--nlchainclean=0` 时输出与此前的 nlchain 逐字节一致。

### 8.3 默认值（已在 Alignathon simTest 上调定，2026-11-01）

"初值"列是最初按量纲设定的起点；"冻结值"列是 §9 第 5 步在 simTest（10 个物种对，UCSC 远缘参数，hspthresh=2200）上调参后定下的默认值，过程和数据见 `docs/nlchain_tune_testset.md`。16 个测试物种对不参与调参。

| 参数 | 初值 | 冻结值 | 量纲依据 |
|---|---|---|---|
| g(d) 表 | `medium` | 不变（loose 差别 ≤ 0.0003） | axtChain 用于人–鼠的表；远缘用 `loose` |
| R | 0 | 不变 | 与现有 `--chain` 默认值一致 |
| maxGap / maxDiag | 50000 / 10000 | 不变（test set 上不起作用） | 默认参数下 indel 超过约 2.6 kb 时跳转已更便宜（§5.1），maxDiag 主要约束 U 前驱和搜索范围 |
| J | = hspthresh（默认 3000） | 不变 | 跳转进入的新段至少要值一个 HSP |
| λ, D0, L | 500, 1000, 8 | **1000, 3000**, 8 | 跨度项最多加 4000；形式沿用 DuetMap |
| T_solo | = 2 × hspthresh | **= 1.25 × hspthresh** | 必须明显高于 hspthresh，否则每个 HSP 都直接是 R，U 状态失效 |
| T_src | 0（关闭） | 不变 | |
| f_ov | 0.5 | 不变 | 与 minimap2 mask_level 一致 |
| H_max | 1000 | 不变 | |
| T_sec / N_sec / r_sec | 2 × hspthresh / 5 / 0.5 | **0（关闭）/ 1 / 0.9**；需要旁系拷贝时推荐 T_sec = 4 × hspthresh | T_sec 理由同 T_solo：单个刚过门槛的 HSP 不应单独成为次级链。N_sec、r_sec 取 minimap2 的 `-N 5`、`-p 0.8` 的量级，r_sec 放宽到 0.5 以保留分化较大的旁系拷贝 |

T_solo、T_sec **不与 gappedthresh 挂钩**：没给 `--gappedthresh` 时 lastz 令它等于 hspthresh（`lastz.c:9347`），Human–Mouse 用的又是 lastz 默认参数，两者都是 3000，挂钩等于没设门槛。hspthresh 为自适应（`--hspthresh=top...`）时，J、T_solo 中的 hspthresh 改用本次收集到的最低 HSP 分（与 `finish_one_strand` 里 `lowAnchorScore` 的口径相同）。

### 8.4 代码要求：转移方程必须一眼能对上公式

目标：读 `nlchain.c` 的人拿着 §4–§5 就能逐行核对每个转移，不需要先弄懂剪枝和数据结构。DuetMap 的 `relax()` 把三类转移、两种来源、剪枝和平局都塞进一个双重循环，是这里要避免的写法。

1. **公式层与搜索层分开。**
   - 公式层：纯函数，只做算术，不知道 K-d 树、池、批次的存在。
   - 搜索层（K-d 树查询、跳转池、待释放堆、批次冻结）只负责枚举候选 j，对每个 j 调用公式层。剪枝条件只出现在搜索层，并在旁边注释写明它依据的上界（例如"`gain ≤ s_i`、`C_N ≥ 0`，见 §6.1"）。
2. **一个公式一个函数，名字对应记号**，注释第一行写出它实现的式子：
   ```c
   // gain(j,i) = s_i - ρ_i·ov            (§4)
   static double edge_gain        (const nlhsp* j, const nlhsp* i);
   // g(d)：代价表上分段线性插值，d = 0 时为 0   (§5.1)
   static double gap_cost         (const nlgaptable* t, unspos d);
   // C_N = g(d) + R·max(0,min(dq,dt)) + ρ_i·max(0, ov_t - ov)   (§5.1)
   static double normal_cost      (const nlhsp* j, const nlhsp* i, const nlchainparams* p);
   // P_J = J + λ·min(L, log2(1 + D/D0))   (§5.1)
   static double jump_cost        (const nlhsp* j, const nlhsp* i, const nlchainparams* p);
   // D(j,i) 的三种情形   (§5.1)
   static double jump_span        (const nlhsp* j, const nlhsp* i);
   // Normal(i) 的判定   (§5.1)
   static int    is_normal_edge   (const nlhsp* j, const nlhsp* i, const nlchainparams* p);
   ```
3. **每条转移一个函数，按方程编号命名**，函数体就是方程本身：
   ```c
   // (E2) UJ[i] = max_j R[j] + gain(j,i) - P_J(j,i)
   static void relax_jump   (nlstate* st, u32 j, u32 i, ...);
   // (E3) RN[i] = max_j max(U[j],R[j]) + gain(j,i) - C_N(j,i)
   static void relax_normal (nlstate* st, u32 j, u32 i, ...);
   // (E6) L[i]  = max_j L[j] + gain(j,i) - C_N(j,i)
   static void relax_local  (nlstate* st, u32 j, u32 i, ...);
   // (E1)(E4) 初始化起点，所有候选处理完后按 T_solo 定 U/R
   static void init_state   (nlstate* st, u32 i, ...);
   static void settle_state (nlstate* st, u32 i, ...);
   ```
   E2、E3、E6 的候选值先算成一个局部变量 `cand`，再交给统一的 `offer()` 更新。不在一行表达式里混合多个状态。
4. **平局只在一处**：`static int better (const nlpath* a, const nlpath* b)` 实现 §5.3 的平局规则，所有状态更新都经过它。
5. **记号一致**：变量和字段名与文档记号一致（`qb`、`qe`、`tin`、`tout`、`sigma`、`dq`、`dt`、`ov`、`U`、`R`、`L`），不另起缩写。状态类别、边类型用具名常量（`stateU`、`stateR`、`edgeNormal`、`edgeJump`），不用裸 0/1。
6. **不用宏实现转移或代价**。lastz 习惯用宏写调试输出，这可以保留，但只用于 `chain_dbgChaining` 一类的打印。
7. **坐标先换算再计算**：§4 的统一坐标在进入 DP 前一次性算好，存进 `nlhsp`（qb、qe、tin、tout、σ、p、s、ρ）。转移函数里不再出现 `pos2`、`Q - pos2 - len`、`rcf_rev` 这类原始 lastz 字段，正反链的分支只出现在换算函数里。
8. **参考实现与正式实现共用公式层**：单元测试里的穷举参考解只替换搜索层（暴力枚举全部 j），公式层完全复用。这样测试检验的是剪枝和搜索是否丢解；公式本身另用手算的小例子单独测。
9. 符合全局规范：注释用中文，标识符英文；沿用 lastz 的缩进和函数头注释格式。

## 9. 实施与验证步骤

| 步 | 内容 | 验证方式 |
|---|---|---|
| 0 | **不写代码，先估上限**：在现有 baseline_plain 的 Human–Mouse MAF 上贪心选块（分数降序，与已选块重叠 < 50% 才接受），分别做 query 端 1:1 和两端 1:1 两种，用无旁系同源版真值评 recall/precision。**对照组也要在同一真值上重评**：现有 `--chain` 和 baseline 的 MAF 直接拿来在无旁系同源版真值上评即可，不必重跑（0.236 是全真值口径，不能直接比）。同时核实无旁系同源版真值是哪一侧 1:1（预计两端都是） | 上限 recall 明显高于 `--chain` 在同一真值上的 recall 才继续，否则停下讨论。两种侧别的结果用来确定第 4 步统一用哪一侧 |
| 1 | 实现 `nlchain.c` 和单元测试：随机小 HSP 集合（n ≤ 12），穷举所有合法路径作为参考，比较最优分和路径；公式层另用手算例子测试 | U/R/L 分数与穷举完全一致（平局按规则）；H_max 不触发时剪枝结果与暴力扫描完全一致；随机测试用较小的 maxGap，让部分 HSP 的长度超过 maxGap（覆盖 §6.1 的 box 边界），并包含两端同时重叠的 HSP 对（检查只撤回一次）；`gap_cost` 在表的各节点处取值与 axtChain 表一致；人造"主拷贝 + 较弱旁系拷贝"用例里，次级链能完整回收旁系拷贝，而不是碎片；在一段长主路径内部放一个短的片段重复拷贝，检查 r_sec 不会把它丢掉 |
| 2 | 接入 `lastz.c` | 不加 `--nlchain` 时，`test_data` 上的各种输出格式与当前版本逐字节一致；人造序列（插入一个倒位和一个跨序列易位）的段表里出现 `I` 和 `T`，坐标正确（同时核对反链 pos2 的坐标系）；正反链各放一个对角线号相同的 HSP，确认它们没有被合并、链方向标记正确；`--nlchain --strand=plus` 和带 `[chores]` 时报错退出 |
| 3 | 本地小数据功能验证：pseudocat vs pseudopig | 能运行、`--stats` 计数自洽；统计断点处的 query 重叠量，决定 §7 的已知问题是否需要处理 |
| 4 | 远端：Alignathon Human–Mouse，target 用 `[multiple]` 全基因组，每条 query 染色体一个作业 | 与 baseline_plain、`--chain`、LASTZ+axtChain/chainNet 比较 recall/precision/F/CPU/内存。`--chain` 的准确度不必用 `[multiple]` 重跑：`try_reduce_to_chain` 本来就按 target 分区分批，结果与按染色体对运行相同，只有 CPU 时间需要按新运行方式重测。1:1 的侧别三方统一（按第 0 步的结论）：nlchain 主路径是 query 端（按 `wga_eval_plan.md` 的运行方式即 simMouse）1:1，而 chainNet 默认是 target 端的 net，需要 netSwap 后取 query 端，或者三方都再做同一种两端 1:1 过滤。分两组报：只看主路径（对无旁系同源版真值），主路径加次级链（对全真值）；同时报次级链带来的 gapped 延伸时间增量 |
| 5 | 在 test set 上调 §8.3 参数并冻结，然后跑完整 16 个物种对 | 按 `wga_eval_plan.md` 的流程 |

## 10. 风险与未知

- **recall 上限未知**：`--chain` 损失的 0.41 recall 中有多少来自重排/倒位（由主路径救回），多少来自旁系同源（由次级链救回），目前没有数据。第 0 步就是为了回答这个问题。
- **次级链和速度互相制约**：次级链放得越宽，recall 越高，但 gapped 延伸越接近 baseline 的工作量。需要在 test set 上找 T_sec/N_sec/r_sec 的折中，并且速度和准确度要一起报。
- **罚分在不同进化距离上能否通用**：近缘和远缘的 HSP 密度和分数分布差别很大，很可能需要两套参数（`wga_eval_plan.md` 本来就规划了近缘/远缘两套）。
- **重复区假跳转**：T_solo 和 J 是主要防线（T_solo 必须高于 hspthresh 才起作用，见 §8.3）；在重复区分区评测里单独看 precision。
- **target 端重复使用**：主路径可以经 D 回跳重复使用同一个 target 位点。如果无旁系同源版真值是两端 1:1，这部分在主路径评测中会算成假阳性，需要在第 0 步的两种侧别对比里看量有多大。
- **与 chainNet 比**：chainNet 本身就能表示倒位和嵌套，nlchain 的优势只在于它发生在 gapped 延伸之前（省掉最多约 32% 的时间），并且是单个程序完成。如果准确度不如 chainNet，这件事就只剩速度优势，而速度优势本身不大。

## 11. 修订记录

**2026-10-31，按源码评审修订。** 改动和原因如下：

1. **正常边 gap 代价改为 axtChain 式查表**（§5.1、§8.3）：旧的 `O_c + 30·d` 是碱基级罚分，搬到 HSP 层太陡，d≈88 bp 时就超过跳转代价，maxDiag 形同虚设。现在默认用 axtChain `medium` 表，与跳转代价的交点约 2.6 kb。nlchain 因此不再复用 `--chain` 的 G 参数。
2. **重叠只撤回一次**（§4、§5.1）：旧 C_N 的 `ρ_i·max(0,-dt)` 与 gain 里的 `ρ_i·ov` 在两端同时重叠时会重复扣分；改为总撤回 `ρ_i·max(ov, ov_t)`。
3. **T_solo、T_sec 改为 2 × hspthresh**（§8.3）：旧值 = gappedthresh，而默认 gappedthresh = hspthresh，每个 HSP 都达标，U 状态和 T_sec 都会失效。
4. **按链方向分别 merge，接口改为两张单链表**（§2、§4、§7、§8.1、§8.2）：`merge_segments` 不看也不写回 `id`，不能用于双链合表。调用点改为在现有 `split_anchors` 之后、两张表各自 merge 之后。
5. **互斥检查补全**（§8.2）：`--strand=plus/minus`、`[chores]` 会跳过反链代码块，调用点不执行；会打开 `hspImmediate` 的选项不收集 HSP。
6. **K-d 树改用出口坐标 (qe, σ·tout)**（§6.1）：按入口坐标 qb 截断会漏掉长 HSP 前驱。
7. **跳转边禁止 query 重叠，池剪枝改为与目前最优真实候选比**（§4、§5.1、§6.2）：入池条件变成精确的 `qe_j ≤ qb_i`，剪枝不再依赖 ov_max。
8. **r_sec 改为与同一 query 区间内的主路径分数比**（§5.4）：与整段主路径比，会丢掉长段内部的片段重复。
9. **评测口径**（§0、§9、§10）：0.236 是全真值口径，第 0 步的对照组要在无旁系同源版真值上重评；1:1 的侧别（query 端 / target 端 / 两端）三方要统一；`--chain` 的准确度不必用 `[multiple]` 重跑。
10. **记录现有 lastz 缺陷**（§1）：`chain.c:963-964` 中 `best_predecessor` 的参数顺序传反，G 或 R 非 0 时会漏前驱。
11. 小修正：§7 的 `σ·dt` 改为 `dt`；定义 `ov_max_i`；内存估算改为约 150 字节/HSP；行号更正。

**2026-11-01，第 5 步调参后修订 §8.3。** 在 Alignathon simTest（10 个物种对）上调参，详见 `docs/nlchain_tune_testset.md`：

1. λ 改为 1000、D0 改为 3000、T_solo 改为 1.25 × hspthresh：noparalogy 真值上的平均 F 从 0.7060 升到 0.7085，10 个物种对全部提高。
2. 次级链默认关闭（T_sec = 0），N_sec、r_sec 的默认值改为 1、0.9：次级链在每个物种对上都降低 noparalogy F，ancestor F 的收益只有 +0.001 ~ +0.003。T_sec 不再按 hspthresh 自动换算。
3. 其余参数不变；maxGap、maxDiag 在 test set 上不起作用。与 chainNet 的剩余差距来自 gapped 延伸越出段边界，调参解决不了。
