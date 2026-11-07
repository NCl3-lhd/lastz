# 实施方案：`--realign` 的后验解码模式（`--mea`）

2026-11-07。依据：`docs/realign_ksw2_eval.md` §7 第 4 条、`docs/baseline_fn_attribution.md` §3e。方案经确认后实施（路线 A、目标 (a)、接口 `--mea`）。

## 1. 动机与预期

- 3e：近距离错位（simHuman–simMouse 上占真值 18.3%）里 99% 是"按打分预测路径不比真值差"（①），分差中位数 −320。
- HoxD55 在均匀背景下的 Karlin-Altschul λ ≈ 0.00889（每分），−320 分的路径权重比约 e^−2.84 ≈ 0.058。只把 Viterbi 换成后验解码，翻不过多少 ①。
- MEA 能多做的两件事：
  1. 后验分散（多条路径分数接近）时，取期望正确对数最多的路径；
  2. 用 θ 直接按 precision/recall 取舍，后验低的碱基对改为不比对。mafComparator 只按碱基对计数，这个目标和评测口径一致。

## 2. 模型

- 框架同 `--realign`：切段规则、两端固定、以主对角线为中心的带都不变。
- 路径分布 P(path) ∝ exp(λ·S(path))，S 与 `banded_global` 的打分相同（`scoring->sub` + `--realign` 的 O、E）。λ = 倍数 × λ_KA，λ_KA 按打分矩阵 ACGT 部分、均匀背景求（二分法）。λ 很大时分布集中到 Viterbi 路径。
- 前向-后向：sum-product 版 Gotoh，状态按最后一列分 M/X/Y；X 由 M、Y 开（×go）或由 X 延伸（×ge），Y 同理；段首 (0,0) 视为 M。
- 数值：
  - 每条全局路径 2·(对角线列数) + (gap 列数) = lt + lq，所以每个对角线列乘 e^{−λ smax}、每个 gap 列乘 e^{−λ smax/2}，分布不变，所有单步权重 ≤ 1；
  - 前向、后向按行缩放，记对数缩放因子；
  - 后验 = 前向 M × 后向 B_M / Z。
- 解码目标：最大化 Σ_匹配列 (P(i,j) − θ)，gap 列记 0。θ = c/(1+c) 等价于最大化 E[TP] − c·E[FP]。
- 带宽：先按原规则做 Viterbi（含加宽），在同一带宽里做前向-后向与 MEA；MEA 路径碰到带边就把 w 改为 2w+16 重算。
- 平局：与 `banded_global` 相同（right：删除 > 插入 > 对角线；left：对角线 > 插入 > 删除）。
- 保留 Viterbi 路径（计入 `mea fallback`）的情况：
  - (lt+1)(2w+1) > 2^26；
  - 某行全部下溢；
  - 前向与后向的总权重相差 > 1e-6（对数）。λ 很大时，与行最大值相差 e^−709 以上的格下溢为 0，两个方向丢掉的质量不同。
- 校验（违反即退出）：MEA 目标值 ≥ Viterbi 路径的目标值（Viterbi 路径在带内，是可行解）。

## 3. 接口

- `--mea=<λ倍数>,<θ>`，要求同时给出 `--realign`；λ 倍数 > 0，0 ≤ θ < 1。gap 模型和平局方向取自 `--realign`。
- 不加 `--mea` 时输出逐字节不变。
- 输出分数仍按 lastz 自身打分重算（同 realign D4）。

## 4. 文件

| 文件 | 改动 |
|---|---|
| `src/realign.h` | `realignparams` 加 `mea`、`meaLambdaScale`、`meaTheta`；统计加 `meaChunks`、`meaWidened`、`meaFallback`、`meaChanged` |
| `src/realign.c` | `mea_lambda_ka`、`mea_prepare`（expSub 表与 go/ge，按参数缓存）、`mea_posteriors`、`mea_path_obj`、`mea_decode`、`mea_chunk`；`realign_one` 在 Viterbi 之后调 `mea_chunk` |
| `src/realign_test.c` | (5) 穷举对照，(6) λ 极限，(6b) MEA 平局，(7) `realign_aligns` 整体性质 |
| `src/lastz.c` | 初始化、`--mea=` 解析、与 `--realign` 的依赖检查、`--help` |

## 5. 步骤与验证

1. 单元测试（本地 VM，gcc 11 `-Werror`）：
   - 长度 ≤ 6 的小例子枚举带内全部全局路径，精确后验与前向-后向误差 ≤ 1e-9；MEA 目标与穷举最优相同；路径消耗两侧、不出带。
   - λ×30、θ=0.5 时，MEA 的对齐碱基对集合与 Viterbi 相同（左右对齐 Viterbi 不同的例子是平局，跳过）。
   - 变异检验。
2. 接入 lastz：`base_tests`；不加 `--mea` 时 37 组选项与 HEAD 逐字节一致；报错用例；冒烟（块数、坐标、去 gap 序列不变，各输出格式，`lastz_D`）。
3. 服务器 simTest：输入 G3a + `--inner=1000 --gappedthresh=1000 --gap=500,30 --realign=400,8`（当前推荐配置）；网格 λ 倍数 ∈ {1, 2, 3} × θ ∈ {0, 0.25, 0.5}，加同批对照。本地冒烟显示 λ×0.5 时后验整体低于 0.5（θ=0.5 几乎去掉全部碱基对），所以没有放进网格。
4. 归因：`fn_shift` 在 simHuman–simMouse 上量错位变化。
5. 收尾：`docs/mea_eval.md`。

## 6. 不包括

- Baum-Welch 估计参数；
- 按后验拆块（只改块内路径，不改块的首尾）；
- 输出每列的后验。
