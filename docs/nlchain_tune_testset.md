# nlchain 第 5 步：在 Alignathon simTest 上调参

日期：2026-11-01。对应设计 `docs/nlchain_design.md` §8.3、§9 第 5 步；进度总览见 `docs/nlchain_progress.md`。

## 1. 数据：simTest

Alignathon 的 Simulation test set，与哺乳类集用同一棵树（同样 5 个物种，叶到叶距离与哺乳类集的 10 个物种对一一对应）。根基因组是 hg18 chr6 的 10 Mb，每个物种 2 条染色体，约 11 Mb。`wga_eval_plan.md` 规定调参只用这个集，16 个最终测试物种对不参与调参。

| 文件 | 大小 | md5 |
|---|---|---|
| `simTest.seqs.tar.gz` | 18,348,289 B | `579fbbfd2dfc8b07667eb0047d8a90da` |
| `simTest.ancestor.maf.tar.gz`（取 `ancestor.maf`） | 37,782,565 B | `5c73438007aef3cd255403e5ad4f5483` |
| `simTest.noparalogyMafs.tar.gz`（取 `ancestor.noparalogies.maf`） | 84,974,095 B | `28fcdd2181ca095b66a73d674d355ae9` |

**物种对**（target 取名字里的第一个；距离按树计算，与 `wga_eval_plan.md` 一致）：

| 物种对 | 距离 | ancestor 真值碱基对 | noparalogy 真值碱基对 |
|---|---|---|---|
| Mouse–Rat | 0.176 | 12,634,550 | 9,843,332 |
| Cow–Dog | 0.352 | 10,139,165 | 8,596,264 |
| Dog–Human | 0.361 | 9,934,190 | 8,502,334 |
| Cow–Human | 0.387 | 9,839,928 | 8,311,918 |
| Human–Mouse | 0.501 | 9,996,054 | 7,673,702 |
| Human–Rat | 0.508 | 9,860,122 | 7,711,629 |
| Dog–Mouse | 0.573 | 9,742,637 | 7,365,917 |
| Dog–Rat | 0.580 | 9,610,904 | 7,403,233 |
| Cow–Mouse | 0.599 | 9,658,677 | 7,214,166 |
| Cow–Rat | 0.606 | 9,529,731 | 7,231,306 |

- 投影用 `project_truth.py`，20 份投影的碱基对数都与 `--check` 的逐列独立计数一致。
- noparalogy 真值在 10 个物种对上两侧的重复使用都是 0，即两端 1:1，与第 0 步在 simMammals 上的结论相同。
- **simTest 没有近缘物种对**（最近的是 Mouse–Rat 0.176）。所以这次只定了一套中远缘参数，灵长类距离（0.013–0.037）在 test set 上调不了。

## 2. 方法

- **运行方式**与第 4 步相同：target 为 `<target>.fa[multiple]`，每条 query 染色体一个作业。lastz 的公共参数一律用 UCSC 远缘参数（`--hspthresh=2200 --gappedthresh=6000 --ydrop=3400 --masking=50 --scores=HoxD55.q --ambiguous=iupac --format=maf`）。下文 h = hspthresh = 2200。
- **参数显式写全**：每个 nlchain 参数组都写全全部参数，不依赖二进制默认值。二进制 `bin/lastz_nl` 用的是 7939842。
- **对照**：baseline_plain、`--chain`，以及 LASTZ+chainNet medium/3000（流程与第 4 步的 `chainnet.sh` 相同，MAF↔axt 往返校验 0 差异）。
- **评测**：mafComparator `--samples=10000000 --seed=1`，每组在 noparalogy、ancestor 两份真值上各评一次，用原始输出，不做后过滤。
- **汇总口径**：R、P、F 都取 10 个物种对的算术平均（F 先逐对算，再平均）。CPU 为 20 个作业 user+sys 之和；chainNet 的 CPU 另加链化和建 net 各步的时间。
- **选择规则**（开始前与用户约定）：
  - 主路径按 noparalogy 真值上的平均 F（F_np）选；
  - 次级链原定"ancestor 真值上的平均 F（F_anc）不低于只跑主路径、且大多数物种对不低，才默认开启"。实际结果与这条规则有冲突，最后由用户决定，见 §4。
- **可复现性**：不同阶段里参数相同的组，输出的指标完全一致（如 `nlp2_lam500_solo4400_medium_D01000` 与 `nlp_default`、`nls_T4400_r0.5_N5` 与 `nls_default`）。

## 3. 结果

### 3.1 参照组

| 组 | R_np | P_np | F_np | R_anc | P_anc | F_anc | CPU s | target 重复使用 | query 重复使用 |
|---|---|---|---|---|---|---|---|---|---|
| baseline_plain | 0.759 | 0.539 | 0.622 | 0.693 | 0.623 | 0.645 | 2,864 | 18.8% | 18.1% |
| `--chain` | 0.503 | 0.691 | 0.580 | 0.431 | 0.746 | 0.543 | 1,455 | 3.7% | 3.5% |
| chainNet medium/3000 | 0.7505 | **0.6930** | **0.7205** | 0.6696 | 0.7780 | **0.7191** | 2,884 | 0 | 1.6% |
| nlp_default（§8.3 初值，次级链关） | 0.7508 | 0.6668 | 0.7060 | 0.6709 | 0.7496 | 0.7069 | 1,481 | 4.7% | 3.8% |
| nls_default（§8.3 初值，次级链开） | 0.7519 | 0.6410 | 0.6915 | 0.6804 | 0.7305 | 0.7029 | 1,523 | 7.1% | 6.3% |

- test set 上的 baseline precision 比全量 simHuman–simMouse 高很多（0.539 对 0.228）。原因是基因组小，重复序列少。
- 在 test set 上，nlp_default 落后 chainNet 0.0145（F_np），比全量数据上的 0.004 大。

### 3.2 阶段 1：主路径单参数扫描（次级链关，其余取初值）

| 参数 | 取值 → F_np / F_anc |
|---|---|
| 初值 | 0.7060 / 0.7069 |
| J（初值 h） | 0.5h 0.7053/0.7081；0.75h 0.7058/0.7075；1.5h 0.7054/0.7044；2h 0.7042/0.7018 |
| λ（初值 500） | 0 0.7010/0.7091；250 0.7042/0.7086；1000 **0.7072**/0.7034 |
| D0（初值 1000） | 300 0.7052/0.7060；3000 0.7063/0.7073 |
| T_solo（初值 2h） | 1.25h 0.7063/0.7074；1.5h 0.7062/0.7072；3h 0.7059/0.7067；4h 0.7059/0.7066 |
| gap 表（初值 medium） | loose 0.7063/0.7072 |
| maxGap（初值 50000） | 20000、100000：指标与初值完全相同 |
| maxDiag（初值 10000） | 5000、20000：指标与初值完全相同 |
| f_ov（初值 0.5） | 0.3 0.7061/0.7068；0.8 与初值相同 |

- 所有单参数变化下，F 的变化都不超过 ±0.005。
- 主要的取舍来自 λ 和 J：跳转越贵，precision 越高、recall 越低。所以 noparalogy 真值偏好贵的跳转，ancestor 真值偏好便宜的跳转。

### 3.3 阶段 2：λ × T_solo × gap 表 × D0 组合（3×2×2×2 = 24 组）

| 组合 | F_np | F_anc | F_np 不低于初值的物种对 |
|---|---|---|---|
| **λ=1000, T_solo=1.25h, medium, D0=3000** | **0.7085** | 0.7054 | 10/10 |
| λ=1000, T_solo=1.25h, loose, D0=3000 | 0.7084 | 0.7057 | 10/10 |
| λ=1000, T_solo=2h, loose, D0=3000 | 0.7082 | 0.7054 | 10/10 |
| λ=750, T_solo=1.25h, medium, D0=3000 | 0.7079 | 0.7067 | 10/10 |
| λ=500, T_solo=1.25h, medium, D0=3000 | 0.7068 | 0.7081 | 7/10 |

- 在 λ ≥ 750 时，D0=3000 一致地好于 1000。
- medium 与 loose 的差别不超过 0.0003，所以保留 medium。

### 3.4 阶段 3、4：次级链

**阶段 3**（主路径取初值，T_sec ∈ {2,3,4,6}h × r_sec ∈ {0.5,0.7,0.9} × N_sec ∈ {1,5}）：N_sec=1 一律好于 5，r_sec=0.9 一律最好。T_sec 越低，F_anc 越高、F_np 越低。

**阶段 4**（两个候选主路径 + N_sec=1、r_sec=0.9），下表是相对同一主路径只跑主路径时的变化：

| 主路径 | T_sec | ΔF_np（平均 / 不降的物种对数） | ΔF_anc（平均 / 不降的物种对数） |
|---|---|---|---|
| λ=750 | 2h | −0.0085 / 0 | +0.0024 / 6 |
| λ=750 | 3h | −0.0052 / 0 | +0.0021 / 10 |
| λ=750 | 4h | −0.0032 / 0 | +0.0017 / 10 |
| λ=750 | 6h | −0.0013 / 0 | +0.0008 / 10 |
| λ=1000 | 2h | −0.0093 / 0 | +0.0035 / 7 |
| λ=1000 | 3h | −0.0061 / 0 | +0.0031 / 10 |
| λ=1000 | 4h | −0.0040 / 0 | +0.0026 / 10 |
| λ=1000 | 6h | −0.0020 / 0 | +0.0014 / 10 |

- 次级链在每一个物种对上都降低 F_np。
- 收紧后（T_sec ≥ 3h、N_sec=1、r_sec=0.9），10 个物种对的 F_anc 都会提高，平均 +0.001 ~ +0.003。
- 不开次级链、只把 λ 降到 0，也能得到 F_anc 0.7091（F_np 0.7010），与开次级链的最好结果相当。

### 3.5 和 chainNet 的剩余差距（诊断，不属于评测口径）

最好的主路径组合仍比 chainNet 低 0.012（F_np）：recall 相同（0.752 对 0.751），差距全在 precision（0.671 对 0.693），10 个物种对都是这样。

为了把差距拆开，对 nlp_default 和 chainNet 的输出都做了一次两端 1:1 贪心过滤（第 0 步的 `maf_1to1.py greedy --mode both`），再在 noparalogy 真值上评：

| 组 | 过滤前 R / P | 过滤后 R / P |
|---|---|---|
| nlp_default | 0.7508 / 0.6668 | 0.7294 / 0.6812 |
| chainNet medium/3000 | 0.7505 / 0.6930 | 0.7483 / 0.7025 |

- 过滤后 precision 的差距从 0.0262 缩到 0.0213，碱基重复使用只解释了约 1/5。
- 其余约 4/5 来自比对内容本身。最可能的原因是 nlchain 的 gapped 延伸越出了段边界：第 4 步测得 35% 的块超出所属段的 query 区间，而 anchor 变少以后，y-drop 延伸不再被相邻比对挡住。
- 这些都不是 §8.3 的参数能调的，要改算法（例如按段边界裁剪 gapped 延伸）。
- **后续结论（2026-11-01）**：上面的"越出段边界"推测已被否定，按段框裁剪会让 F 下降。差距来自 nlchain 缺少 chainNet 流程中"拆块重新串链、丢掉低分碎片链"这一步，见 `docs/nlchain_precision_gap.md`。
- 另外，过滤让 nlchain 的 recall 降了 0.021，而 chainNet 只降 0.002。说明 nlchain 重复使用的碱基里，也有一部分是真同源。

### 3.6 冻结后的复核

用改了默认值的源码在服务器上重新构建（gcc 13，`bin/lastz_nl_step5`），在 10 个物种对上只给 `--nlchain` 运行。全部 20 个作业的 MAF（去掉 `#` 注释行）与调参组 `nlp2_lam1000_solo2750_medium_D03000` 逐字节一致。本地在 pseudocat vs pseudopig 上，h=3000 和 h=2200 两种设置下，默认输出也都与显式写全参数的输出一致。

## 4. 决定（方案 A，用户确认）

| 参数 | 原初值 | 冻结值 |
|---|---|---|
| λ（spanCoef） | 500 | **1000** |
| D0（spanScale） | 1000 | **3000** |
| T_solo | 2 × hspthresh | **1.25 × hspthresh** |
| T_sec | 2 × hspthresh | **0（次级链默认关闭）** |
| N_sec | 5 | **1**（仅在开启次级链时生效） |
| r_sec | 0.5 | **0.9**（仅在开启次级链时生效） |
| J、L、g(d) 表、R、maxGap、maxDiag、f_ov、T_src、H_max | 不变 | 不变 |

- 在 test set 上，F_np 从 0.7060 升到 0.7085，10/10 个物种对都提高；代价是 F_anc 从 0.7069 降到 0.7054。CPU 为 1,533 s，chainNet 为 2,884 s。
- **次级链默认关闭的理由**：按原判据，T_sec ≥ 3h、N_sec=1、r_sec=0.9 的几组本来可以默认开启。但开启后每个物种对的 F_np 都下降（−0.002 ~ −0.006），F_anc 的收益只有 +0.001 ~ +0.003。nlchain 的定位是替代 chainNet 的 1:1 输出，主指标是 F_np，所以默认关闭。
- **需要旁系拷贝时的推荐设置**：`--nlchainsecondary=<4h>,1,0.9`（h=2200 时为 `8800,1,0.9`）。
- **量纲说明**：T_solo 和 J 随 hspthresh 缩放；λ 和 D0 是绝对值，按 hspthresh=2200 定。lastz 默认的 hspthresh=3000 没有在 test set 上调过。

## 5. 未解决

1. **与 chainNet 的 F_np 差距（test set 上 0.012，全量 simHuman–simMouse 上 0.004）调参补不上**。根因与方案见 `docs/nlchain_precision_gap.md`：不是越出段边界，而是缺少 gapped 后的"重新串链 + 链分下限"；推荐方案 A 尚待确认。
2. **近缘参数**：test set 里没有近缘物种对，灵长类距离要么沿用这套参数，要么另找调参数据。
3. **README.lastz.html 和 `lastz --help`** 还没有写入 nlchain 选项。

## 6. 复现

服务器 `MSA_Dedicated_Server_Remote:/mnt/sda/liuhaodong/workspace/experiment/lastz_profiling`，设置 `PROJECT_ROOT` 为该目录：

```
bash scripts/fetch.sh                      # config/data_manifest.tsv 新增 simTest 三项
bash scripts/prepare_test.sh               # 基因组 + 20 份投影真值 → data/derived/simTest/
python3 -I scripts/gen_tune_groups.py > config/tune_groups.tsv
THREADS=64 bash scripts/run_tune.sh <组,...> [物种对,...]     # runs/tune/<组>/<物种对>/
bash scripts/chainnet_pair.sh <物种对> medium 3000          # runs/tune/chainnet_medium3000/<物种对>/
THREADS=48 bash scripts/eval_tune.sh <组,...> [物种对,...]  # results/tune/eval/<组>/
python3 -I scripts/tune_report.py [--groups ...] [--per-pair F_np] [--sort np]   # results/tune/long.tsv
bash scripts/diag_1to1.sh <组> <物种对>                      # §3.5 诊断
```

- 各阶段的汇总表：`results/tune/stage1_report.txt`、`stage23_report.txt`、`stage4_report.txt`；
- 全部 组 × 物种对 的长表：`results/tune/long.tsv`；
- 驱动脚本：`logs/tune_stage*.sh`；
- 复核用的二进制是 `bin/lastz_nl_step5`（由未提交的工作树构建，sha256 `209e8a8b…1c4c`）。
