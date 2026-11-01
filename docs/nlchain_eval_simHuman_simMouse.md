# nlchain 评测：Alignathon simHuman–simMouse

日期：2026-11-01。本文说明所用数据集与评测流程，并整理有旁系同源（ancestor）和无旁系同源（noparalogy）两份真值上的结果。实现细节与过程记录见 `docs/nlchain_progress.md` §6。

## 1. 数据集

**来源**：Alignathon 模拟哺乳类集（simMammals），由 EVOLVER 模拟得到，真实同源关系已知。simHuman 与 simMouse 的距离约为 0.501 替换/位点，属于远缘。原始文件从 [Alignathon 数据页](https://cgl.gi.ucsc.edu/data/alignathon/set_mammal.html) 下载：

| 文件 | 用途 |
|---|---|
| `simMammals.seqs.tar.gz` | 各物种序列 |
| `simMammals.ancestor.maf.gz` | 有旁系真值（对齐到 MRCA，含内部节点；md5 `4bab2832a972a26a9a43af150096295e`） |
| `simMammals.noparalogyMafs.tar.gz` | 无旁系真值，取其中的 `ancestor.noparalogies.maf`（md5 `bffc18321f937a1eac183c946049d190`） |

**序列**

| 物种 | 角色 | 染色体数 | 总长 | 小写（软屏蔽）比例 | N |
|---|---|---|---|---|---|
| simHuman | target | 6 | 190,795,996 bp | 11.1% | 0 |
| simMouse | query | 7 | 198,921,416 bp | 12.1% | 0 |

**两份真值**

服务器上没有 mafFilter，因此用 `scripts/project_truth.py` 把多物种真值投影到 simHuman–simMouse 这一对。对真值块里的每一对（人行，鼠行），只要两行有共同的比对列，就输出一个两行块；不产生物种内的配对。投影得到的碱基对数与逐列独立计数（Σ 人行数 × 鼠行数）一致。

| 真值 | 块数 | 同源碱基对 | 侧别 | 用于评什么 |
|---|---|---|---|---|
| 有旁系（ancestor） | 3,461,384 | 157,142,590 | 多对多：simHuman 有 1810 万个、simMouse 有 903 万个碱基参与多个碱基对 | 含旁系拷贝的完整同源关系 |
| 无旁系（noparalogy） | 1,979,154 | 99,958,410 | 两端严格 1:1：两个物种的每个碱基都只出现一次 | 1:1 的直系同源比对 |

## 2. 比较的方法

所有方法都在同一台服务器上运行。target 为 `simHuman.fa[multiple]`（全基因组），query 每次为一条 simMouse 染色体，共 7 个作业，每个作业单核运行。lastz 的公共参数为 UCSC 远缘参数：

```
--hspthresh=2200 --gappedthresh=6000 --ydrop=3400 --masking=50
--scores=HoxD55.q --ambiguous=iupac --format=maf
```

| 方法 | 说明 |
|---|---|
| baseline | 公共参数，不做链化 |
| `--chain` | lastz 自带的共线链化：gapped 之前只保留一条共线链 |
| nlchain_primary | `--nlchain --nlchainsecondary=0`：只保留主路径 |
| nlchain | `--nlchain`：主路径加次级链，其余参数取默认值 |
| LASTZ+chainNet | 对 baseline 的输出按 UCSC doBlastzChainNet 标准流程链化、建 net：`axtChain -linearGap=medium -minScore=3000`，然后 `chainAntiRepeat → chainMergeSort → chainPreNet → chainNet -minSpace=1 → netSyntenic → netToAxt → axtSort → axtToMaf`。另用 `loose/5000` 做对照。 |

注意：`--masking` 是动态掩蔽，只影响同一次运行中后续的 query。作业拆分方式会改变结果，所以所有方法都必须用同一种拆分方式运行。

## 3. 评测流程

1. **拼接预测**：每个方法 7 个作业的 MAF 拼成一个文件；chainNet 直接用 `chainnet.maf`。三类方法都用**原始输出**，不做任何 1:1 后过滤，也不做 netSwap。
2. **比较**：使用 mafComparator 0.9，参数为 `--samples=10000000 --near=0 --seed=1`。
   - recall：从真值中抽约 1000 万个同源碱基对，看其中多大比例在预测中出现；
   - precision：反过来，从预测中抽样，看多大比例在真值中出现；
   - `--near=0` 要求两端位置完全一致；F = 2PR / (P + R)。
3. **两份真值各评一次**。
   - 无旁系真值两端 1:1，任何一侧的重复使用都会被计为假阳性，适合评价 1:1 输出（nlchain_primary、chainNet）。
   - 有旁系真值把旁系配对也算作真同源，适合评价带次级链的输出（nlchain）。
   - 为便于横向比较，所有方法在两份真值上都评了。
4. **抽样误差**：baseline 在无旁系真值上换种子 2 重评，recall 为 0.70914 / 0.70906，precision 为 0.22765 / 0.22781，误差在 10⁻⁴ 量级。
5. **辅助统计**：用 `scripts/maf_1to1.py stats` 统计预测中每个物种有多少碱基参与多于一个碱基对（重复使用），用于解释 precision 的差异。

脚本都在服务器 `MSA_Dedicated_Server_Remote:/mnt/sda/liuhaodong/workspace/experiment/lastz_profiling/scripts/` 下，分别是 `run_nl.sh`、`chainnet.sh`、`eval_nl.sh`、`nl_report.py`。完整汇总见 `results/runs/<RUN_ID>/report_nl.md`，其中 RUN_ID 为 `nl_simHuman-simMouse_baseline_plain+chain+nlchain+nlchain_primary+nlchain_prof+nlchain_primary_prof_all_20261031T104720Z`。

## 4. 结果

### 4.1 无旁系真值（noparalogy，两端 1:1）

| 方法 | recall | precision | F | 预测碱基对 |
|---|---|---|---|---|
| baseline | 0.709 | 0.228 | 0.345 | 311,373,857 |
| `--chain` | 0.290 | 0.551 | 0.380 | 52,628,036 |
| **nlchain_primary** | 0.700 | 0.511 | **0.591** | 136,738,379 |
| nlchain | 0.701 | 0.453 | 0.550 | 154,668,147 |
| **chainNet medium/3000** | 0.702 | 0.516 | **0.595** | 135,820,463 |
| chainNet loose/5000 | 0.700 | 0.516 | 0.594 | 135,636,016 |

### 4.2 有旁系真值（ancestor，多对多）

| 方法 | recall | precision | F |
|---|---|---|---|
| baseline | 0.642 | 0.324 | 0.431 |
| `--chain` | 0.234 | 0.700 | 0.351 |
| **nlchain_primary** | 0.619 | 0.712 | **0.662** |
| nlchain | 0.628 | 0.637 | 0.632 |
| **chainNet medium/3000** | 0.620 | 0.718 | **0.666** |
| chainNet loose/5000 | 0.619 | 0.717 | 0.665 |

### 4.3 资源与碱基重复使用

CPU 为 7 个作业 user+sys 之和。chainNet 的 CPU = baseline 的 lastz CPU + 链化/建 net 各步（共约 75 s）。"重复使用"指参与多于一个碱基对的碱基，占该方法预测碱基对数的比例。

| 方法 | CPU s | 相对 baseline 提速 | 最慢作业墙钟 s | 峰值内存 GB | simHuman 重复使用 | simMouse 重复使用 |
|---|---|---|---|---|---|---|
| baseline | 55,151 | 1.00× | 22,412 | 1.65 | 54.9% | 54.2% |
| `--chain` | 32,445 | 1.70× | 11,473 | 1.38 | 4.2% | 3.7% |
| nlchain_primary | 32,762 | 1.68× | 11,580 | 1.73 | 4.9% | 3.1% |
| nlchain | 32,484 | 1.70× | 11,681 | 1.73 | 13.9% | 11.9% |
| chainNet medium/3000 | 55,226 | 1.00× | — | 1.65 | 0 | 3.6% |
| chainNet loose/5000 | 55,228 | 1.00× | — | 1.65 | 0 | 3.4% |

chainNet 的 net 建在 target（simHuman）上，所以 simHuman 端的重复使用为 0。

## 5. 解读

1. **nlchain_primary 与 LASTZ+chainNet 基本持平，CPU 只用 1/1.68**。
   - F 值：无旁系真值上 0.591 对 0.595，有旁系真值上 0.662 对 0.666。
   - 差距主要来自 precision，约低 0.005，大于抽样误差，但幅度很小。
   - 两者的 recall 几乎相同。
2. **`--chain` 的 recall 损失主要来自"只留一条共线链"**。在 CPU 几乎相同的情况下，nlchain_primary 在无旁系真值上把 recall 从 0.290 提到 0.700，约 2.4 倍。
3. **两份真值的差别**：
   - 有旁系真值上，各方法的 precision 普遍更高，因为预测里的旁系配对也被算作真；
   - 1:1 类方法的 recall 普遍更低，因为真值里的旁系配对在 1:1 输出中本来就无法全部给出。
   - baseline 的重复使用超过 50%，在两份真值上的 precision 都很低。
4. **次级链在两份真值上都降低了 F**：
   - 无旁系真值：recall 加 0.001，precision 减 0.058，F 减 0.041；
   - 有旁系真值：本应对它有利，但 recall 只加 0.009，precision 减 0.075，F 仍减 0.030。
   - 次级链让 simHuman 端的重复使用从 4.9% 升到 13.9%，新增的配对大多不是真值认可的同源对。按当前默认参数，次级链不值得开启；默认值留待第 5 步在 test set 上决定。
5. **chainNet 的参数不敏感**：medium/3000 与 loose/5000 的 F 只差 0.001。
