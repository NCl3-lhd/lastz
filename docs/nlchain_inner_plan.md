# `--inner` 插值按串联模式选择串链方式：实施方案

> **已被 `docs/nlchain_prechain_inner_plan.md` 取代，未实施**（2026-11-03）。

日期：2026-11-03。状态：D1–D5 已按推荐确认（2026-11-03），未动代码。

## 1. 目标

`--inner` 插值时，每个窗口内用 7-mer 种子找到的 HSP 要先串链，再做 gapped 延伸（`src/tweener.c` 第 817 行）。串链方式跟随主流程：

| 主流程 | 窗口内串链 | 说明 |
|---|---|---|
| `--chain` | 原版 `reduce_to_chain`，用 `--chain` 的 diag/anti | 行为不变 |
| 都不给 | 原版 `reduce_to_chain`，diag = anti = 0 | 行为不变 |
| `--nlchain` | `nl_reduce_to_chains`，参数见 §3 | 新行为 |

`--chain` 和 `--nlchain` 本来就互斥（`lastz.c` 中 `--nlchain can't be used with --chain`），三种情况不会重叠，因此不需要新增命令行选项。

## 2. 已知事实与预期（实验前须知）

- 窗口只取当前链方向（tweener 在 `finish_one_strand` 里按方向分别调用），所以 nlchain 的倒位跳转（换向）用不上。能起作用的是：
  - 同一方向上的跳转：target 上的局部易位、重复和回跳；
  - g(d) gap 罚分；
  - U/R 状态带来的单 HSP 门槛 T_solo；
  - query 上的部分重叠（f_ov）。
- 现状下，`--nlchain` 时窗口内的 `reduce_to_chain` 两个罚分都是 0（`chainDiag`/`chainAnti` 只能由 `--chain` 设置）。也就是说，窗口内所有在两条序列上都单调的 HSP 都会被串起来，间隔不扣分。
- 换成 nlchain 后：
  - 共线但离得远、或偏离对角线的弱 HSP，会因为 g(d) 被丢掉；
  - 单个 HSP 自成一段要求 s ≥ T_solo；
  - target 上次序互换的 HSP 可以通过跳转边保留下来。

  因此 precision 和 recall 往哪个方向变，要靠 §5 第 4 步的实验确认。
- **本方案不解决插值比对绕开 `--nlchainclean`（T_chain）的问题。** 重新串链（`nl_rechain_aligns`）仍在插值之前执行，插值出的比对不受 T_chain 约束。这个问题如果要解决，属于另一个方案（把重新串链挪到插值之后），与本方案独立，可以叠加。

## 3. 设计决策（推荐项待确认）

| # | 问题 | 推荐 | 备选 |
|---|---|---|---|
| D1 | 什么时候启用 | `--nlchain` 和 `--inner` 同时给就启用，不加新开关 | 新增开关 |
| D2 | 窗口内的 J、T_solo 怎么取 | 未显式指定时：J = innerThreshold，T_solo = ⌊5·innerThreshold/4⌋（与外层"按 hspthresh 推导"的规则同形）；显式给了 `--nlchainjump`/`--nlchainsolo` 时照用 | 一律沿用外层取值 |
| D3 | 其余 nlchain 参数（gap 表、maxGap、maxDiag、f_ov、λ、D0、L、T_src、H_max、次级链） | 全部沿用外层 `currParams->nl`；次级链默认关闭，开了就跟外层一致 | 窗口内固定关闭次级链 |
| D4 | 窗口内的调用会累加 `nlchainStats` | 调用前用 memcpy 保存，调用后恢复（只在 `collect_stats` 下），外层统计不受影响 | 不处理，在文档里说明 |
| D5 | 构造测试的生成脚本放哪里 | 放进仓库 `test_data/nlchain_inner/`，避免像上次那样放在 scratchpad 里丢失 | 放 scratchpad |

D2 的理由：外层的 J 默认等于 hspthresh = 3000。窗口内的 HSP 从 innerThreshold（例如 2000）起算，如果沿用 3000，跳入一个 2000–3000 分的 HSP 总是得不偿失，单 HSP 段也几乎进不来，nlchain 就退化成"只做 g(d) 过滤"。

## 4. 代码改动

共 3 个源文件，不新增源文件，Makefile 不用改（`src/Makefile` 的依赖里已经列了所有头文件）。

### 4.1 `src/tweener.h`

```c
#include "nlchain.h"			// non-linear chaining stuff

alignel* tweener_interpolate (alignel* a, seq* seq1, seq* seq2,
                              ...
                              int scale, chainer connect, u32 windowSize,
                              const nlchainparams* nlParams);
```

新参数放在最后。传 NULL 时走原版 `reduce_to_chain`。

### 4.2 `src/tweener.c`

1. 文件头注释加一句：`nlParams` 不为 NULL 时，窗口内改用非线性串联（见本文档）。
2. `#include "nlchain.h"`（通过 tweener.h 引入即可，不重复写）。
3. 私有全局变量：

   ```c
   static const nlchainparams* nlParams;	// NULL => 原版线性串链
   static segtable* nlNoMinus = NULL;		// 空的反链表（窗口内只有当前方向）
   ```

4. `tweener_interpolate`：
   - 参数表末尾加 `const nlchainparams* _nlParams`；
   - setup 阶段：`nlParams = _nlParams;`，非 NULL 时 `nlNoMinus = new_segment_table (1, 0);`；
   - cleanup 阶段：`free_segment_table (nlNoMinus); nlNoMinus = NULL;`（仅当非 NULL）。
5. 把第 815–818 行的串链段改成：

   ```c
   	// chain：--nlchain 时用非线性串联，否则用原版线性串联
   	// （docs/nlchain_inner_plan.md）

   	if (nlParams == NULL)
   		reduce_to_chain (innerAnchors, diagPen, antiPen, scale, konnekt);
   	else
   		nl_chain_inner_hsps ();
   	sort_segments   (innerAnchors, qSegmentsByPos1);
   ```

6. 新增私有函数 `nl_chain_inner_hsps`：

   ```c
   // nl_chain_inner_hsps--
   //	对窗口内的 HSP 做非线性串联。窗口子序列是从当前方向的序列上截出来的，
   //	HSP 都按正链处理，反链表为空；子序列没有分区。innerAnchors 原地压缩为
   //	主路径（和次级链）的 HSP。

   static void nl_chain_inner_hsps (void)
   	{
   	nlchainresult	nlResult = { 0, 0, NULL };
   #ifdef collect_stats
   	u8				savedStats[sizeof(nlchainStats)];
   	memcpy (savedStats, &nlchainStats, sizeof(nlchainStats));
   #endif // collect_stats

   	nl_reduce_to_chains (tweenSeq1, tweenSeq2, nlNoMinus, innerAnchors,
   	                     nlParams, &nlResult);
   	nl_free_result (&nlResult);

   #ifdef collect_stats
   	memcpy (&nlchainStats, savedStats, sizeof(nlchainStats));
   #endif // collect_stats
   	}
   ```

   - 坐标：`convert_hsp` 在 `sigma > 0` 时直接用 `pos1`/`pos2`；`new_sequence` 把 `partition.p` 置为 NULL（`sequences.c` 第 1690 行），走"无分区"分支，`part2 = 0`。与 tweener 子序列的坐标一致。
   - `nl_reduce_to_chains` 在 HSP 为 0 时直接返回，不用另外判断。
   - 窗口内不写 `--nlchainreport`。
   - `nlchainStats` 是匿名结构体，所以用字节缓冲保存（D4）。

### 4.3 `src/lastz.c`（第 3546 行的插值调用处）

```c
	if ((alignList != NULL) && (currParams->innerThreshold > 0))
		{
		nlchainparams	nlInner;
		nlchainparams*	nlInnerP = NULL;

		// --nlchain 时窗口内也用非线性串联；J、T_solo 未指定时按 innerThreshold
		// 推导，与外层按 hspthresh 推导同形（docs/nlchain_inner_plan.md D2）

		if (currParams->nlchain)
			{
			nlInner = currParams->nl;
			if (nlInner.jumpPen   < 0) nlInner.jumpPen   = currParams->innerThreshold;
			if (nlInner.soloScore < 0) nlInner.soloScore = 5 * currParams->innerThreshold / 4;
			nlInnerP = &nlInner;
			}

		dbg_timing_sub (debugClockInterpolation);
		alignList = tweener_interpolate
		              (..., chainScale, chain_connect_penalty, currParams->innerWindow,
		               nlInnerP);
		dbg_timing_add (debugClockInterpolation);
		}
```

`currParams->nl` 里的 `jumpPen`/`soloScore` 在 `nlchain_both_strands` 中只改局部副本，所以这里读到的 -1 仍然表示"未指定"。

### 4.4 文档

- `docs/nlchain_design.md`：在 §8.2 加一条"`--inner` 窗口内的串链方式"，在 §11 修订记录里加一行。
- `docs/nlchain_progress.md`：新增一节，记录实现、验证结果和实验数字。
- `README.lastz.html`：nlchain 目前还没写进 README（这是已列出的待办），本次不改，留到统一补写时一并加上。

## 5. 步骤与验证

1. **改代码**（§4.1–4.3）。
   验证：在本地 OrbStack VM 上用 `-Werror` 编译通过；打开 `collect_stats` 的构建也要编译通过。
2. **回归：不该变的输出必须逐字节一致**，基准为 d9db595。
   - `make base_tests`；
   - 不加 `--nlchain` 的 22 组选项（其中包括 `--inner`），再加一组 `--chain --inner=2000`：stdout、stderr、退出码都一致；
   - `--nlchain` 不加 `--inner` 的 8 组选项一致（命令行回显已经相同，不需要去掉）。
3. **构造测试，确认新行为**（D5：脚本放在 `test_data/nlchain_inner/`）。
   - 生成一对序列：外层是两段约 5 kb、高相似度的共线块，中间夹约 10 kb 空隙。空隙里放：
     - (i) 两段约 300 bp 的同源片段，它们在 target 上的次序与在 query 上相反；
     - (ii) 一个偏离主对角线几 kb、分数在 innerThreshold 和 T_solo 之间的孤立弱 HSP。
   - 运行参数：`--nlchain --hspthresh=6000 --inner=2000`。外层阈值调高，是为了让 (i) 和 (ii) 只能在插值阶段找到。片段的分数用 `dbgShowInnerHsps` 校准。
   - 预期：
     - 新版：(i) 两段都保留，(ii) 被丢掉；
     - d9db595：(i) 只留一段，(ii) 只要在两条序列上都单调就会被保留。
4. **在服务器上跑 simTest 20 个作业**（不在本地跑），三组：
   - `--nlchain`（现有默认，含 T_chain = 5000）；
   - `--nlchain --inner=2000`，旧二进制 d9db595；
   - `--nlchain --inner=2000`，新二进制。

   比较 F_np、F_anc、precision、recall 和插值耗时。评测流程同 `nlchain_progress.md` §10.3。判断标准：新版的 F 要同时高于另外两组，才算有效。如果 precision 仍然下降，再讨论是否叠加"重新串链挪到插值之后"的方案。
5. **记录结果**：更新 `docs/nlchain_progress.md` 和 `docs/nlchain_design.md`。提交前先给用户看 diff 和 commit message。

## 6. 风险

- 窗口很小（≤ 20 kb），D0 = 3000、λ = 1000 的跨度项最多约 2.7k，跳转代价主要由 J 决定。D2 的取值会直接影响跳转是否发生，可能需要在第 4 步里另加一组 J 的取值来对比。
- `nlchain_dbgChaining` 打开时，窗口内的调用也会输出调试信息，输出会变多。这只影响调试，不处理。
- 插值耗时：每个窗口的 HSP 数不多，nlchain 的 K-d 树开销预计可以忽略。第 4 步会看插值耗时来确认。
