#!/usr/bin/env python3
"""运行并检查 nlchain --inner 的构造测试（docs/nlchain_prechain_inner_plan.md §7 第 3 步）。

先用 make_inner_test.py 生成数据，再运行：
  check_inner_test.py <lastz> <数据目录> [--base <d9db595 的 lastz>] [--stats <collect_stats 构建的 lastz>]
全部通过时退出码为 0。
"""

import argparse
import json
import os
import subprocess
import sys

PARAMS = ("--nlchain --hspthresh=50000 --gappedthresh=3000"
          " --nlchainjump=3000 --nlchainsolo=3750 --inner=2000").split()
FMT = "--format=general:name1,strand1,zstart1,end1,name2,strand2,zstart2+,end2+,score"
# 第二次串链默认沿用第一次的 J、T_solo（D6 取 a）；X 的分数按 b2 档校准，
# 要让 X 进入主路径需显式给 b2 档的值
B2 = ["--nlchaininner=1000,2500"]

failures = 0


def check(cond, what):
    global failures
    print("%s  %s" % ("PASS" if cond else "FAIL", what))
    if not cond:
        failures += 1


def run(lastz, d, target, query, extra, report=True):
    """返回 (比对行集合, 段表行列表, 去掉 # 行的 stdout)。"""
    rep = os.path.join(d, "_report.tsv")
    if os.path.exists(rep):
        os.remove(rep)
    args = [lastz, os.path.join(d, target), os.path.join(d, query)] + PARAMS + [FMT] + extra
    if report:
        args.append("--nlchainreport=" + rep)
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("lastz 失败：%s\n%s" % (" ".join(args), p.stderr))
    body = [l for l in p.stdout.splitlines() if not l.startswith("#")]
    aligns = set()
    for l in body:
        f = l.split("\t")
        aligns.add((f[0], f[1], int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7])))
    segs = []
    if report:
        with open(rep) as f:
            for l in f:
                if l.startswith("#"):
                    continue
                f2 = l.rstrip("\n").split("\t")
                # query kind chain seg jumpIn target strand qBeg qEnd tBeg tEnd score hsps
                segs.append((f2[0], f2[4], f2[5], f2[6], int(f2[7]), int(f2[8]),
                             int(f2[9]), int(f2[10]), int(f2[12])))
    return aligns, segs, "\n".join(body)


def has_align(aligns, strand, t, q, tname="target", qname="query"):
    return ("%s" % tname, "+", t[0], t[1], qname, strand, q[0], q[1]) in aligns


def overlaps_W(aligns, e):
    # 只覆盖 W 而不属于 A/B 外层比对的行
    w0, w1 = e["W_t"]
    return any(a[2] <= w0 and w1 <= a[3] and a[3] - a[2] < 1000 for a in aligns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lastz")
    ap.add_argument("data")
    ap.add_argument("--base")
    ap.add_argument("--stats")
    o = ap.parse_args()
    d = o.data
    with open(os.path.join(d, "expected.json")) as f:
        e = json.load(f)
    Q = e["Q"]
    X_t, X_q, W_t, W_q = e["X_t"], e["X_q"], e["W_t"], e["W_q"]

    # 1. b2 档（J2 = 1000，T_solo2 = 2500）
    al, segs, body_b2 = run(o.lastz, d, "target.fa", "query.fa", B2)
    exp_segs = [
        ("query", "S", "target", "+", e["A_q"][0], W_q[1], e["A_t"][0], W_t[1], 2),
        ("query", "I", "target", "-", X_q[0], X_q[1], X_t[0], X_t[1], 1),
        ("query", "I", "target", "+", e["B_q"][0], e["B_q"][1], e["B_t"][0], e["B_t"][1], 1),
    ]
    check(segs == exp_segs, "b2：段表为 S + A..W、I - X、I + B，坐标与构造一致")
    check(has_align(al, "-", X_t, X_q), "b2：输出中有 X 的反链比对")
    check(not overlaps_W(al, e), "b2：W 的 gapped 分 < 3000，不输出")
    segs_b2 = segs

    # 2. 降 gapped 门槛后 W 出现
    al, _, _ = run(o.lastz, d, "target.fa", "query.fa", B2 + ["--gappedthresh=2000"])
    check(has_align(al, "+", W_t, W_q), "--gappedthresh=2000：W 出现在输出中")

    # 3. 默认（D6a）：第二次与第一次取同样的 J、T_solo
    al, segs, body_a = run(o.lastz, d, "target.fa", "query.fa", [])
    check(segs == [("query", "S", "target", "+", e["A_q"][0], e["B_q"][1],
                    e["A_t"][0], e["B_t"][1], 3)],
          "默认：只有一段 S + A..B（含 W，3 个 HSP），没有 X")
    check(not has_align(al, "-", X_t, X_q) and not overlaps_W(al, e),
          "默认：输出中没有 X、没有 W")
    segs_a = segs
    al, segs, body = run(o.lastz, d, "target.fa", "query.fa", ["--nlchaininner=3000,3750"])
    check(body == body_a and segs == segs_a, "--nlchaininner=3000,3750：与默认逐字节一致")

    # 4. 只给 J2：T_solo2 沿用第一次的 T_solo
    al, segs, body = run(o.lastz, d, "target.fa", "query.fa", ["--nlchaininner=1000"])
    _, segs2, body2 = run(o.lastz, d, "target.fa", "query.fa", ["--nlchaininner=1000,3750"])
    check(body == body2 and segs == segs2, "--nlchaininner=1000：与 --nlchaininner=1000,3750 逐字节一致")

    # 5. d9db595：旧 tweener
    if o.base:
        al, _, _ = run(o.base, d, "target.fa", "query.fa", [], report=False)
        check(not has_align(al, "-", X_t, X_q), "d9db595：没有 X")
        check(has_align(al, "+", W_t, W_q), "d9db595：W 由旧 tweener 输出")

    # 6. query 整体反向互补
    al_rc, segs, _ = run(o.lastz, d, "target.fa", "query_rc.fa", B2)
    flip = {"+": "-", "-": "+"}
    exp_rc = []
    for k, s in enumerate(reversed(exp_segs)):
        exp_rc.append((s[0], "S" if k == 0 else "I", s[2], flip[s[3]], Q - s[5], Q - s[4],
                       s[6], s[7], s[8]))
    check(segs == exp_rc, "query 反向互补：段表 σ 对调、坐标与正向运行一致")
    al_fw, _, _ = run(o.lastz, d, "target.fa", "query.fa", B2, report=False)
    mapped = set((a[0], a[1], a[2], a[3], a[4], flip[a[5]], Q - a[7], Q - a[6]) for a in al_fw)
    check(al_rc == mapped, "query 反向互补：比对坐标与正向运行一致")

    # 7. target [multiple]（构造在第二条序列上）
    al, segs, body = run(o.lastz, d, "target_multi.fa[multiple]", "query.fa", B2)
    check(segs == segs_b2 and body == body_b2, "target [multiple]：段表与比对与单条 target 一致")

    # 8. query [multiple]（构造在第二条 query 上）
    al, segs, body = run(o.lastz, d, "target.fa", "query_multi.fa[multiple]", B2)
    check(segs == segs_b2 and body == body_b2, "query [multiple]：段表与比对与单独跑这条 query 一致")

    # 9. 自适应 hspthresh（不显式给 --gappedthresh，gapped 门槛取 lowAnchorScore）。
    #    主搜索在自适应下本来就以很低的门槛收集 HSP，W 不加 --inner 也会输出，
    #    所以这里检查的是：加 --inner 只多出 X，其余比对不变，即插值 HSP 没有
    #    改变 gapped 门槛（D9）。lowAnchorScore 本身是否不变用临时调试输出确认过，
    #    见 docs/nlchain_progress.md。
    saved = PARAMS[:]
    PARAMS[:] = [p for p in saved if not p.startswith(("--hspthresh", "--gappedthresh"))] \
                + ["--hspthresh=top50%"]
    al, segs, _ = run(o.lastz, d, "target.fa", "query.fa", B2)
    PARAMS[:] = [p for p in PARAMS if not p.startswith("--inner")]
    al0, _, _ = run(o.lastz, d, "target.fa", "query.fa", [], report=False)
    PARAMS[:] = saved
    check(segs[0][3] == "+" and segs[0][5] == W_q[1] and segs[0][8] >= 2,
          "--hspthresh=top50%：W 在主路径第一段的末尾")
    x_rows = set(a for a in al if a[5] == "-" and (a[6], a[7]) == tuple(X_q))
    check(len(x_rows) == 1 and al - x_rows == al0,
          "--hspthresh=top50%：加 --inner 只多出 X，其余比对与不加 --inner 相同")

    # 10. collect_stats 构建的计数
    if o.stats:
        st = os.path.join(d, "_stats.txt")
        subprocess.run([o.stats, os.path.join(d, "target.fa"), os.path.join(d, "query.fa")]
                       + PARAMS + B2 + [FMT, "--stats=" + st], capture_output=True, check=True)
        txt = open(st).read()
        want = ["nlchain inner windows: 1", "windows prescreened: 0",
                "inner HSPs -/+: 1/1", "pass-1 HSPs dropped: 0", "nlchain HSPs: 4"]
        for w in want:
            check(w in txt, "统计：%s" % w)

    print("failures=%d" % failures)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
