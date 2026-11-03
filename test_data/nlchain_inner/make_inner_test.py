#!/usr/bin/env python3
"""生成 nlchain --inner 的构造测试数据（docs/nlchain_prechain_inner_plan.md §7 第 3 步）。

target = 随机 5 kb + A(3 kb) + G(8 kb) + B(3 kb) + 随机 5 kb
query  = A' + G' + B'，A'/B' 为 A/B 的 90% 一致拷贝（只有替换），G' 与 G 等长
G' 中放两段（其余为与 G 无关的随机序列）：
  X：G 正中一段的反向互补（平衡倒位），长度按 HSP 分数校准到约 6500；
  W：一段约 25 bp 的精确同向拷贝，与 A/B 同对角线，HSP 分数在 2000–2500。
     8 kb 空隙中 X 居中时，W 不可能距 X 和两端都 ≥ 2 kb，这里取距左端约
     1.95 kb、距 X 约 2 kb。

用法：make_inner_test.py <输出目录>
输出：target.fa、query.fa、query_rc.fa、target_multi.fa、query_multi.fa、
      expected.json（构造坐标、校准后的 HSP 分数与跳转代价）
"""

import json
import math
import os
import random
import sys

SEED = 20261103
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# HoxD70（lastz 默认）与默认 x-drop
HOXD70 = {
    "A": {"A": 91, "C": -114, "G": -31, "T": -123},
    "C": {"A": -114, "C": 100, "G": -125, "T": -31},
    "G": {"A": -31, "C": -125, "G": 100, "T": -114},
    "T": {"A": -123, "C": -31, "G": -114, "T": 91},
}
XDROP = 910

# 测试参数（与 check_inner_test.py 一致）
J1, TSOLO1 = 3000, 3750          # 第一次串链：--nlchainjump=3000 --nlchainsolo=3750
INNER = 2000                     # --inner=2000
J2, TSOLO2 = min(J1, INNER // 2), min(TSOLO1, 5 * INNER // 4)
LAMBDA, D0, LCAP = 1000, 3000, 8  # nlchain 默认 λ、D0、L


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, s, rate):
    out = []
    for c in s:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != c]))
        else:
            out.append(c)
    return "".join(out)


def mismatch_base(c):
    # 与 c 配对得分最低的碱基，用来挡住 x-drop 延伸
    return min("ACGT", key=lambda b: HOXD70[c][b])


def xdrop_hsp(t, q, t0, q0, n):
    """[t0,t0+n) 与 [q0,q0+n) 为一段全匹配的核心；按 x-drop 向两侧延伸，
    返回 (分数, 左延长度, 右延长度)。与 lastz 的无 gap 延伸同口径。"""
    core = sum(HOXD70[t[t0 + i]][q[q0 + i]] for i in range(n))
    best = run = 0
    ext_r = 0
    i = 0
    while t0 + n + i < len(t) and q0 + n + i < len(q):
        run += HOXD70[t[t0 + n + i]][q[q0 + n + i]]
        i += 1
        if run > best:
            best, ext_r = run, i
        if run < best - XDROP:
            break
    right = best
    best = run = 0
    ext_l = 0
    i = 1
    while t0 - i >= 0 and q0 - i >= 0:
        run += HOXD70[t[t0 - i]][q[q0 - i]]
        if run > best:
            best, ext_l = run, i
        if run < best - XDROP:
            break
        i += 1
    return core + right + best, ext_l, ext_r


def jump_pen(j, d):
    return j + LAMBDA * min(LCAP, math.log2(1 + d / D0))


def write_fa(path, recs):
    with open(path, "w") as f:
        for name, s in recs:
            f.write(">%s\n" % name)
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    rng = random.Random(SEED)

    lenR, lenA, lenG, lenB = 5000, 3000, 8000, 3000
    R1, A, G, B, R2 = (rand_seq(rng, n) for n in (lenR, lenA, lenG, lenB, lenR))
    A2, B2 = mutate(rng, A, 0.10), mutate(rng, B, 0.10)
    G2base = rand_seq(rng, lenG)

    tA, tG, tB = lenR, lenR + lenA, lenR + lenA + lenG
    qA, qG, qB = 0, lenA, lenA + lenG
    target = R1 + A + G + B + R2

    def build(lx, lw, w):
        g2 = list(G2base)
        m = lenG // 2 - lx // 2
        g2[m:m + lx] = revcomp(G[m:m + lx])
        # X 两侧：反链比对中 query 的 m+lx 位（互补后）对 target 的 m-1 位，
        # query 的 m-1 位（互补后）对 target 的 m+lx 位
        g2[m + lx] = COMP[mismatch_base(G[m - 1])]
        g2[m - 1] = COMP[mismatch_base(G[m + lx])]
        g2[w:w + lw] = G[w:w + lw]
        g2[w - 1] = mismatch_base(G[w - 1])
        g2[w + lw] = mismatch_base(G[w + lw])
        return "".join(g2), m

    w = 1950
    # 校准 W：分数取最接近 2250 的长度，须落在 [2050, 2450]
    # （≥ innerThreshold，< T_solo2，gapped < 3000）
    best = None
    for lw in range(20, 32):
        g2, m = build(70, lw, w)
        sW, wl, wr = xdrop_hsp(target, A2 + g2 + B2, tG + w, qG + w, lw)
        if best is None or abs(sW - 2250) < abs(best[1] - 2250):
            best = (lw, sW)
    lw, sW = best
    if not (2050 <= sW <= 2450):
        sys.exit("W 校准失败：s_W=%d" % sW)

    # 校准 X：分数尽量接近 6500，且落在两档跳转代价之间
    best = None
    for lx in range(56, 80, 2):
        g2, m = build(lx, lw, w)
        query = A2 + g2 + B2
        qrc = revcomp(query)
        Q = len(query)
        # 反链上：target[tG+m, +lx) 对 query_rc[Q-(qG+m+lx), +lx)
        sX, xl, xr = xdrop_hsp(target, qrc, tG + m, Q - (qG + m + lx), lx)
        if best is None or abs(sX - 6500) < abs(best[1] - 6500):
            best = (lx, sX, xl, xr, g2, m)
    lx, sX, xl, xr, g2, m = best
    query = A2 + g2 + B2
    Q = len(query)
    sW, wl, wr = xdrop_hsp(target, query, tG + w, qG + w, lw)

    # 构造中各段的 HSP 范围（target 全局坐标、正链 query 坐标，半开）
    X_t = (tG + m - xl, tG + m + lx + xr)
    X_q = (qG + m - xr, qG + m + lx + xl)     # 反链上的左右延伸在正链 query 上对调
    W_t = (tG + w - wl, tG + w + lw + wr)
    W_q = (qG + w - wl, qG + w + lw + wr)

    # §2.1 的公式：路径 ...W →(I) X →(I) B，与直连 W → B（同对角线，C_N = 0）比较
    D_WX = abs(X_t[1] - W_t[1])          # tin_X = te_X，tout_W = te_W
    D_XB = abs(tB - X_t[0])              # tin_B = tb_B，tout_X = tb_X
    cost_b2 = jump_pen(J2, D_WX) + jump_pen(J2, D_XB)
    cost_a = jump_pen(J1, D_WX) + jump_pen(J1, D_XB)
    if not (cost_b2 < sX < cost_a and sX >= TSOLO2):
        sys.exit("X 校准失败：s_X=%d，J2 档代价 %.0f，J 档代价 %.0f" % (sX, cost_b2, cost_a))

    write_fa(os.path.join(out, "target.fa"), [("target", target)])
    write_fa(os.path.join(out, "query.fa"), [("query", query)])
    write_fa(os.path.join(out, "query_rc.fa"), [("query", revcomp(query))])
    write_fa(os.path.join(out, "target_multi.fa"),
             [("chrRand", rand_seq(rng, 7000)), ("target", target)])
    write_fa(os.path.join(out, "query_multi.fa"),
             [("qRand", rand_seq(rng, 6000)), ("query", query)])

    exp = {
        "Q": Q, "T": len(target),
        "A_t": [tA, tA + lenA], "A_q": [qA, qA + lenA],
        "B_t": [tB, tB + lenB], "B_q": [qB, qB + lenB],
        "X_t": list(X_t), "X_q": list(X_q), "X_len": lx, "s_X": sX,
        "W_t": list(W_t), "W_q": list(W_q), "W_len": lw, "s_W": sW,
        "J2": J2, "T_solo2": TSOLO2,
        "D_WX": D_WX, "D_XB": D_XB,
        "jump_cost_J2": round(cost_b2, 1), "jump_cost_J": round(cost_a, 1),
    }
    with open(os.path.join(out, "expected.json"), "w") as f:
        json.dump(exp, f, indent=1)
    print(json.dumps(exp, indent=1))


if __name__ == "__main__":
    main()
