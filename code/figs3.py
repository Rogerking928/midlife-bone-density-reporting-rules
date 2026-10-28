"""第三版圖（第二輪外部審閱）：正文 Figure 1–3，補充 Figure S2–S5（S1 流程圖在 figs.py）。只存 PNG。

顏色語意（全篇固定）：琥珀＝不一致、紅＝骨質疏鬆、藍綠＝停經狀態結構上無法判定、灰＝一致正常、深藍＝模型曲線／中性。
族裔、部位等類別一律用其他顏色（深藍、紫、綠、棕、梅），不借用上面四個語意色。
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle

import figstyle as fs
from figstyle import bare, panel, save, NAVY, RED, AMBER, GREY, FAINT, INK, TEAL, PURPLE, GREEN, BROWN, PLUM
from figs2 import box_text, arrow_between

R1 = json.load(open(f"{fs.ROOT}/out/results.json"))
R2 = json.load(open(f"{fs.ROOT}/out/results2.json"))
R = json.load(open(f"{fs.ROOT}/out/results3.json"))
R4 = json.load(open(f"{fs.ROOT}/out/results4.json"))
P = pd.read_csv(f"{fs.ROOT}/out/women45_55_v3.csv")
P40 = pd.read_csv(f"{fs.ROOT}/out/women40_59_v2.csv")
DARK_AMBER, LIGHT_TEAL = "#8a5a12", "#9fc9d1"
FN_LOW, FN_OP, TM, TS = 0.74, 0.56, 0.86, 0.12
RACE_COL = {3: NAVY, 4: PURPLE, 1: GREEN, 2: BROWN, 5: "#8a8f3c"}
RACE5 = {3: "Non-Hispanic White", 4: "Non-Hispanic Black", 1: "Mexican American", 2: "Other Hispanic", 5: "Other, including multiracial"}
mn = lambda s: s.replace("-", "−")


def ci_h(ax, y, lo, hi, p, color, marker="o", ms=5):
    ax.plot([lo, hi], [y, y], color=color, lw=1.3, solid_capstyle="butt")
    ax.plot(p, y, marker, color=color, ms=ms, mec="white", mew=0.6, zorder=3)


def spread(ys, gap):
    """把標籤的 y 拉開到至少 gap，保持順序。"""
    order = np.argsort(ys); out = np.array(ys, float)
    for _ in range(50):
        moved = False
        for i in range(1, len(order)):
            a, b = order[i - 1], order[i]
            if out[b] - out[a] < gap:
                d = (gap - (out[b] - out[a])) / 2; out[a] -= d; out[b] += d; moved = True
        if not moved:
            break
    return out


# ------------------------------------------------------------------ Figure 1
def fig1():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 6.0))
    top = fig.add_gridspec(1, 1, left=0.02, right=0.98, bottom=0.73, top=0.97)
    bot = fig.add_gridspec(2, 2, width_ratios=[1.1, 0.9], height_ratios=[1.0, 0.30], wspace=0.85, hspace=0.55,
                           left=0.10, right=0.93, bottom=0.07, top=0.64)
    # A
    ax = fig.add_subplot(top[0]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off"); panel(ax, "A", left=True)
    bmd = 0.70; T = (bmd - TM) / TS
    zm, zs = R2["zmean"]["3|50"], R2["zsd"]["3|50"]; Z = (bmd - zm) / zs
    a = box_text(ax, 13, 50, f"Woman aged 50 years,\nnon-Hispanic White\nFemoral neck BMD\n{bmd:.2f} g/cm²")
    b = box_text(ax, 48, 80, "Postmenopausal rule\nT-score vs women aged 20–29 y")
    c = box_text(ax, 48, 20, "Premenopausal rule\nZ-score vs women aged 50–54 y")
    d = box_text(ax, 84, 80, mn(f"T-score {T:.1f}\nLow bone mass"), fc="#f6ead6", ec=AMBER, color=DARK_AMBER)
    e = box_text(ax, 84, 20, mn(f"Z-score {Z:.1f}\nWithin expected range for age"), fc="#eef1f4", ec=GREY)
    for p_, q_ in ((a, b), (a, c), (b, d), (c, e)):
        arrow_between(ax, p_, q_)
    # B：相圖（每一歲做加權 KDE，再與相鄰年齡平均）
    ax = bare(fig.add_subplot(bot[:, 0])); panel(ax, "B")
    yg = np.linspace(0.36, 1.20, 169); bw = 0.03
    D = np.zeros((20, len(yg)))
    for i, a_ in enumerate(range(40, 60)):
        t = P40[P40.RIDAGEYR == a_]
        k = np.exp(-0.5 * ((yg[None, :] - t.DXXNKBMD.values[:, None]) / bw) ** 2) * t.w.values[:, None]
        D[i] = k.sum(axis=0) / k.sum()
    Ds = np.array([D[max(0, i - 1):i + 2].mean(axis=0) for i in range(20)])
    ax.pcolormesh(np.arange(39.5, 60, 1), np.r_[yg - 0.0025, yg[-1] + 0.0025], Ds.T, cmap="Greys", vmin=0, vmax=Ds.max() * 1.6,
                  shading="flat", rasterized=True)
    xs = np.repeat(np.arange(40, 60), 2) + np.tile([-0.5, 0.5], 20)
    zc = {r: np.array([R2["zcut"][f"{r}|{a_ // 5 * 5}"] for a_ in range(40, 60) for _ in (0, 1)]) for r in (3, 4, 1)}
    ax.fill_between(xs, zc[3], FN_LOW, color=AMBER, alpha=0.35, lw=0)
    # 50 歲以上：白人 Z −2.0 低於 T −2.5 → 骨質疏鬆也會被判同齡範圍內
    m50 = xs >= 49.5
    ax.fill_between(xs[m50], zc[3][m50], FN_OP, where=zc[3][m50] < FN_OP, color=RED, alpha=0.45, lw=0)
    ax.plot([39.5, 59.5], [FN_LOW] * 2, color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.plot([39.5, 59.5], [FN_OP] * 2, color=INK, lw=0.8, ls=(0, (1, 1.5)))
    for r, ls, lw in ((3, "-", 1.2), (4, (0, (4, 1.5)), 0.8), (1, (0, (1, 1)), 0.9)):
        ax.plot(xs, zc[r], color=RED, lw=lw, ls=ls)
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0.40, 1.20)
    ax.set_xticks([40, 45, 50, 55]); ax.set_yticks([0.4, 0.6, 0.8, 1.0, 1.2])
    ax.set_xlabel("Age, years"); ax.set_ylabel("Femoral neck BMD, g/cm²")
    # 右緣直接標線名
    labs = [("T −1.0", FN_LOW, INK), ("T −2.5", FN_OP, INK), ("White", zc[3][-1], RED), ("Black", zc[4][-1], RED),
            ("Mex. Am.", zc[1][-1], RED)]
    ys = spread([l[1] for l in labs], 0.036)
    for (s, y0, c_), y1 in zip(labs, ys):
        ax.annotate(s, xy=(59.5, y0), xytext=(60.4, y1), fontsize=8, color=c_, va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=c_, lw=0.5, shrinkA=0, shrinkB=0))
    ax.text(40.2, 0.535, "Red band (≥50 y):\nT osteoporosis yet\nZ within expected range", fontsize=8, color=RED, va="top", ha="left")
    # C：Sankey（股骨頸，百萬人）
    ax = fig.add_subplot(bot[0, 1]); ax.axis("off"); panel(ax, "C")
    X = R["xtab"]["FN"]; L = ["Normal", "Low bone mass", "Osteoporosis"]; Rn = [("within", "Within\nexpected\nrange"), ("below", "Below\nexpected\nrange")]
    tot = X["N"]; gap = 0.06 * tot
    lv = {k: sum(X[f"{k}|{z}"]["mil"]["t"] for z, _ in Rn) for k in L}; rv = {z: sum(X[f"{k}|{z}"]["mil"]["t"] for k in L) for z, _ in Rn}
    ytop = tot + 2 * gap
    def stack(vals, order):
        y, pos = ytop, {}
        for k in order:
            pos[k] = [y - vals[k], y]; y -= vals[k] + gap
        return pos
    lp, rp = stack(lv, L), stack(rv, [z for z, _ in Rn])
    x0, x1, bw_ = 0.33, 0.66, 0.04
    col = {"Normal": GREY, "Low bone mass": AMBER, "Osteoporosis": RED}
    lc = {k: lp[k][1] for k in L}; rc = {z: rp[z][1] for z in rv}
    for k in L:
        for z, _ in Rn:
            v = X[f"{k}|{z}"]["mil"]["t"]
            if v <= 0:
                continue
            ya1, ya0, yb1, yb0 = lc[k], lc[k] - v, rc[z], rc[z] - v; lc[k] -= v; rc[z] -= v
            xa, xb = x0 + bw_, x1; xm = (xa + xb) / 2
            verts = [(xa, ya1), (xm, ya1), (xm, yb1), (xb, yb1), (xb, yb0), (xm, yb0), (xm, ya0), (xa, ya0), (xa, ya1)]
            codes = [MPath.MOVETO] + [MPath.CURVE4] * 3 + [MPath.LINETO] + [MPath.CURVE4] * 3 + [MPath.CLOSEPOLY]
            alpha = 0.55 if (z == "within" and k != "Normal") else 0.30
            ax.add_patch(PathPatch(MPath(verts, codes), fc=col[k], alpha=alpha, lw=0))
    for k in L:
        y0, y1 = lp[k]; ax.add_patch(Rectangle((x0, y0), bw_, y1 - y0, fc=col[k], lw=0))
    lab_y = spread([(lp[k][0] + lp[k][1]) / 2 for k in L], 0.13 * tot)
    for k, y in zip(L, lab_y):
        ax.text(x0 - 0.03, y, f"{k}\n" + (f"{lv[k]:.1f}M" if lv[k] >= 1 else f"{lv[k]:.2f}M"), ha="right", va="center", fontsize=8)
    for z, lab in Rn:
        y0, y1 = rp[z]; ax.add_patch(Rectangle((x1, y0), bw_, y1 - y0, fc="#5b6b7a", lw=0))
    rlab = spread([(rp[z][0] + rp[z][1]) / 2 for z, _ in Rn], 0.20 * tot)
    for (z, lab), y in zip(Rn, rlab):
        ax.text(x1 + bw_ + 0.03, y, f"{lab}\n" + (f"{rv[z]:.1f}M" if rv[z] >= 1 else f"{rv[z]:.2f}M"), ha="left", va="center", fontsize=8)
    ax.text(x0 + bw_ / 2, ytop + 0.03 * tot, "T-score", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.text(x1 + bw_ / 2, ytop + 0.03 * tot, "Z-score", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(-0.05 * tot, ytop + 0.12 * tot)
    # D：骨質疏鬆 → 同齡範圍內，股骨頸 vs 三部位
    ax = bare(fig.add_subplot(bot[1, 1])); panel(ax, "D")
    rows = [("Femoral neck", R["xtab"]["FN"]["op"]["mil"]), ("Lowest of 3 sites", R["xtab"]["LOW3"]["op"]["mil"])]
    for i, (lab, v) in enumerate(rows):
        y = 1 - i
        ci_h(ax, y, v["lo"], v["hi"], v["t"], RED, marker="s", ms=4.5)
        ax.text(v["hi"] + 0.03, y, f"{v['t']:.2f}M", va="center", fontsize=8)
    ax.set_yticks([1, 0]); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlim(0, 1.0); ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0", "0.5", "1.0"]); ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Osteoporosis → within range, millions", fontsize=8.5)
    save(fig, "Figure_1")


# ------------------------------------------------------------------ Figure 2
STAT6 = [("Menstruating", "Menstruating"), ("No periods, menopause", "No periods, natural menopause"),
         ("Bilateral oophorectomy", "Bilateral oophorectomy"), ("Hysterectomy, ovaries conserved", "Hysterectomy, ovaries conserved"),
         ("No periods, other reason", "No periods, other reason"), ("Not ascertained", "Status not ascertained")]
INDET = {"Hysterectomy, ovaries conserved", "No periods, other reason", "Not ascertained"}
COMP_COL = {"Menstruating": "#b7c4d6", "Natural menopause": NAVY, "Bilateral oophorectomy": PURPLE, "Structurally indeterminate": TEAL}


def fig2():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 6.6))
    left = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.22, left=0.10, right=0.47, bottom=0.56, top=0.97)
    right = fig.add_gridspec(1, 1, left=0.60, right=0.98, bottom=0.56, top=0.97)
    bot = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.6], wspace=0.08, left=0.43, right=0.98, bottom=0.11, top=0.44)
    # A：組成
    ax = bare(fig.add_subplot(left[0])); panel(ax, "A")
    C = R["composition"]; x = np.array(C["x"]); base = np.zeros(len(x))
    for k, c in COMP_COL.items():
        v = np.array(C[k]); ax.fill_between(x, base, base + v, color=c, lw=0); base += v
    ax.set_xlim(40, 59); ax.set_ylim(0, 100); ax.set_xticks([40, 45, 50, 55]); ax.set_xticklabels([])
    ax.set_ylabel("Status, %")
    # 標籤放在各層在 x = 57.5 的中點（白字放深色層、深字放淺色層）
    xi = int(np.argmin(abs(x - 57.5))); cum = 0.0; mids = {}
    for k in COMP_COL:
        v = C[k][xi]; mids[k] = cum + v / 2; cum += v
    ax.text(41, 40, "Menstruating", fontsize=8, color=INK, va="center")
    ax.text(58.6, mids["Natural menopause"], "Natural\nmenopause", fontsize=8, color="white", ha="right", va="center")
    ax.text(58.6, mids["Bilateral oophorectomy"], "Oophorectomy", fontsize=8, color="white", ha="right", va="center")
    ax.text(58.6, mids["Structurally indeterminate"], "Indeterminate", fontsize=8, color="white", ha="right", va="center")
    # B：不一致 spline
    ax = bare(fig.add_subplot(left[1])); panel(ax, "B")
    d = R["spline_disc"]
    ax.fill_between(d["x"], d["lo"], d["hi"], color=AMBER, alpha=0.2, lw=0)
    ax.plot(d["x"], d["p"], color=AMBER, lw=1.5)
    A1 = R1["age1"]; ages = sorted(int(a) for a in A1)
    ax.plot(ages, [A1[str(a)]["any"]["p"] for a in ages], "o", color=AMBER, ms=2.2, alpha=0.45)
    ax.set_xlim(40, 59); ax.set_ylim(0, 70); ax.set_xticks([40, 45, 50, 55]); ax.set_yticks([0, 20, 40, 60])
    ax.set_xlabel("Age, years"); ax.set_ylabel("Discordant, %")
    # C：上限與期望值
    ax = bare(fig.add_subplot(right[0])); panel(ax, "C")
    Cn = R["counts"]; bands = [(a, a + 1) for a in range(40, 60, 2)]; xc = np.array([a + 0.5 for a, _ in bands])
    up = [Cn[f"{a}–{b}|upper"] for a, b in bands]; ex = [Cn[f"{a}–{b}|expected"] for a, b in bands]
    ax.bar(xc, [u["t"] for u in up], width=1.6, facecolor="none", edgecolor=TEAL, lw=0.9, label="Upper bound")
    ax.bar(xc, [e["t"] for e in ex], width=1.6, color=TEAL, label="Expected, premenopausal")
    for v, xo, c in ((up, -0.4, TEAL), (ex, 0.4, INK)):
        ax.errorbar(xc + xo, [u["t"] for u in v], yerr=[[u["t"] - u["lo"] for u in v], [u["hi"] - u["t"] for u in v]],
                    fmt="none", ecolor=c, elinewidth=0.7, capsize=1.5, capthick=0.7)
    ax.set_xlim(39, 60); ax.set_xticks([40, 45, 50, 55]); ax.set_ylim(0, 0.62); ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("Age, years (2-year bands)"); ax.set_ylabel("Women, millions")
    ax.legend(loc="upper left", frameon=False, fontsize=8, handlelength=1.0, borderaxespad=0.1)
    # D：狀態 forest
    H = R["hetero"]["Menstrual and surgical status"]
    for j, (key, xmax, ticks, lab, c, mk) in enumerate((("any", 70, [0, 20, 40, 60], "T abnormal,\nZ within range, %", AMBER, "o"),
                                                        ("op", 10, [0, 5, 10], "T osteoporosis,\nZ within range, %", RED, "s"))):
        ax = bare(fig.add_subplot(bot[j]))
        if j == 0: panel(ax, "D", left=True)
        for i, (k, lab_) in enumerate(STAT6):
            y = len(STAT6) - 1 - i; v = H[k][key]
            ci_h(ax, y, v["lo"], v["hi"], v["p"], TEAL if (k in INDET and key == "any") else c, marker=mk)
        ax.set_yticks(range(len(STAT6)))
        ax.set_yticklabels([f"{lab_} (n = {H[k]['any']['n']:,}; {H[k]['share']['p']:.1f}%)" for k, lab_ in STAT6][::-1] if j == 0 else [], fontsize=8)
        ax.set_xlim(0, xmax); ax.set_xticks(ticks); ax.set_ylim(-0.6, len(STAT6) - 0.4)
        ax.set_xlabel(lab, fontsize=8.5)
    save(fig, "Figure_2")


# ------------------------------------------------------------------ Figure 3：T 對 Z 正面比較
def fig3():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 6.0))
    gA = fig.add_gridspec(1, 1, left=0.10, right=0.42, bottom=0.58, top=0.96)
    gB = fig.add_gridspec(1, 1, left=0.66, right=0.97, bottom=0.58, top=0.96)
    gC = fig.add_gridspec(1, 1, left=0.30, right=0.49, bottom=0.08, top=0.44)
    ax = bare(fig.add_subplot(gA[0])); panel(ax, "A", left=True)
    C = R2["frac_curve"]; lo_, hi_ = C["T_range_obs"]; Tg = np.array(C["T"]); m = (Tg >= lo_) & (Tg <= hi_)
    ax.axvspan(lo_, -1, color=AMBER, alpha=0.16, lw=0)
    ax.fill_between(Tg[m], np.array(C["lo"])[m], np.array(C["hi"])[m], color=NAVY, alpha=0.18, lw=0)
    ax.plot(Tg[m], np.array(C["p"])[m], color=NAVY, lw=1.5)
    ax.set_xlim(lo_, hi_); ax.set_ylim(0, 32); ax.set_yticks([0, 10, 20, 30])
    ax.set_xlabel("Femoral neck T-score"); ax.set_ylabel("Prior fracture, age-standardized %")
    ax.text(-0.9, 31, "T below −1.0\n(discordant zone)", ha="left", va="top", fontsize=8, color=DARK_AMBER)
    # B：盛行率比
    ax = bare(fig.add_subplot(gB[0])); panel(ax, "B")
    PR = R["pr_tz"]
    rows = [("T-score, per 1 lower", "T-score, per 1 unit lower", NAVY), ("Z-score, per 1 lower", "Z-score, per 1 unit lower", NAVY),
            ("T low bone mass vs normal", "T-score category: low bone mass vs normal", AMBER),
            ("T osteoporosis vs normal", "T-score category: osteoporosis vs normal", RED),
            ("Z below vs within range", "Z-score: below vs within expected range", NAVY)]
    for i, (lab, k, c) in enumerate(rows):
        y = len(rows) - 1 - i - (0.4 if i >= 2 else 0); v = PR[k]
        ci_h(ax, y, v["lo"], v["hi"], v["pr"], c)
    ax.set_xscale("log"); ax.set_xlim(0.4, 8); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.5", "1", "2", "4", "8"]); ax.minorticks_off()
    ax.plot([1, 1], [-0.9, len(rows) - 0.3], color=INK, lw=0.6)
    ax.set_yticks([len(rows) - 1 - i - (0.4 if i >= 2 else 0) for i in range(len(rows))]); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_ylim(-0.9, len(rows) - 0.3)
    ax.set_xlabel("Age-adjusted prevalence ratio")
    # C：C 統計量（第四版：補 T ≤ −1.0 二分、同陽性率 Z、Z ≤ −1.0）
    ax = bare(fig.add_subplot(gC[0])); panel(ax, "C", left=True)
    A = R4["auc"]; zc = R4["z_matched_cut"]
    groups = [[("T-score, continuous", "T-score, continuous", NAVY), ("Z-score, continuous", "Z-score, continuous", NAVY)],
              [("T-score category (3 levels)", "T-score category (3 levels)", AMBER)],
              [("T-score ≤ −1.0 (2 levels)", "T-score ≤ −1.0", AMBER),
               ("Z-score at matched positivity (2 levels)", "Z-score ≤ " + f"{zc:.2f}".replace("-", "−") + " (same positivity)", NAVY),
               ("Z-score ≤ −1.0 (2 levels)", "Z-score ≤ −1.0", NAVY), ("Z-score ≤ −2.0 (2 levels)", "Z-score ≤ −2.0 (rule)", GREY)]]
    ys, labs, y = [], [], 0.0
    for g in groups:
        for k, lab, c in g:
            v = A[k]; ci_h(ax, -y, v["lo"], v["hi"], v["c"], c)
            ax.text(1.04, -y, f"{v['c']:.3f} ({v['lo']:.3f}–{v['hi']:.3f})", transform=ax.get_yaxis_transform(), va="center", fontsize=8)
            ys.append(-y); labs.append(lab); y += 1
        y += 0.5
    ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=8)
    ax.plot([0.5, 0.5], [-y + 0.9, 0.7], color=INK, lw=0.6)
    ax.set_xlim(0.45, 0.66); ax.set_xticks([0.5, 0.6]); ax.set_xticklabels(["0.50", "0.60"]); ax.set_ylim(-y + 0.9, 0.7)
    ax.set_xlabel("C-statistic for prior fracture")
    fmt = lambda d: f"{d['c']:.3f} ({d['lo']:.3f} to {d['hi']:.3f})".replace("-0", "−0")
    dl = [("Differences (95% CI)", None),
          ("Continuous T vs Z", A["T-score, continuous minus Z-score, continuous"]),
          ("T ≤ −1.0 vs Z, same positivity", A["T-score ≤ −1.0 (2 levels) minus Z-score at matched positivity (2 levels)"]),
          ("Z, same positivity vs Z ≤ −2.0", A["Z-score at matched positivity (2 levels) minus Z-score ≤ −2.0 (2 levels)"])]
    tx = fig.add_axes([0.755, 0.08, 0.24, 0.36]); tx.axis("off"); tx.set_xlim(0, 1); tx.set_ylim(0, 7)
    for i, (lab, d) in enumerate(dl):
        tx.text(0, 6.6 - 1.75 * i, lab if d is None else f"{lab}\n{fmt(d)}", fontsize=8, va="top", fontweight="bold" if d is None else "normal")
    save(fig, "Figure_3")


# ------------------------------------------------------------------ Figure S2：規格曲線
SITE_COL = {"Femoral neck": NAVY, "Total hip": GREEN, "Lumbar spine": BROWN, "Lower of femoral neck and total hip": PURPLE,
            "Lowest of femoral neck, total hip and lumbar spine": PLUM}


def figS2():
    sp = pd.DataFrame(R2["specs"]).sort_values("any").reset_index(drop=True)
    dims = [("site", list(SITE_COL)), ("tref", ["NHANES III cut-points", "2005–2010 reference", "2005–2010 reference (spine)"]),
            ("zref", ["Race-specific, 5-year age bands", "All races, 5-year age bands", "Race-specific, 10-year age bands"]),
            ("weighting", ["Weighted", "Unweighted"]), ("hormone", ["All women", "Never-users only"])]
    labels = [(d, l) for d, ls in dims for l in ls]
    fig = plt.figure(figsize=(fs.WIDTH_IN, 6.4))
    ax = bare(fig.add_axes([0.42, 0.62, 0.56, 0.35])); panel(ax, "A", left=True)
    x = np.arange(len(sp))
    ax.vlines(x, sp.any_lo, sp.any_hi, color=FAINT, lw=1.0)
    ax.scatter(x, sp["any"], s=8, c=[SITE_COL[s] for s in sp.site], zorder=3, lw=0)
    k = int(np.flatnonzero(sp.primary)[0])
    ax.scatter([k], [sp["any"][k]], s=40, facecolors="none", edgecolors=INK, lw=0.9, zorder=4)
    ax.annotate(f"Primary analysis\n({R['spec_primary_pct_rank']:.0f}th percentile)", xy=(k, sp["any"][k] + 1.2), xytext=(k - 18, 55), fontsize=8, ha="right", va="top",
                arrowprops=dict(arrowstyle="-", color=INK, lw=0.6, shrinkA=1, shrinkB=2))
    ax.set_xlim(-1, len(sp)); ax.set_xticks([]); ax.set_ylim(0, 58)
    ax.set_ylabel("T abnormal, Z within range, %")
    ax2 = bare(fig.add_axes([0.42, 0.22, 0.56, 0.37]), keep=())
    for j, (d, l) in enumerate(labels):
        y = len(labels) - 1 - j; on = (sp[d] == l).values
        ax2.scatter(x[on], np.full(on.sum(), y), s=4, c=[SITE_COL[s] for s in sp.site[on]] if d == "site" else INK, marker="s", lw=0)
    ax2.set_xlim(-1, len(sp)); ax2.set_ylim(-0.6, len(labels) - 0.4)
    ax2.set_yticks(range(len(labels))); ax2.set_yticklabels([l for _, l in labels][::-1], fontsize=8)
    ax2.tick_params(axis="y", length=0); ax2.set_xticks([])
    panel(ax2, "B", left=True)
    ax3 = bare(fig.add_axes([0.42, 0.06, 0.40, 0.115])); panel(ax3, "C", left=True)
    an = R["spec_anova"]; keys = [("site", "Site"), ("tref", "T reference"), ("zref", "Z reference"), ("weighting", "Weighting"), ("hormone", "Hormone users"), ("Residual", "Residual")]
    ax3.barh(range(len(keys))[::-1], [100 * an[k] for k, _ in keys], color=NAVY, height=0.6)
    for i, (k, _) in enumerate(keys):
        ax3.text(100 * an[k] + 1.5, len(keys) - 1 - i, f"{100 * an[k]:.1f}%", va="center", fontsize=8)
    ax3.set_yticks(range(len(keys))[::-1]); ax3.set_yticklabels([l for _, l in keys], fontsize=8)
    ax3.set_xlim(0, 100); ax3.set_xlabel("Share of variance across specifications, %", fontsize=8.5)
    save(fig, "Figure_S2")


# ------------------------------------------------------------------ Figure S3：異質性（不含狀態）
def figS3():
    H = R["hetero"]; rows = []
    for sec in ("Race and ethnicity", "Body mass index, kg/m²", "Female hormone use", "Current smoking"):
        rows.append(("head", sec, None))
        for k, v in H[sec].items():
            rows.append(("row", k, v))
    n = len(rows); hgt = 0.205 * n + 0.75
    fig = plt.figure(figsize=(fs.WIDTH_IN, hgt))
    ax = bare(fig.add_axes([0.38, 0.75 / hgt, 0.27, 1 - 1.0 / hgt]))
    for i, (kind, lab, v) in enumerate(rows):
        y = n - 1 - i
        if kind == "head":
            ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5, fontweight="bold"); continue
        ax.text(-0.02, y, f"{lab} (n = {v['any']['n']:,})", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8)
        ci_h(ax, y, v["any"]["lo"], v["any"]["hi"], v["any"]["p"], AMBER, ms=4.5)
        m = v["mil"]
        ax.text(1.03, y, f"{v['any']['p']:.1f} ({v['any']['lo']:.1f}–{v['any']['hi']:.1f})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
        ax.text(1.62, y, f"{m['t']:.2f} ({m['lo']:.2f}–{m['hi']:.2f})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
    ax.text(1.03, n - 0.1, "%, 95% CI", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8, fontweight="bold")
    ax.text(1.62, n - 0.1, "Millions, 95% CI", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8, fontweight="bold")
    p = R["xtab"]["FN"]["primary"]["p"]; ax.plot([p, p], [-0.6, n - 0.5], color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 70); ax.set_xticks([0, 20, 40, 60]); ax.set_ylim(-0.7, n + 0.3)
    ax.set_xlabel("T abnormal, Z within range, %")
    save(fig, "Figure_S3")


# ------------------------------------------------------------------ Figure S4：T 對 Z，依族裔
def figS4():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 3.6))
    ax = bare(fig.add_axes([0.10, 0.13, 0.55, 0.82]))
    for code, lab in RACE5.items():
        t = P[P.RIDRETH1 == code]
        ax.scatter(t["T"], t.Z, s=2.5, c=RACE_COL[code], lw=0, alpha=0.6, rasterized=True, label=lab)
    x0, x1, y0, y1 = -3.6, 5.4, -3.6, 4.2  # 第五輪審閱：涵蓋全部 1,929 點（T 最高 5.0、Z 最高 3.9）
    ax.plot([-1, -1], [y0, y1], color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.plot([-2.5, -2.5], [y0, y1], color=INK, lw=0.7, ls=(0, (1, 1.5)))
    ax.plot([x0, x1], [-2, -2], color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_xticks(range(-3, 6)); ax.set_yticks(range(-3, 5))
    ax.set_xlabel("Femoral neck T-score"); ax.set_ylabel("Femoral neck Z-score")
    fig.legend(*ax.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.68, 0.95), frameon=False, fontsize=8, markerscale=3,
               title="Z-score reference group", title_fontsize=8)
    fig.text(0.68, 0.42, "Dashed lines: T-score −1.0\nand Z-score −2.0\nDotted line: T-score −2.5", fontsize=8, va="top")
    save(fig, "Figure_S4")


# ------------------------------------------------------------------ Figure S5：逐歲 vs spline
def figS5():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 2.8))
    ax = bare(fig.add_axes([0.10, 0.16, 0.86, 0.80]))
    A1 = R1["age1"]; ages = np.array(sorted(int(a) for a in A1))
    p = np.array([A1[str(a)]["any"]["p"] for a in ages]); lo = np.array([A1[str(a)]["any"]["lo"] for a in ages]); hi = np.array([A1[str(a)]["any"]["hi"] for a in ages])
    ax.vlines(ages, lo, hi, color=AMBER, lw=1.0, alpha=0.6)
    ax.plot(ages, p, "o", color=AMBER, ms=3.5, label="Single year of age (95% CI)")
    d = R["spline_disc"]
    ax.fill_between(d["x"], d["lo"], d["hi"], color=NAVY, alpha=0.15, lw=0)
    ax.plot(d["x"], d["p"], color=NAVY, lw=1.5, label="Restricted cubic spline (95% CI)")
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0, 80); ax.set_xticks(range(40, 60, 5))
    ax.set_xlabel("Age, years"); ax.set_ylabel("T abnormal, Z within range, %")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    save(fig, "Figure_S5")


if __name__ == "__main__":
    import sys
    for f in (sys.argv[1:] or ["fig1", "fig2", "fig3", "figS2", "figS3", "figS4", "figS5"]):
        globals()[f]()
