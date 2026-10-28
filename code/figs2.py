"""第二版圖（回應外部審閱）：Figure 1–4 ＋ Figure S2–S4（S1 流程圖沿用 figs.py）。只存 PNG。"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle

import figstyle as fs
from figstyle import bare, panel, save, NAVY, RED, AMBER, GREY, FAINT, INK, TEAL, PURPLE, GREEN

R1 = json.load(open(f"{fs.ROOT}/out/results.json"))
R = json.load(open(f"{fs.ROOT}/out/results2.json"))
P = pd.read_csv(f"{fs.ROOT}/out/women45_55_v2.csv")
P40 = pd.read_csv(f"{fs.ROOT}/out/women40_59_v2.csv")
DARK_AMBER = "#8a5a12"
FN_LOW, FN_OP, TM, TS = 0.74, 0.56, 0.86, 0.12


def box_text(ax, x, y, s, **kw):
    style = dict(boxstyle="round,pad=0.45", fc=kw.pop("fc", "white"), ec=kw.pop("ec", "#5b6b7a"), lw=0.8)
    return ax.text(x, y, s, ha=kw.pop("ha", "center"), va="center", fontsize=kw.pop("fontsize", 8.5), bbox=style, **kw)


def arrow_between(ax, a, b, gap=0.8, color="#5b6b7a"):
    """從框 a 的右緣到框 b 的左緣（量出框的實際邊界，箭頭不進框）。"""
    fig = ax.figure; fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    (ax0, ay0), (ax1, ay1) = inv.transform(a.get_bbox_patch().get_window_extent(r))
    (bx0, by0), (bx1, by1) = inv.transform(b.get_bbox_patch().get_window_extent(r))
    ya, yb = (ay0 + ay1) / 2, (by0 + by1) / 2
    ax.annotate("", xy=(bx0 - gap, yb), xytext=(ax1 + gap, ya),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=9,
                                connectionstyle="arc3,rad=0"))
    assert bx0 - gap > ax1 + gap, "arrow reversed"


# ------------------------------------------------------------------ Figure 1
def fig1():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 5.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.40, 1.0], width_ratios=[1.1, 1.0], hspace=0.30, wspace=0.38,
                          left=0.10, right=0.98, bottom=0.08, top=0.96)
    # A：概念圖
    ax = fig.add_subplot(gs[0, :]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off"); panel(ax, "A", left=True)
    bmd, age = 0.70, 50
    T = (bmd - TM) / TS
    zm, zs = R["zmean"]["3|50"], R["zsd"]["3|50"]; Z = (bmd - zm) / zs
    a = box_text(ax, 15, 50, f"Woman aged {age} years,\nnon-Hispanic White\nFemoral neck BMD\n{bmd:.2f} g/cm²")
    b = box_text(ax, 50, 78, "Postmenopausal rule\nT-score vs women aged 20–29 y")
    c = box_text(ax, 50, 22, "Premenopausal rule\nZ-score vs women aged 50–54 y")
    d = box_text(ax, 85, 78, f"T-score {T:.1f}\nLow bone mass".replace("-", "−"), fc="#f6ead6", ec=AMBER, color=DARK_AMBER)
    e = box_text(ax, 85, 22, f"Z-score {Z:.1f}\nWithin expected range for age".replace("-", "−"), fc="#eef1f4", ec=GREY)
    for p, q in ((a, b), (a, c), (b, d), (c, e)):
        arrow_between(ax, p, q)
    # B：相圖（年齡 × BMD）
    ax = bare(fig.add_subplot(gs[1, 0])); panel(ax, "B")
    yb = np.arange(0.36, 1.26, 0.02)
    H = np.zeros((20, len(yb) - 1))
    for i, a_ in enumerate(range(40, 60)):
        t = P40[P40.RIDAGEYR == a_]
        H[i], _ = np.histogram(t.DXXNKBMD, bins=yb, weights=t.w)
        H[i] /= H[i].sum()
    ax.pcolormesh(np.arange(39.5, 60, 1), yb, H.T, cmap="Greys", vmin=0, vmax=H.max() * 1.25, shading="flat", rasterized=True)
    # 不一致區（非西裔白人參考）：Z 截點 < BMD ≤ 0.74
    xs, zc = [], []
    for a_ in range(40, 60):
        xs += [a_ - 0.5, a_ + 0.5]; zc += [R["zcut"][f"3|{a_ // 5 * 5}"]] * 2
    xs, zc = np.array(xs), np.array(zc)
    ax.fill_between(xs, zc, FN_LOW, color=AMBER, alpha=0.38, lw=0)
    ax.plot(xs, zc, color=RED, lw=1.1)
    ax.plot([39.5, 59.5], [FN_LOW] * 2, color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.plot([39.5, 59.5], [FN_OP] * 2, color=INK, lw=0.8, ls=(0, (1, 1.5)))
    for r, ls in ((4, (0, (4, 1.5))), (1, (0, (1, 1)))):
        z = np.array([R["zcut"][f"{r}|{a_ // 5 * 5}"] for a_ in range(40, 60) for _ in (0, 1)])
        ax.plot(xs, z, color=RED, lw=0.7, ls=ls, alpha=0.8)
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0.36, 1.44)
    ax.set_xticks([40, 45, 50, 55]); ax.set_yticks([0.4, 0.6, 0.8, 1.0, 1.2])
    ax.set_xlabel("Age, years"); ax.set_ylabel("Femoral neck BMD, g/cm²")
    ax.text(40.0, 1.43, "Black dashed: T-score −1.0; black dotted: T-score −2.5", ha="left", va="top", fontsize=8)
    ax.text(40.0, 1.365, "Red: Z-score −2.0 (solid, non-Hispanic White;\ndashed, Black; dotted, Mexican American)", ha="left", va="top", fontsize=8, color=RED)
    # C：Sankey（百萬人）
    ax = fig.add_subplot(gs[1, 1]); ax.axis("off"); panel(ax, "C")
    X = R["xtab2"]["FN"]
    L = [("Normal", "Normal"), ("Low bone mass", "Low bone mass"), ("Osteoporosis", "Osteoporosis")]
    Rn = [("within", "Within\nexpected\nrange"), ("below", "Below\nexpected\nrange")]
    tot = X["N"]["t"]; gap = 0.06 * tot; scale = 1.0
    lv = {k: sum(X[f"{k}|{z}"]["mil"]["t"] for z, _ in Rn) for k, _ in L}
    rv = {z: sum(X[f"{k}|{z}"]["mil"]["t"] for k, _ in L) for z, _ in Rn}
    ytop = tot + gap * 2
    def stack(vals, order):
        y, pos = ytop, {}
        for k in order:
            pos[k] = [y - vals[k], y]; y -= vals[k] + gap
        return pos
    lp, rp = stack(lv, [k for k, _ in L]), stack(rv, [z for z, _ in Rn])
    x0, x1, bw = 0.30, 0.60, 0.035
    col = {"Normal": "#9aa7b3", "Low bone mass": AMBER, "Osteoporosis": RED}
    lcur = {k: lp[k][1] for k in lv}; rcur = {z: rp[z][1] for z in rv}
    for k, _ in L:
        for z, _ in Rn:
            v = X[f"{k}|{z}"]["mil"]["t"]
            if v <= 0:
                continue
            ya1, ya0 = lcur[k], lcur[k] - v; yb1, yb0 = rcur[z], rcur[z] - v
            lcur[k] -= v; rcur[z] -= v
            xa, xb = x0 + bw, x1; xm = (xa + xb) / 2
            verts = [(xa, ya1), (xm, ya1), (xm, yb1), (xb, yb1), (xb, yb0), (xm, yb0), (xm, ya0), (xa, ya0), (xa, ya1)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MPath(verts, codes), fc=col[k], alpha=0.45, lw=0))
    for k, lab in L:
        y0, y1 = lp[k]; ax.add_patch(Rectangle((x0, y0), bw, y1 - y0, fc=col[k], lw=0))
        ax.text(x0 - 0.02, (y0 + y1) / 2, f"{lab}\n{lv[k]:.1f}M", ha="right", va="center", fontsize=8)
    for z, lab in Rn:
        y0, y1 = rp[z]; ax.add_patch(Rectangle((x1, y0), bw, y1 - y0, fc="#5b6b7a", lw=0))
        ax.text(x1 + bw + 0.02, (y0 + y1) / 2, f"{lab}\n{rv[z]:.1f}M", ha="left", va="center", fontsize=8)
    ax.text(x0 + bw / 2, ytop + 0.03 * tot, "T-score", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax.text(x1 + bw / 2, ytop + 0.03 * tot, "Z-score", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    m = X["Low bone mass|within"]["mil"]
    ax.text(0.45, -0.16 * tot, f"Low bone mass → within range:\n{m['t']:.1f} million (95% CI {m['lo']:.1f} to {m['hi']:.1f})",
            ha="center", va="top", fontsize=8, color=DARK_AMBER)
    ax.set_xlim(0, 1); ax.set_ylim(-0.32 * tot, ytop + 0.1 * tot)
    save(fig, "Figure_1")


# ------------------------------------------------------------------ Figure 2
STAT6 = [("Menstruating", "Menstruating"), ("No periods, menopause", "No periods, natural menopause"),
         ("Bilateral oophorectomy", "Bilateral oophorectomy"), ("Hysterectomy, ovaries conserved", "Hysterectomy, ovaries conserved"),
         ("No periods, other reason", "No periods, other reason"), ("Not ascertained", "Status not ascertained")]
INDET_COL = {"Hysterectomy, ovaries conserved": TEAL, "No periods, other reason": PURPLE, "Not ascertained": "#7fa7b0"}


def ci_bar(ax, y, d, color, marker="o", ms=5, key="p"):
    ax.plot([d["lo"], d["hi"]], [y, y], color=color, lw=1.3, solid_capstyle="butt")
    ax.plot(d[key], y, marker, color=color, ms=ms, mec="white", mew=0.6, zorder=3)


def fig2():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 6.2))
    top = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.38, left=0.10, right=0.98, bottom=0.57, top=0.96)
    # A：曲線
    ax = bare(fig.add_subplot(top[0])); panel(ax, "A")
    S = R["spline_age"]
    for k, c, lab, ls in (("disc", AMBER, "Discordant category", "-"), ("indet", TEAL, "Status structurally indeterminate", "-"),
                          ("menstruating", GREY, "Still menstruating", (0, (3, 2)))):
        d = S[k]
        ax.fill_between(d["x"], d["lo"], d["hi"], color=c, alpha=0.16, lw=0)
        ax.plot(d["x"], d["p"], color=c, lw=1.5, ls=ls, label=lab)
    A1 = R1["age1"]; ages = sorted(int(a) for a in A1)
    ax.plot(ages, [A1[str(a)]["any"]["p"] for a in ages], "o", color=AMBER, ms=2.2, alpha=0.45)
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0, 100); ax.set_xticks([40, 45, 50, 55])
    ax.set_xlabel("Age, years"); ax.set_ylabel("Women, %")
    ax.legend(loc="upper right", frameon=False, fontsize=8, handlelength=1.8, borderaxespad=0.1)
    # B：人數
    ax = bare(fig.add_subplot(top[1])); panel(ax, "B")
    AF = R["affected"]; bands = ["40–44", "45–49", "50–54", "55–59"]
    bottom = np.zeros(4)
    for s, c in INDET_COL.items():
        v = np.array([AF[b][s]["t"] for b in bands])
        ax.bar(range(4), v, bottom=bottom, color=c, width=0.62, label={"Not ascertained": "Status not ascertained"}.get(s, s))
        bottom += v
    ax.set_xticks(range(4)); ax.set_xticklabels(bands)
    ax.set_xlabel("Age, years"); ax.set_ylabel("Women with discordant category\nand indeterminate status, millions")
    ax.set_ylim(0, max(bottom) * 1.95)
    ax.legend(loc="upper left", frameon=False, fontsize=8, handlelength=1.0, borderaxespad=0.1)
    t = AF["total 45–55"]
    ax.text(3.35, max(bottom) * 1.32, f"Aged 45–55 y: {t['t']:.2f} million\n(95% CI {t['lo']:.2f} to {t['hi']:.2f})",
            ha="right", va="bottom", fontsize=8, color=INK)
    # C：狀態 forest，兩個軸
    bot = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.6], wspace=0.08, left=0.43, right=0.98, bottom=0.08, top=0.43)
    H = R["hetero"]["Menstrual and surgical status"]
    for j, (key, xmax, ticks, lab, c, mk) in enumerate((("any", 70, [0, 20, 40, 60], "T abnormal,\nZ within range, %", AMBER, "o"),
                                                        ("op", 10, [0, 5, 10], "T osteoporosis,\nZ within range, %", RED, "s"))):
        ax = bare(fig.add_subplot(bot[j]))
        if j == 0: panel(ax, "C", left=True)
        for i, (k, lab_) in enumerate(STAT6):
            y = len(STAT6) - 1 - i
            ci_bar(ax, y, H[k][key], TEAL if (k in INDET_COL and key == "any") else c, marker=mk)
        ax.set_yticks(range(len(STAT6)))
        if j == 0:
            ax.set_yticklabels([f"{lab_} (n = {H[k]['any']['n']:,}; {H[k]['share']['p']:.1f}%)" for k, lab_ in STAT6][::-1], fontsize=8)
        else:
            ax.set_yticklabels([])
        ax.set_xlim(0, xmax); ax.set_xticks(ticks); ax.set_ylim(-0.6, len(STAT6) - 0.4)
        ax.set_xlabel(lab, fontsize=8.5)
    save(fig, "Figure_2")


# ------------------------------------------------------------------ Figure 3
def fig3():
    H = R["hetero"]
    rows = []
    for sec, items in H.items():
        rows.append(("head", sec, None))
        for k, v in items.items():
            rows.append(("row", {"Not ascertained": "Status not ascertained", "No periods, menopause": "No periods, natural menopause"}.get(k, k), v))
    n = len(rows)
    fig = plt.figure(figsize=(fs.WIDTH_IN, 0.205 * n + 0.75))
    ax = bare(fig.add_axes([0.38, 0.75 / (0.205 * n + 0.75), 0.27, 1 - 1.0 / (0.205 * n + 0.75)]))
    for i, (kind, lab, v) in enumerate(rows):
        y = n - 1 - i
        if kind == "head":
            ax.text(-0.02, y, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5, fontweight="bold")
            continue
        ax.text(-0.02, y, f"{lab} (n = {v['any']['n']:,})", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8)
        ci_bar(ax, y, v["any"], AMBER, ms=4.5)
        m = v["mil"]
        ax.text(1.03, y, f"{v['any']['p']:.1f} ({v['any']['lo']:.1f}–{v['any']['hi']:.1f})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
        ax.text(1.62, y, f"{m['t']:.2f} ({m['lo']:.2f}–{m['hi']:.2f})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8)
    ax.text(1.03, n - 0.1, "%, 95% CI", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8, fontweight="bold")
    ax.text(1.62, n - 0.1, "Millions, 95% CI", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8, fontweight="bold")
    p = R["xtab2"]["FN"]["primary"]["p"]
    ax.plot([p, p], [-0.6, n - 0.5], color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 70); ax.set_xticks([0, 20, 40, 60]); ax.set_ylim(-0.7, n + 0.3)
    ax.set_xlabel("T abnormal, Z within range, %")
    save(fig, "Figure_3")


# ------------------------------------------------------------------ Figure 4
def fig4():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 2.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.95, left=0.10, right=0.97, bottom=0.21, top=0.9)
    ax = bare(fig.add_subplot(gs[0])); panel(ax, "A")
    C = R["frac_curve"]; lo_, hi_ = C["T_range_obs"]
    Tg = np.array(C["T"]); m = (Tg >= lo_) & (Tg <= hi_)
    ax.axvspan(max(lo_, -3), -1, color=AMBER, alpha=0.16, lw=0)
    ax.fill_between(Tg[m], np.array(C["lo"])[m], np.array(C["hi"])[m], color=NAVY, alpha=0.18, lw=0)
    ax.plot(Tg[m], np.array(C["p"])[m], color=NAVY, lw=1.5)
    ax.set_xlim(lo_, hi_); ax.set_ylim(0, 30); ax.set_yticks([0, 10, 20, 30])
    ax.set_xlabel("Femoral neck T-score"); ax.set_ylabel("Prior hip, wrist or spine\nfracture, age-standardized %")
    X = R["xtab2"]["FN"]
    ab = sum(X[f"{c}|{z}"]["p"] for c in ("Low bone mass", "Osteoporosis") for z in ("within", "below"))
    share = 100 * sum(X[f"{c}|within"]["p"] for c in ("Low bone mass", "Osteoporosis")) / ab
    ax.text(-0.95, 28.5, f"T-score below −1.0\n({share:.1f}% of these women\nhave Z within range)", ha="left", va="top", fontsize=8, color=DARK_AMBER)
    # B：盛行率比
    ax = bare(fig.add_subplot(gs[1])); panel(ax, "B")
    PR = R["pr"]; rows = []
    for o, ol in (("frac", "Prior fracture"), ("told", "Told of osteoporosis")):
        for k, kl in (("Crude", "Crude"), ("Age-adjusted", "Age-adjusted"), ("Adjusted for age, BMI, race and ethnicity, hormone use", "Fully adjusted")):
            rows.append((f"{ol}, {kl.lower()}", PR[o][k]))
    for i, (lab, v) in enumerate(rows):
        y = len(rows) - 1 - i - (0.5 if i >= 3 else 0)
        ax.plot([v["lo"], v["hi"]], [y, y], color=AMBER, lw=1.3)
        ax.plot(v["pr"], y, "o", color=AMBER, ms=5, mec="white", mew=0.6)
    ax.set_xscale("log"); ax.set_xlim(0.8, 5); ax.set_xticks([1, 2, 4]); ax.set_xticklabels(["1", "2", "4"]); ax.minorticks_off()
    ax.plot([1, 1], [-1, len(rows)], color=INK, lw=0.6)
    ax.set_yticks([len(rows) - 1 - i - (0.5 if i >= 3 else 0) for i in range(len(rows))])
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_ylim(-1.0, len(rows) - 0.3)
    ax.set_xlabel("Prevalence ratio (discordant\nvs concordant normal)")
    save(fig, "Figure_4")


# ------------------------------------------------------------------ Figure S2：規格曲線
def figS2():
    sp = pd.DataFrame(R["specs"]).sort_values("any").reset_index(drop=True)
    dims = [("site", ["Femoral neck", "Total hip", "Lumbar spine", "Lower of femoral neck and total hip", "Lowest of femoral neck, total hip and lumbar spine"]),
            ("tref", ["NHANES III cut-points", "2005–2010 reference", "2005–2010 reference (spine)"]),
            ("zref", ["Race-specific, 5-year age bands", "All races, 5-year age bands", "Race-specific, 10-year age bands"]),
            ("weighting", ["Weighted", "Unweighted"]), ("hormone", ["All women", "Excluded ever-users"])]
    labels = [(d, l) for d, ls in dims for l in ls]
    fig = plt.figure(figsize=(fs.WIDTH_IN, 5.6))
    ax = bare(fig.add_axes([0.40, 0.56, 0.58, 0.40])); panel(ax, "A")
    x = np.arange(len(sp))
    ax.vlines(x, sp.any_lo, sp.any_hi, color=FAINT, lw=1.0)
    ax.scatter(x, sp["any"], s=7, c=np.where(sp.primary, RED, NAVY), zorder=3, lw=0)
    ax.set_xlim(-1, len(sp)); ax.set_xticks([]); ax.set_ylim(0, 60)
    ax.set_ylabel("T abnormal, Z within range, %")
    k = int(np.flatnonzero(sp.primary)[0])
    ax.annotate("Primary analysis", xy=(k, sp["any"][k]), xytext=(k - 14, 54), fontsize=8, color=RED, ha="right",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.6, shrinkA=1, shrinkB=3))
    ax2 = bare(fig.add_axes([0.40, 0.04, 0.58, 0.48]), keep=())
    for j, (d, l) in enumerate(labels):
        y = len(labels) - 1 - j
        on = (sp[d] == l).values
        ax2.scatter(x[on], np.full(on.sum(), y), s=4, c=np.where(sp.primary[on], RED, INK), marker="s", lw=0)
    ax2.set_xlim(-1, len(sp)); ax2.set_ylim(-0.6, len(labels) - 0.4)
    ax2.set_yticks(range(len(labels))); ax2.set_yticklabels([l for _, l in labels][::-1], fontsize=8)
    ax2.tick_params(axis="y", length=0); ax2.set_xticks([])
    panel(ax2, "B")
    save(fig, "Figure_S2")


# ------------------------------------------------------------------ Figure S3：T 對 Z，依參考組上色
def figS3():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 3.6))
    ax = bare(fig.add_axes([0.10, 0.13, 0.55, 0.82]))
    cols = {"Non-Hispanic White": NAVY, "Non-Hispanic Black": AMBER, "Hispanic": GREEN, "Other, including multiracial": PURPLE}
    for r, c in cols.items():
        t = P[P.race4 == r]
        ax.scatter(t["T"], t.Z, s=2.5, c=c, lw=0, alpha=0.6, rasterized=True, label=r)
    ax.plot([-1, -1], [-3.6, 3.6], color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.plot([-2.5, -2.5], [-3.6, 3.6], color=INK, lw=0.7, ls=(0, (1, 1.5)))
    ax.plot([-3.6, 2.6], [-2, -2], color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.set_xlim(-3.6, 2.6); ax.set_ylim(-3.6, 3.6); ax.set_xticks([-3, -2, -1, 0, 1, 2]); ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel("Femoral neck T-score"); ax.set_ylabel("Femoral neck Z-score")
    fig.legend(*ax.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.68, 0.95), frameon=False, fontsize=8, markerscale=3,
               title="Z-score reference group", title_fontsize=8)
    fig.text(0.68, 0.45, "Dashed lines: T-score −1.0\nand Z-score −2.0\nDotted line: T-score −2.5", fontsize=8, va="top")
    save(fig, "Figure_S3")


# ------------------------------------------------------------------ Figure S4：逐歲 vs spline
def figS4():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 2.8))
    ax = bare(fig.add_axes([0.10, 0.16, 0.86, 0.80]))
    A1 = R1["age1"]; ages = np.array(sorted(int(a) for a in A1))
    p = np.array([A1[str(a)]["any"]["p"] for a in ages]); lo = np.array([A1[str(a)]["any"]["lo"] for a in ages]); hi = np.array([A1[str(a)]["any"]["hi"] for a in ages])
    ax.vlines(ages, lo, hi, color=AMBER, lw=1.0, alpha=0.6)
    ax.plot(ages, p, "o", color=AMBER, ms=3.5, label="Single year of age (95% CI)")
    d = R["spline_age"]["disc"]
    ax.fill_between(d["x"], d["lo"], d["hi"], color=NAVY, alpha=0.15, lw=0)
    ax.plot(d["x"], d["p"], color=NAVY, lw=1.5, label="Restricted cubic spline (95% CI)")
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0, 80); ax.set_xticks(range(40, 60, 5))
    ax.set_xlabel("Age, years"); ax.set_ylabel("T abnormal, Z within range, %")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    save(fig, "Figure_S4")


if __name__ == "__main__":
    import sys
    todo = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "figS2", "figS3", "figS4"]
    for f in todo:
        globals()[f]()
