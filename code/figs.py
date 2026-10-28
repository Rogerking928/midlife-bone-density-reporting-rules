"""Figure 1–2 與 Figure S1（流程圖）。只存 PNG（700 dpi），版式見 figstyle.py。"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import figstyle as fs
from figstyle import bare, panel, save, NAVY, RED, AMBER, GREY, FAINT, INK, TEAL

R = json.load(open(f"{fs.ROOT}/out/results.json"))
P = pd.read_csv(f"{fs.ROOT}/out/women45_55.csv")
f1 = lambda v: f"{v:.1f}"


def ci_bar(ax, y, d, color, marker="o", dy=0.0, ms=5.5, label=None):
    ax.plot([d["lo"], d["hi"]], [y + dy] * 2, color=color, lw=1.3, solid_capstyle="butt")
    ax.plot(d["p"], y + dy, marker, color=color, ms=ms, mec="white", mew=0.6, label=label, zorder=3)


# ------------------------------------------------------------------ Figure 1
def fig1():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 3.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.42, left=0.09, right=0.99, bottom=0.15, top=0.93)
    # A：T 對 Z（45–55 歲）
    ax = bare(fig.add_subplot(gs[0])); panel(ax, "A")
    col = np.where(P.Z <= -2, RED, np.where(P.tcat >= 1, AMBER, "#9aa7b3"))
    ax.scatter(P["T"], P.Z, s=2.2, c=col, lw=0, alpha=0.55, rasterized=True)
    xl, yl = (-3.6, 2.6), (-3.6, 3.6)
    ax.set_xlim(*xl); ax.set_ylim(*yl)
    ax.set_xticks([-3, -2, -1, 0, 1, 2]); ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
    ytop = 2.15                                   # 參考線只畫到這裡，上方留給標籤
    ax.plot([-1, -1], [yl[0], ytop], color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.plot([-2.5, -2.5], [yl[0], ytop], color=INK, lw=0.7, ls=(0, (1, 1.5)))
    ax.plot(xl, [-2, -2], color=INK, lw=0.7, ls=(0, (3, 2)))
    X = R["xtab"]
    pct = lambda *ks: sum(X[k]["p"] for k in ks)
    labs = [ax.text(-1.75, 3.5, f"T abnormal, Z within range\n{pct('Low bone mass|within', 'Osteoporosis|within'):.1f}%",
                    ha="center", va="top", fontsize=8, color="#8a5a12"),
            ax.text(2.5, -2.75, f"Concordant normal\n{pct('Normal|within'):.1f}%", ha="right", va="top", fontsize=8, color=INK),
            ax.text(-0.85, -2.2, f"Z ≤ −2.0: {pct('Normal|below', 'Low bone mass|below', 'Osteoporosis|below'):.1f}%",
                    ha="left", va="top", fontsize=8, color=RED)]
    # 自檢：標籤框內不可有散點、不可被參考線穿過
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    for t in labs:
        (x0, y0), (x1, y1) = inv.transform(t.get_window_extent(r))
        inside = ((P["T"] > x0) & (P["T"] < x1) & (P.Z > y0) & (P.Z < y1)).sum()
        assert inside == 0, f"label over points: {t.get_text()!r} ({inside})"
        for xv in (-1, -2.5):
            assert not (x0 < xv < x1 and y0 < ytop), f"vertical line through {t.get_text()!r}"
        assert not (y0 < -2 < y1), f"horizontal line through {t.get_text()!r}"
    ax.set_xlabel("T-score (postmenopausal rule)")
    ax.set_ylabel("Z-score (premenopausal rule)")
    # B：依單一年齡的不一致比例
    ax = bare(fig.add_subplot(gs[1])); panel(ax, "B")
    A = R["age1"]; ages = np.array(sorted(int(a) for a in A))
    for k, c, lab in (("any", AMBER, "T abnormal, Z within range"), ("op", RED, "T osteoporosis, Z within range")):
        p = np.array([A[str(a)][k]["p"] for a in ages]); lo = np.array([A[str(a)][k]["lo"] for a in ages]); hi = np.array([A[str(a)][k]["hi"] for a in ages])
        ax.fill_between(ages, lo, hi, color=c, alpha=0.18, lw=0)
        ax.plot(ages, p, "-o", color=c, ms=3, lw=1.2, label=lab)
    ax.axvspan(44.5, 55.5, color=FAINT, alpha=0.35, lw=0, zorder=0)
    ax.text(50, 76, "Primary age range", ha="center", va="top", fontsize=8, color=INK)
    ax.set_xlim(39.5, 59.5); ax.set_ylim(0, 78)
    ax.set_xticks(range(40, 60, 5))
    ax.set_xlabel("Age, years"); ax.set_ylabel("Women with discordant category, %")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.70), frameon=False, handlelength=1.6, fontsize=8)
    save(fig, "Figure_1")


# ------------------------------------------------------------------ Figure 2
STAT = [("Menstruating", "Menstruating"), ("No periods, menopause", "No periods, natural menopause"),
        ("Hysterectomy, ovaries conserved", "Hysterectomy, ovaries conserved"), ("Bilateral oophorectomy", "Bilateral oophorectomy")]


def fig2():
    fig = plt.figure(figsize=(fs.WIDTH_IN, 3.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=1.0, left=0.31, right=0.93, bottom=0.22, top=0.80)
    ax = bare(fig.add_subplot(gs[0])); panel(ax, "A")
    S, SH = R["status"], R["status_share"]
    for i, (k, lab) in enumerate(STAT):
        y = len(STAT) - 1 - i
        ci_bar(ax, y, S[k]["any"], AMBER, dy=0.13, label="T abnormal, Z within range" if i == 0 else None)
        ci_bar(ax, y, S[k]["op"], RED, marker="s", dy=-0.13, ms=4.5, label="T osteoporosis, Z within range" if i == 0 else None)
    ax.set_yticks(range(len(STAT)))
    ax.set_yticklabels([f"{lab}\n(n = {S[k]['any']['n']:,}; {SH[k]['p']:.1f}%)" for k, lab in STAT][::-1], fontsize=8)
    ax.set_xlim(0, 60); ax.set_ylim(-0.6, len(STAT) - 0.25)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xlabel("Women with discordant category, %")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=8, ncol=1, handlelength=1.2,
              borderaxespad=0.1)
    # B：骨折史（年齡標準化）
    ax = bare(fig.add_subplot(gs[1])); panel(ax, "B")
    O, OS = R["outcome"]["frac"], R["outcome_std"]["frac"]
    G = [("T normal", "Concordant normal", GREY), ("T abnormal, Z within range", "T abnormal,\nZ within range", AMBER)]
    for i, (g, lab, c) in enumerate(G):
        y = len(G) - 1 - i
        ci_bar(ax, y, OS[g], c)
        ax.text(OS[g]["hi"] + 0.6, y, f"{OS[g]['p']:.1f}%", va="center", fontsize=8, color=INK)
    ax.set_yticks(range(len(G)))
    ax.set_yticklabels([f"{lab}\n(n = {O[g]['n']:,})" for g, lab, _ in G][::-1], fontsize=8)
    ax.set_title("Hip, wrist or spine", fontsize=8.5, loc="left", pad=4)
    ax.set_xlim(0, 24); ax.set_ylim(-0.6, len(G) - 0.4)
    ax.set_xticks([0, 5, 10, 15, 20])
    ax.set_xlabel("Prior fracture, %\n(age-standardized)")
    d = OS["T abnormal, Z within range"]
    ax.text(23.5, -0.45, f"Difference {d['d']:.1f} points\n(95% CI {d['dlo']:.1f} to {d['dhi']:.1f})",
            ha="right", va="bottom", fontsize=8, color=INK)
    save(fig, "Figure_2")


# ------------------------------------------------------------------ Figure S1：流程圖
def figS1():
    F = dict(R["flow"])
    R3 = json.load(open(f"{fs.ROOT}/out/results3.json")); FX = R3["flow_extra"]
    fig = plt.figure(figsize=(fs.WIDTH_IN, 5.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    main = [
        (93, f"Women in NHANES 2005–2010, 2013–2014 and 2017–2018\nn = {F['Women, 5 cycles']:,}"),
        (78, f"Aged 40–59 years\nn = {F['Aged 40–59 y']:,}"),
        (55, f"Valid femoral neck DXA, aged 40–59 years\nn = {F['Analysed, aged 40–59 y']:,}\nUsed in Figures 1B, 2A–2C and S5"),
        (33, f"Primary sample, aged 45–55 years\nn = {F['Primary sample, aged 45–55 y']:,}\nUsed in Tables 1–2, S2–S4; Figures 1C, 1D, 2D, 3, S2–S4"),
        (8, f"Valid lumbar spine DXA, aged 45–55 years\nn = {FX['ls_valid']:,}\nUsed for lumbar spine and lowest-of-three-sites results"),
    ]
    lm = FX["ls_missing"]
    side = [
        (66.5, "Excluded\n"
               f"Aged 40–49 y in 2017–2018\n  (femur DXA offered ≥50 y only): {F['Excluded: 2017–2018 femur DXA limited to ≥50 y']:,}\n"
               f"Not scanned, pregnant: {F['Excluded: not scanned, pregnant']:,}\n"
               f"Not scanned, weight over DXA table limit\n  (300 lb in 2005–2010, 450 lb in 2013–2018): {F['Excluded: not scanned, weight over DXA table limit']:,}\n"
               f"Not scanned, other reason: {F['Excluded: not scanned, other reason']:,}\n"
               f"Scan completed, data invalid: {F['Excluded: scan invalid']:,}\n"
               f"No femur scan record: {F['Excluded: no femur scan record']:,}"),
        (44, f"Not in primary age range\nAged 40–44 or 56–59 y: {FX['excluded_age']:,}"),
        (20.5, f"No valid lumbar spine result: {sum(lm.values()):,}\n"
               f"One or more vertebrae invalid: {lm['one or more vertebrae invalid']:,}\n"
               f"Not scanned, other reason: {lm['not scanned, other reason']:,}\nClassified from hip sites only"),
    ]
    box = dict(boxstyle="round,pad=0.5", fc="white", ec="#5b6b7a", lw=0.8)
    T = [ax.text(28, y, s_, ha="center", va="center", fontsize=8.5, bbox=box) for y, s_ in main]
    E = [ax.text(78, y, s_, ha="center", va="center", fontsize=8, bbox=box, linespacing=1.3) for y, s_ in side]
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    bb = lambda t: inv.transform(t.get_bbox_patch().get_window_extent(r))
    GAP = 1.0
    arrow = dict(arrowstyle="-|>", color="#5b6b7a", lw=0.8, shrinkA=0, shrinkB=0, mutation_scale=8)
    for a, b in zip(T[:-1], T[1:]):
        (_, y0a), _ = bb(a); _, (_, y1b) = bb(b)
        ax.annotate("", xy=(28, y1b + GAP), xytext=(28, y0a - GAP), arrowprops=arrow)
    for (y, _), e in zip(side, E):
        (ex0, _), _ = bb(e)
        ax.annotate("", xy=(ex0 - GAP, y), xytext=(28, y), arrowprops=arrow)
        assert ex0 - GAP > max(bb(t)[1][0] for t in T), "side arrow tip inside main box column"
        for t in T:                                      # 側向箭頭的 y 不可穿過主欄的框
            (_, ty0), (_, ty1) = bb(t)
            assert not (ty0 - GAP < y < ty1 + GAP), f"side arrow at y={y} crosses a main box"
    save(fig, "Figure_S1")


if __name__ == "__main__":
    figS1()          # 10-27 第二版：Figure 1–4 改由 figs2.py 產生，這裡只產流程圖
