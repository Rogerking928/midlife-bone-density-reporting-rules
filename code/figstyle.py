"""本專案共用的繪圖設定。

沿用 Bcc 生合成島圖譜／墨西哥剖腹產那篇的版式：DejaVu Sans、細軸線、**無格線**，
以印出來的寬度作畫（180 mm = 7.09 in），存 700 dpi PNG。
**只存 PNG，不存 PDF**（使用者明確要求過兩次）。
"""
import os
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out", os.environ.get("FIG_OUT", ""))   # FIG_OUT：另存到 out/ 下的子資料夾

WIDTH_IN = 6.50          # 165 mm；圖內字最小 8 pt
DPI = int(os.environ.get("FIG_DPI", 700))   # FIG_DPI 可改解析度

NAVY = "#31485d"
RED = "#a4303f"
AMBER = "#b8791f"
GREEN = "#4f7d5b"
PURPLE = "#6b5b95"
TEAL = "#2f7d8c"
BROWN = "#7a5c3e"
PLUM = "#9c6b8a"
GREY = "#8a939c"
FAINT = "#c9d1d8"
INK = "#1f2933"

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.unicode_minus": True,      # 負號用 −，不用連字號
    "axes.edgecolor": "#5b6b7a",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": "#5b6b7a",
    "ytick.color": "#5b6b7a",
    "xtick.labelcolor": INK,      # 刻度字用深色（外部審閱：灰字比資料還淡）
    "ytick.labelcolor": INK,
    "figure.facecolor": "white",
    "hatch.linewidth": 0.6,
})


def bare(ax, keep=("bottom", "left")):
    for side, sp in ax.spines.items():
        sp.set_visible(side in keep)
    ax.grid(False)
    return ax


def panel(ax, letter, dx=None, dy=None, left=None):
    """分格字母：統一放在該軸（含刻度字與軸標題）左上角的外面，存檔前才定位。
    left=True：強制對齊左邊界（上下疊放、但軸本身往右縮的分格，例如 Figure 1B）。"""
    ax._panel_letter = letter
    ax._panel_left = left


def place_letters(fig):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    axs = [ax for ax in fig.axes if getattr(ax, "_panel_letter", None)]
    # 10-26 外部審閱：左欄（軸左緣在畫布左 40% 內）的字母全部對齊同一個左邊界；其他欄對齊自己的最左側
    left = [ax for ax in axs if (ax.get_position().x0 < 0.4 if getattr(ax, "_panel_left", None) is None else ax._panel_left)]
    x_left = min(ax.get_tightbbox(r).x0 for ax in left) if left else None
    for ax in axs:
        L = ax._panel_letter
        x0 = x_left if ax in left else ax.get_tightbbox(r).x0
        top = ax.get_window_extent(r).y1
        fig.text(x0 / W, (top + 8 * fig.dpi / 72) / H, L, fontsize=11, fontweight="bold",
                 ha="left", va="bottom")


def check(fig, stem):
    """畫完後的自動檢查：撐寬是無聲失敗，必須擋下來（字會整頁變小）。"""
    probs = []
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    if bb.width > WIDTH_IN + 0.02:
        probs.append(f"tight bbox {bb.width:.2f} in > {WIDTH_IN} in")
    r = fig.canvas.get_renderer()
    items = list(fig.texts)
    for ax in fig.axes:
        items += [t for t in ax.texts]
        if ax.axison:
            items += ax.get_xticklabels() + ax.get_yticklabels() + [ax.xaxis.label, ax.yaxis.label]
        items += [ax.title]
        if ax.get_legend(): items += ax.get_legend().get_texts()
    for lg in fig.legends: items += lg.get_texts()
    boxes = [(t, t.get_window_extent(r)) for t in items if t.get_visible() and t.get_text().strip()]
    pt = fig.dpi / 72
    for t, bb_ in boxes:   # 圖內字最小 8 pt
        if t.get_fontsize() < 7.99:
            probs.append(f"font < 8 pt: {t.get_text()!r} ({t.get_fontsize():.1f})")
    fams = {t.get_fontname() for t, _ in boxes}
    if fams - {"Arial"}:
        probs.append(f"non-Arial text: {fams}")
    for t, bb_ in boxes:   # 字不可超出畫布（存檔時會被切掉）
        if t not in fig.texts and (bb_.x0 < 0 or bb_.x1 > fig.bbox.width or bb_.y0 < 0 or bb_.y1 > fig.bbox.height):   # 10-27：存檔會切在畫布邊，容差改 0
            probs.append(f"outside canvas: {t.get_text()!r}")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            w = min(a.x1, b.x1) - max(a.x0, b.x0); h = min(a.y1, b.y1) - max(a.y0, b.y0)
            if w > 0.6 * pt and h > 0.6 * pt:
                probs.append(f"overlap: {boxes[i][0].get_text()!r} / {boxes[j][0].get_text()!r}")
    return probs


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    place_letters(fig)
    probs = check(fig, stem)
    if probs:
        raise SystemExit(f"{stem}: " + "; ".join(probs[:5]))
    # 10-26：把圖上印出來的所有文字（含刻度字）倒成文字檔，給 check_manuscript.py 當「圖上看得到的數字」
    os.makedirs(os.path.join(OUT, "fig_text"), exist_ok=True)
    fig.canvas.draw()
    txt = [t.get_text() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    open(os.path.join(OUT, "fig_text", f"{stem}.txt"), "w").write("\n".join(txt))
    png = os.path.join(OUT, f"{stem}.png")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    from PIL import Image          # 存成 RGB（無透明通道）
    im = Image.open(png); im.load()
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, "white"); bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        bg.save(png, dpi=(DPI, DPI))
    print(f"  fig -> {os.path.basename(png)}", flush=True)
