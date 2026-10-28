"""第四版分析（2026-10-27，第三輪外部審閱）：全部為事後分析，登記在 Table S1。

1. Figure 3C 補列：Z 二分改用與 T ≤ −1.0 同陽性率的切點、Z ≤ −1.0，以及 T ≤ −1.0 二分（與 Z 二分同為兩級，比較才公平）。
2. 期望人數的敏感度：保留卵巢的子宮切除者卵巢功能可能較早衰退 → P(停經前 | 年齡) 取 a + 1、a + 2 歲的值。
3. 過渡晚期：只有 2005–2006 問「最後一次月經是幾個月前」（RHQ051）；仍有月經但 ≥ 3 個月前（STRAW 過渡晚期以 ≥ 60 天無月經定義，問卷最接近的分界）
   → ISCD 也允許兩種規則；估此群中不一致者的人數（2005–2006 母群）。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

import analysis as A
import analysis2 as B
import analysis3 as C
from analysis import Design, OUT, D


def main():
    w, flow = B.load_all()
    des = Design(w); q = stats.t.ppf(0.975, des.df)
    ages = w.RIDAGEYR.values.astype(float)
    a40 = (ages >= 40) & (ages <= 59); prim = (ages >= 45) & (ages <= 55)
    FN = B.site_scores(w, "FN"); w["Z"], w["tcat"] = FN; w["T"] = (w.DXXNKBMD - A.T_MEAN) / A.T_SD
    F = B.flips(*FN); w["disc"] = F["any"]
    indet = w.status.isin(B.INDET).values
    R = {}

    # ---- 1. C 統計量補列
    y = w.frac.values; Z, T = w.Z.values, w["T"].values
    m = prim & ~np.isnan(y)
    tpos = des.prop((w.tcat.values >= 1).astype(float), prim)["p"] / 100
    zz, ww = Z[prim], des.w[prim]; o = np.argsort(zz); cw = np.cumsum(ww[o]) / ww.sum()
    zcut = float(zz[o][np.searchsorted(cw, tpos)])
    R["z_matched_cut"] = zcut; R["t_positivity"] = 100 * tpos
    R["z_matched_positivity"] = des.prop((Z <= zcut).astype(float), prim)["p"]
    R["z_m1_positivity"] = des.prop((Z <= -1).astype(float), prim)["p"]
    scores = {"T-score, continuous": -T, "Z-score, continuous": -Z,
              "T-score category (3 levels)": w.tcat.values, "Z-score ≤ −2.0 (2 levels)": (Z <= -2).astype(float),
              "T-score ≤ −1.0 (2 levels)": (w.tcat.values >= 1).astype(float),
              "Z-score at matched positivity (2 levels)": (Z <= zcut).astype(float),
              "Z-score ≤ −1.0 (2 levels)": (Z <= -1).astype(float)}
    aucs = lambda wt: {k: C.wauc(s[m], y[m], wt[m]) for k, s in scores.items()}
    a0 = aucs(des.w); rng = np.random.default_rng(7)
    bt = pd.DataFrame([aucs(des.w * C.rao_wu(des, rng)) for _ in range(C.NBOOT)])
    R["auc"] = {k: {"c": v, "lo": float(np.percentile(bt[k], 2.5)), "hi": float(np.percentile(bt[k], 97.5))} for k, v in a0.items()}
    for a, b in (("T-score category (3 levels)", "Z-score ≤ −2.0 (2 levels)"), ("T-score, continuous", "Z-score, continuous"),
                 ("T-score ≤ −1.0 (2 levels)", "Z-score at matched positivity (2 levels)"), ("Z-score at matched positivity (2 levels)", "Z-score ≤ −2.0 (2 levels)")):
        d = bt[a] - bt[b]
        R["auc"][f"{a} minus {b}"] = {"c": a0[a] - a0[b], "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5))}
    R["corr_TZ_45_55"] = float(np.corrcoef(T[prim], Z[prim])[0, 1])

    # ---- 2. 機率平移
    intact = np.isin(w.status.values, ["Menstruating", "No periods, menopause"]) & a40
    yprem = (w.status.values == "Menstruating").astype(float)
    kp = np.percentile(ages[intact], [5, 35, 65, 95])
    dis = np.nan_to_num(w.disc.values)
    def counts(wt):
        b = C.fit_logit(B.rcs(ages[intact], kp), yprem[intact], wt[intact])
        out = {}
        for sh in (0, 1, 2):
            p = 1 / (1 + np.exp(-(np.column_stack([np.ones(len(ages)), B.rcs(ages + sh, kp)]) @ b)))
            out[f"shift{sh}"] = (wt * dis * prim * indet * p).sum() / 1e6
        return out
    est = counts(des.w); rng = np.random.default_rng(20261028)
    bt = pd.DataFrame([counts(des.w * C.rao_wu(des, rng)) for _ in range(C.NBOOT)])
    R["shift"] = {k: {"t": v, "lo": float(np.percentile(bt[k], 2.5)), "hi": float(np.percentile(bt[k], 97.5))} for k, v in est.items()}

    # ---- 3. 過渡晚期（2005–2006）
    r = pd.read_sas(D / "RHQ_D.xpt")[["SEQN", "RHQ051"]]
    w2 = w.merge(r, on="SEQN", how="left")
    assert (w2.SEQN.values == w.SEQN.values).all()
    dD = Design(w.assign(w=w.WTMEC2YR))              # 單一週期：用原始兩年權重
    inD = (w.cycle == "D").values
    men = inD & prim & (w.status.values == "Menstruating") & w2.RHQ051.isin([1, 2, 3, 4, 5]).values
    late = men & w2.RHQ051.isin([3, 4, 5]).values
    R["late"] = {"n_menstruating": int(men.sum()), "n_late": int(late.sum()),
                 "share_late": dD.prop(w2.RHQ051.isin([3, 4, 5]).astype(float).values, men),
                 "disc_late": dD.prop(w.disc.values, late),
                 "count_late_disc": C_total(dD, w.disc.values, late),
                 "count_indet_disc_D": C_total(dD, w.disc.values, inD & prim & indet),
                 "N_45_55_D": float((dD.w * (inD & prim)).sum() / 1e6)}
    # 2026-10-28 補進程式（先前是另外算好寫進 results4.json）：未加權人數、相對標準誤、自由度
    dis = np.nan_to_num(w.disc.values)
    for key, dom in (("late_disc", late), ("indet_disc_D", inD & prim & indet)):
        z = dD.w * dom * dis
        R["late"][f"n_{key}"] = int((dom & (dis == 1)).sum())
        R["late"][f"rse_{key}"] = float(100 * np.sqrt(dD._var(z)) / z.sum())
    R["late"]["df_D"] = int(dD.df)
    json.dump(R, open(OUT / "results4.json", "w"), indent=1, default=float, ensure_ascii=False)
    print("T positivity %.1f; matched Z cut %.2f (pos %.1f); Z≤−1 pos %.1f; corr %.3f" % (R["t_positivity"], zcut, R["z_matched_positivity"], R["z_m1_positivity"], R["corr_TZ_45_55"]))
    for k, v in R["auc"].items():
        print(f"  {k:70s} {v['c']:.3f} ({v['lo']:.3f}–{v['hi']:.3f})")
    print("shift", {k: (round(v["t"], 2), round(v["lo"], 2), round(v["hi"], 2)) for k, v in R["shift"].items()})
    L = R["late"]
    print("late: n men %d, late %d, share %.1f%% (%.1f–%.1f), disc %.1f%%, count %.2fM (%.2f–%.2f); indet disc in D %.2fM; N D %.1fM" % (
        L["n_menstruating"], L["n_late"], L["share_late"]["p"], L["share_late"]["lo"], L["share_late"]["hi"], L["disc_late"]["p"],
        L["count_late_disc"]["t"], L["count_late_disc"]["lo"], L["count_late_disc"]["hi"], L["count_indet_disc_D"]["t"], L["N_45_55_D"]))


def C_total(des, y, dom):
    """加權總數（百萬）＋ 線性化 CI，下限截 0。"""
    return B.total(des, y, dom)


if __name__ == "__main__":
    main()
