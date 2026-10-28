"""第三版分析（2026-10-27，第二輪外部審閱）：全部為事後分析，登記在 Table S1。

1. 期望誤分類人數：不一致 × 停經狀態結構上無法判定 × P(仍在停經前 | 年齡)；P 取子宮與卵巢完整、狀態已知的女性中仍有月經的比例（年齡 spline）。
   2 歲年齡帶；上限（不乘 P）與期望值都用 Rao–Wu PSU bootstrap（每次重配 spline）求 CI。
2. 百萬人數一律 = 該格占母群的比例（Korn–Graubard CI）× 加權總人數，零格的 CI 不再退化。
3. T 與 Z 正面比較（曾經骨折）：年齡校正下每低 1 單位的盛行率比；T 分類、Z 分類與連續 T、Z 的加權 C 統計量（bootstrap CI 與差值）。
4. 族裔一律用 Z 參考的五組（RIDRETH1）。
5. 規格曲線的變異分解（各設定的 type II 平方和占比）。
6. 流程圖補：年齡排除、腰椎缺失原因。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

import analysis as A
import analysis2 as B
from analysis import Design, OUT

RACE5 = {3: "Non-Hispanic White", 4: "Non-Hispanic Black", 1: "Mexican American", 2: "Other Hispanic", 5: "Other, including multiracial"}
INDET = B.INDET
BANDS2 = [(a, a + 1) for a in range(40, 60, 2)]
NBOOT = 1000


def rao_wu(des, rng):
    """每層抽 n_h − 1 個 PSU（有放回），權重 × n_h/(n_h − 1) × 次數。"""
    t = pd.DataFrame({"h": des.h, "p": des.psu})
    mult = {}
    for h, g in t.drop_duplicates().groupby("h"):
        ps = g.p.values; n = len(ps)
        draw = rng.choice(ps, size=n - 1, replace=True)
        for p in ps:
            mult[p] = (draw == p).sum() * n / (n - 1)
    return np.array([mult[p] for p in des.psu])


def fit_logit(X, y, w):
    X1 = np.column_stack([np.ones(len(y)), X]); b = np.zeros(X1.shape[1])
    for _ in range(60):
        mu = 1 / (1 + np.exp(-(X1 @ b))); H = (X1 * (w * mu * (1 - mu))[:, None]).T @ X1
        step = np.linalg.solve(H, X1.T @ (w * (y - mu))); b += step
        if np.abs(step).max() < 1e-9:
            break
    return b


def wauc(score, y, w):
    """加權 C 統計量：分數越高、越可能是病例；同分算 0.5。"""
    m = ~np.isnan(score) & ~np.isnan(y)
    s, y, w = score[m], y[m], w[m]
    df = pd.DataFrame({"s": s, "c": w * (y == 1), "k": w * (y == 0)}).groupby("s")[["c", "k"]].sum().sort_index()
    cum_k = df.k.cumsum().shift(fill_value=0)
    num = (df.c * (cum_k + 0.5 * df.k)).sum()
    return num / (df.c.sum() * df.k.sum())


def main():
    w, flow = B.load_all()
    des = Design(w); q = stats.t.ppf(0.975, des.df); rng = np.random.default_rng(20261027)
    ages = w.RIDAGEYR.values.astype(float)
    a40 = (ages >= 40) & (ages <= 59); prim = (ages >= 45) & (ages <= 55)
    S = {s: B.site_scores(w, s) for s in B.SITES}
    FN, LOW3 = S["FN"], B.combine([S["FN"], S["TH"], S["LS"]])
    w["Z"], w["tcat"] = FN; w["T"] = (w.DXXNKBMD - A.T_MEAN) / A.T_SD
    F = B.flips(*FN); w["disc"], w["disc_op"] = F["any"], F["op"]
    indet = w.status.isin(INDET).values
    R = {"df": des.df}

    def mil(y, dom_num, dom_den):
        """y 在 dom_num 為 1 的比例 × dom_den 的加權總人數（CI 用 Korn–Graubard）。"""
        yy = np.where(dom_num, np.asarray(y, float), 0.0)
        p = des.prop(yy, dom_den); N = (des.w * dom_den).sum() / 1e6
        return {"t": p["p"] / 100 * N, "lo": p["lo"] / 100 * N, "hi": p["hi"] / 100 * N, "pct": p}

    # ---- Table 2（百萬人改 KG）
    X2 = {}
    for name, (Z, tc) in (("FN", FN), ("LOW3", LOW3)):
        ok = ~np.isnan(Z) & ~np.isnan(tc); X2[name] = {}
        for c, lab in ((0, "Normal"), (1, "Low bone mass"), (2, "Osteoporosis")):
            for zc, zm in (("within", Z > -2.0), ("below", Z <= -2.0)):
                y = np.where(ok, ((tc == c) & zm).astype(float), np.nan)
                X2[name][f"{lab}|{zc}"] = {**des.prop(y, prim), "mil": mil(y, prim, prim)}
        fl = B.flips(Z, tc)
        X2[name]["primary"] = {**des.prop(fl["any"], prim), "mil": mil(fl["any"], prim, prim)}
        X2[name]["op"] = {**des.prop(fl["op"], prim), "mil": mil(fl["op"], prim, prim)}
        X2[name]["N"] = (des.w * prim).sum() / 1e6
    R["xtab"] = X2

    # ---- 停經狀態組成（年齡 spline、各類別分別配適後正規化）與不一致 spline
    grid = np.linspace(40, 59, 77); kn = np.percentile(ages[a40], [5, 35, 65, 95])
    comp_def = {"Menstruating": w.status.values == "Menstruating", "Natural menopause": w.status.values == "No periods, menopause",
                "Bilateral oophorectomy": w.status.values == "Bilateral oophorectomy", "Structurally indeterminate": indet}
    G = np.column_stack([np.ones(len(grid)), B.rcs(grid, kn)])
    comp = {}
    for k, m in comp_def.items():
        b, V, _ = B.svyglm(des, B.rcs(ages, kn), m.astype(float), a40)
        comp[k] = 1 / (1 + np.exp(-(G @ b)))
    tot = sum(comp.values())
    R["composition"] = {"x": list(grid), **{k: list(100 * v / tot) for k, v in comp.items()}}
    R["spline_disc"] = B.logit_curve(des, ages, w.disc.values, a40, grid, kn)

    # ---- P(停經前 | 年齡)：子宮與卵巢完整、狀態已知
    intact = np.isin(w.status.values, ["Menstruating", "No periods, menopause"]) & a40
    yprem = (w.status.values == "Menstruating").astype(float)
    kp = np.percentile(ages[intact], [5, 35, 65, 95])
    R["p_premenopausal"] = B.logit_curve(des, ages, yprem, intact, grid, kp)

    def counts(wt):
        b = fit_logit(B.rcs(ages[intact], kp), yprem[intact], wt[intact])
        p = 1 / (1 + np.exp(-(np.column_stack([np.ones(len(ages)), B.rcs(ages, kp)]) @ b)))
        dis = np.nan_to_num(w.disc.values)
        out = {}
        for lo, hi in BANDS2:
            m = indet & (ages >= lo) & (ages <= hi)
            out[f"{lo}–{hi}|upper"] = (wt * dis * m).sum() / 1e6
            out[f"{lo}–{hi}|expected"] = (wt * dis * m * p).sum() / 1e6
        for lab, m in (("45–55", prim), ("40–59", a40)):
            out[f"{lab}|upper"] = (wt * dis * m * indet).sum() / 1e6
            out[f"{lab}|expected"] = (wt * dis * m * indet * p).sum() / 1e6
            out[f"{lab}|op_expected"] = (wt * np.nan_to_num(w.disc_op.values) * m * indet * p).sum() / 1e6
        return out
    est = counts(des.w)
    boot = pd.DataFrame([counts(des.w * rao_wu(des, rng)) for _ in range(NBOOT)])
    R["counts"] = {k: {"t": v, "lo": float(np.percentile(boot[k], 2.5)), "hi": float(np.percentile(boot[k], 97.5))} for k, v in est.items()}

    # ---- 狀態 forest（同 2C）與異質性（族裔五組）
    H = {}
    def grp(sec, lab, m):
        H.setdefault(sec, {})[lab] = {"any": des.prop(w.disc.values, prim & m), "op": des.prop(w.disc_op.values, prim & m),
                                      "mil": mil(w.disc.values, prim & m, prim), "share": des.prop(m.astype(float), prim)}
    for s in ["Menstruating", "No periods, menopause", "Bilateral oophorectomy", "Hysterectomy, ovaries conserved", "No periods, other reason", "Not ascertained"]:
        grp("Menstrual and surgical status", s, w.status.values == s)
    for code, lab in RACE5.items():
        grp("Race and ethnicity", lab, w.RIDRETH1.values == code)
    bmi = w.BMXBMI.values
    for lab, m in (("<25", bmi < 25), ("25 to <30", (bmi >= 25) & (bmi < 30)), ("≥30", bmi >= 30)):
        grp("Body mass index, kg/m²", lab, m)
    grp("Female hormone use", "Ever", w.hrt.values == 1); grp("Female hormone use", "Never", w.hrt.values == 0)
    grp("Current smoking", "Yes", w.smoke.values == 1); grp("Current smoking", "No", w.smoke.values == 0)
    R["hetero"] = H
    R["race5_table1"] = {s: {lab: des.prop((w.RIDRETH1.values == code).astype(float), prim & ((w.status.values == s) if s != "All" else True))
                             for code, lab in RACE5.items()} for s in ["All", "Menstruating", "No periods, menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy"]}

    # ---- T 與 Z 正面比較（曾經骨折，45–55）
    y = w.frac.values
    T, Z = w["T"].values, w.Z.values
    PR = {}
    for lab, X in (("T-score, per 1 unit lower", np.column_stack([-T, ages])), ("Z-score, per 1 unit lower", np.column_stack([-Z, ages])),
                   ("T-score category: low bone mass vs normal", None), ("Z-score: below vs within expected range", None)):
        if X is None:
            if lab.startswith("T"):
                d = np.where(w.tcat.values == 0, 0.0, np.where(w.tcat.values == 1, 1.0, np.nan))
            else:
                d = np.where(Z > -2, 0.0, 1.0)
            X = np.column_stack([d, ages])
        b, V, dm = B.svyglm(des, X, y, prim, family="poisson")
        PR[lab] = {"pr": float(np.exp(b[1])), "lo": float(np.exp(b[1] - q * np.sqrt(V[1, 1]))), "hi": float(np.exp(b[1] + q * np.sqrt(V[1, 1]))), "n": int(dm.sum())}
    b, V, dm = B.svyglm(des, np.column_stack([np.where(w.tcat.values == 2, 1.0, np.where(w.tcat.values == 0, 0.0, np.nan)), ages]), y, prim, family="poisson")
    PR["T-score category: osteoporosis vs normal"] = {"pr": float(np.exp(b[1])), "lo": float(np.exp(b[1] - q * np.sqrt(V[1, 1]))), "hi": float(np.exp(b[1] + q * np.sqrt(V[1, 1]))), "n": int(dm.sum())}
    R["pr_tz"] = PR
    scores = {"T-score category (3 levels)": w.tcat.values, "Z-score category (2 levels)": (Z <= -2).astype(float),
              "T-score, continuous": -T, "Z-score, continuous": -Z}
    m = prim & ~np.isnan(y)
    def aucs(wt):
        return {k: wauc(s[m], y[m], wt[m]) for k, s in scores.items()}
    a0 = aucs(des.w)
    rng2 = np.random.default_rng(7)
    bt = pd.DataFrame([aucs(des.w * rao_wu(des, rng2)) for _ in range(NBOOT)])
    R["auc"] = {k: {"c": v, "lo": float(np.percentile(bt[k], 2.5)), "hi": float(np.percentile(bt[k], 97.5))} for k, v in a0.items()}
    for a, b_ in (("T-score category (3 levels)", "Z-score category (2 levels)"), ("T-score, continuous", "Z-score, continuous")):
        d = bt[a] - bt[b_]
        R["auc"][f"{a} minus {b_}"] = {"c": a0[a] - a0[b_], "lo": float(np.percentile(d, 2.5)), "hi": float(np.percentile(d, 97.5))}
    R["n_frac"] = int(m.sum()); R["n_frac_cases"] = int(np.nansum(y[m]))

    # ---- 規格曲線變異分解
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    sp = pd.DataFrame(json.load(open(OUT / "results2.json"))["specs"])
    fit = smf.ols("any ~ C(site) + C(tref) + C(zref) + C(weighting) + C(hormone)", data=sp).fit()
    an = anova_lm(fit, typ=2)
    R["spec_anova"] = {k.replace("C(", "").replace(")", ""): float(v) for k, v in (an.sum_sq / an.sum_sq.sum()).items()}
    R["spec_primary_pct_rank"] = float((sp["any"] < sp.loc[sp.primary, "any"].iloc[0]).mean() * 100)

    # ---- 流程圖補
    an40 = dict(flow)["Analysed, aged 40–59 y"]
    R["flow_extra"] = {"excluded_age": int(an40 - prim.sum()),
                       "ls_valid": int((prim & w.LS.notna().values).sum()),
                       "ls_missing": {lab: int((prim & (w.DXASPNST.values == c)).sum()) for c, lab in
                                      ((2, "one or more vertebrae invalid"), (3, "not scanned, pregnant"), (4, "not scanned, weight over DXA table limit"), (5, "not scanned, other reason"))}}
    R["flow_extra"]["ls_missing"]["no spine scan record"] = int((prim & np.isnan(w.DXASPNST.values)).sum())
    # 腰椎缺失 vs 有：年齡、BMI（是否隨機）
    lsm = prim & w.LS.isna().values; lso = prim & w.LS.notna().values
    R["ls_missing_compare"] = {"age_missing": des.mean(ages, lsm), "age_valid": des.mean(ages, lso), "bmi_missing": des.mean(bmi, lsm), "bmi_valid": des.mean(bmi, lso),
                               "disc_missing": des.prop(w.disc.values, lsm), "disc_valid": des.prop(w.disc.values, lso)}

    w.loc[prim, ["SEQN", "RIDAGEYR", "RIDRETH1", "status", "DXXNKBMD", "T", "Z", "tcat", "w", "frac"]].to_csv(OUT / "women45_55_v3.csv", index=False)
    json.dump(R, open(OUT / "results3.json", "w"), indent=1, default=float, ensure_ascii=False)

    c = R["counts"]
    print("45–55 upper %.2f (%.2f–%.2f); expected %.2f (%.2f–%.2f); op expected %.3f" % (c["45–55|upper"]["t"], c["45–55|upper"]["lo"], c["45–55|upper"]["hi"],
          c["45–55|expected"]["t"], c["45–55|expected"]["lo"], c["45–55|expected"]["hi"], c["45–55|op_expected"]["t"]))
    print("40–59 upper %.2f expected %.2f" % (c["40–59|upper"]["t"], c["40–59|expected"]["t"]))
    print({b: round(c[f"{b[0]}–{b[1]}|expected"]["t"], 3) for b in BANDS2})
    print("PR", {k: (round(v["pr"], 2), round(v["lo"], 2), round(v["hi"], 2)) for k, v in PR.items()})
    print("AUC", {k: (round(v["c"], 3), round(v["lo"], 3), round(v["hi"], 3)) for k, v in R["auc"].items()})
    print("anova", {k: round(v, 3) for k, v in R["spec_anova"].items()}, "rank", R["spec_primary_pct_rank"])
    print("flow", R["flow_extra"]); print("ls missing", {k: (round(v.get("m", v.get("p")), 2)) for k, v in R["ls_missing_compare"].items()})
    print("normal|below mil", X2["FN"]["Normal|below"]["mil"])


if __name__ == "__main__":
    main()
