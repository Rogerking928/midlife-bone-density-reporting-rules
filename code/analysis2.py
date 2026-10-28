"""第二版分析（2026-10-27，回應外部審閱）：全部為事後分析，登記在 Table S1。

1. 腰椎 DXA：T 參考取本資料 2005–2010 年 20–29 歲非西裔白人女性（廠商參考值查不到公開出處，不用）；「三個部位取最差」＝ISCD 臨床做法。
2. 年齡的 restricted cubic spline（調查加權 logistic、設計為本的三明治變異數）。
3. 停經狀態「結構上無法判定」的比例與受影響人數（加權總數＋CI）。
4. 異質性（狀態、族裔、BMI、荷爾蒙、抽菸），附美國人數。
5. 骨折史對 T 分數的 spline（年齡邊際標準化），以及校正後盛行率比（modified Poisson）。
6. 規格曲線：部位 × T 參考 × Z 參考 × 加權 × 是否排除荷爾蒙使用者。
"""
import itertools, json
import numpy as np
import pandas as pd
from scipy import stats

import analysis as A
from analysis import Design, korn_graubard, zref, OUT, D

INDET = ["Hysterectomy, ovaries conserved", "No periods, other reason", "Not ascertained"]


# ------------------------------------------------------------------ 資料
def load_all():
    x, flow = A.load()
    extra = []
    for c in A.CYC:
        s = pd.read_sas(D / f"DXXSPN_{c}.xpt")[["SEQN", "DXASPNST", "DXXOSBMD"]]
        m = pd.read_sas(D / f"SMQ_{c}.xpt")[["SEQN", "SMQ020", "SMQ040"]]
        o = pd.read_sas(D / f"OSQ_{c}.xpt")[["SEQN", "OSQ130", "OSQ170", "OSQ200"]]
        extra.append(s.merge(m, on="SEQN", how="outer").merge(o, on="SEQN", how="outer"))
    x = x.merge(pd.concat(extra), on="SEQN", how="left")
    w = A.prepare(x)
    w["status"] = A.status(w)
    w["LS"] = np.where(w.DXASPNST == 1, w.DXXOSBMD, np.nan)
    w["smoke"] = np.where(w.SMQ040.isin([1, 2]), 1.0, np.where(w.SMQ040 == 3, 0.0, np.where(w.SMQ020 == 2, 0.0, np.nan)))
    w["par_hip"] = np.where(w[["OSQ170", "OSQ200"]].isin([1]).any(axis=1), 1.0, np.where(w[["OSQ170", "OSQ200"]].isin([2]).all(axis=1), 0.0, np.nan))
    w["frac"] = np.where(w[["OSQ010A", "OSQ010B", "OSQ010C"]].isin([1]).any(axis=1), 1.0,
                         np.where(w[["OSQ010A", "OSQ010B", "OSQ010C"]].isin([2]).all(axis=1), 0.0, np.nan))
    w["told"] = np.where(w.OSQ060 == 1, 1.0, np.where(w.OSQ060 == 2, 0.0, np.nan))
    w["hrt"] = np.where(w.RHQ540 == 1, 1.0, np.where(w.RHQ540 == 2, 0.0, np.nan))
    w["race4"] = np.select([w.RIDRETH1 == 3, w.RIDRETH1 == 4, w.RIDRETH1.isin([1, 2])], ["Non-Hispanic White", "Non-Hispanic Black", "Hispanic"], "Other, including multiracial")
    w["age10"] = (w.RIDAGEYR // 10 * 10).astype(int)
    return w.reset_index(drop=True), flow


def internal_t(w, col):
    y = w[(w.RIDRETH1 == 3) & w.RIDAGEYR.between(20, 29) & w[col].notna() & w.cycle.isin(["D", "E", "F"])]
    mu = np.average(y[col], weights=y.w); sd = np.sqrt(np.average((y[col] - mu) ** 2, weights=y.w))
    return mu, sd, len(y)


SITES = {"FN": ("DXXNKBMD", (A.FN_OP, A.FN_LOW)), "TH": ("DXXOFBMD", (A.TH_OP, A.TH_LOW)), "LS": ("LS", None)}


def site_scores(w, site, tref="published", zby=("RIDRETH1", "age5")):
    col, cuts = SITES[site]
    if tref == "internal" or cuts is None:
        mu, sd, _ = internal_t(w, col); cuts = (mu - 2.5 * sd, mu - sd)
    ref = zref(w, col, by=zby)
    j = w.join(ref, on=list(zby))
    Z = ((j[col] - j.m) / j.sd).values
    v = w[col].values
    tcat = np.where(np.isnan(v), np.nan, np.where(v <= cuts[0], 2, np.where(v <= cuts[1], 1, 0)))
    return Z, tcat


def combine(parts):
    """多部位取最差：T 分類取最高、Z 取最低（有資料的部位）。"""
    Zs = np.vstack([p[0] for p in parts]); Ts = np.vstack([p[1] for p in parts])
    with np.errstate(all="ignore"):
        return np.nanmin(Zs, axis=0), np.nanmax(Ts, axis=0)


def flips(Z, tcat):
    ok = ~np.isnan(Z) & ~np.isnan(tcat)
    within = Z > -2.0
    f = lambda m: np.where(ok, m.astype(float), np.nan)
    return {"any": f((tcat >= 1) & within), "op": f((tcat == 2) & within)}


# ------------------------------------------------------------------ 調查設計：總數、GLM
def total(des, y, dom):
    dom = np.asarray(dom, bool) & ~np.isnan(np.asarray(y, float))
    z = des.w * dom * np.nan_to_num(np.asarray(y, float))
    t = z.sum(); se = np.sqrt(des._var(z)); q = stats.t.ppf(0.975, des.df)
    return {"t": t / 1e6, "lo": max(0.0, (t - q * se) / 1e6), "hi": (t + q * se) / 1e6}


def rcs(x, knots):
    """Harrell restricted cubic spline 基底（不含截距與線性項以外的常數）。"""
    k = np.asarray(knots, float); K = len(k)
    p = lambda u: np.clip(u, 0, None) ** 3
    cols = [x]
    for j in range(K - 2):
        cols.append((p(x - k[j]) - p(x - k[K - 2]) * (k[K - 1] - k[j]) / (k[K - 1] - k[K - 2])
                     + p(x - k[K - 1]) * (k[K - 2] - k[j]) / (k[K - 1] - k[K - 2])) / (k[K - 1] - k[0]) ** 2)
    return np.column_stack(cols)


def svyglm(des, X, y, dom, family="logit"):
    """加權 GLM（logit 或 log-link Poisson），變異數用 strata × PSU 的線性化三明治；域外的列分數為 0，保留完整設計。"""
    y = np.asarray(y, float); dom = np.asarray(dom, bool) & ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    X1 = np.column_stack([np.ones(len(y)), X]); Xd = np.where(dom[:, None], X1, 0.0); yd = np.where(dom, y, 0.0)
    w = des.w * dom
    b = np.zeros(X1.shape[1])
    if family == "poisson":
        b[0] = np.log(max((w * yd).sum() / w.sum(), 1e-6))
    for _ in range(100):
        eta = Xd @ b
        mu = 1 / (1 + np.exp(-eta)) if family == "logit" else np.exp(eta)
        v = mu * (1 - mu) if family == "logit" else mu
        H = (Xd * (w * v)[:, None]).T @ Xd
        U = Xd.T @ (w * (yd - mu))
        step = np.linalg.solve(H, U); b += step
        if np.abs(step).max() < 1e-10:
            break
    eta = Xd @ b; mu = 1 / (1 + np.exp(-eta)) if family == "logit" else np.exp(eta)
    Ui = Xd * (w * (yd - mu))[:, None]
    Hinv = np.linalg.inv(H)
    t = pd.DataFrame(Ui); t["h"] = des.h; t["p"] = des.psu
    S = np.zeros((X1.shape[1],) * 2)
    for _, g in t.groupby("h"):
        P = g.groupby("p").sum(numeric_only=True).drop(columns=["h"], errors="ignore").values
        n = len(P)
        if n > 1:
            dP = P - P.mean(axis=0); S += n / (n - 1) * dP.T @ dP
    return b, Hinv @ S @ Hinv, dom


def logit_curve(des, x, y, dom, grid, knots):
    b, V, _ = svyglm(des, rcs(x, knots), y, dom)
    G = np.column_stack([np.ones(len(grid)), rcs(grid, knots)])
    eta = G @ b; se = np.sqrt(np.einsum("ij,jk,ik->i", G, V, G)); q = stats.t.ppf(0.975, des.df)
    ex = lambda e: 100 / (1 + np.exp(-e))
    return {"x": list(map(float, grid)), "p": list(ex(eta)), "lo": list(ex(eta - q * se)), "hi": list(ex(eta + q * se))}


# ------------------------------------------------------------------ 主程式
def main():
    w, flow = load_all()
    des = Design(w); q = stats.t.ppf(0.975, des.df)
    a40 = w.RIDAGEYR.between(40, 59).values; prim = w.RIDAGEYR.between(45, 55).values
    R = {"df": des.df}

    # 部位分數
    S = {s: site_scores(w, s) for s in SITES}
    FN, LOW3 = S["FN"], combine([S["FN"], S["TH"], S["LS"]])
    w["Z"], w["tcat"] = FN
    w["T"] = (w.DXXNKBMD - A.T_MEAN) / A.T_SD
    F = flips(*FN); w["disc"] = F["any"]; w["disc_op"] = F["op"]
    F3 = flips(*LOW3)
    R["tref_internal"] = {s: dict(zip(("mean", "sd", "n"), internal_t(w, SITES[s][0]))) for s in SITES}
    R["n_ls_45_55"] = int((prim & w.LS.notna().values).sum())

    # Table 2 精簡版：% 與百萬人（含 CI），股骨頸 vs 三部位取最差
    X2 = {}
    for name, (Z, tc) in (("FN", FN), ("LOW3", LOW3)):
        ok = ~np.isnan(Z) & ~np.isnan(tc)
        X2[name] = {}
        for c, lab in ((0, "Normal"), (1, "Low bone mass"), (2, "Osteoporosis")):
            for zc, zm in (("within", Z > -2.0), ("below", Z <= -2.0)):
                y = np.where(ok, ((tc == c) & zm).astype(float), np.nan)
                X2[name][f"{lab}|{zc}"] = {**des.prop(y, prim), "mil": total(des, y, prim)}
        fl = flips(Z, tc)
        X2[name]["primary"] = {**des.prop(fl["any"], prim), "mil": total(des, fl["any"], prim)}
        X2[name]["op"] = {**des.prop(fl["op"], prim), "mil": total(des, fl["op"], prim)}
        X2[name]["N"] = total(des, np.where(ok, 1.0, np.nan), prim)
    R["xtab2"] = X2

    # 年齡 spline（40–59）
    ages = w.RIDAGEYR.values.astype(float); grid = np.linspace(40, 59, 77)
    kn = np.percentile(ages[a40], [5, 35, 65, 95])
    R["spline_age"] = {"disc": logit_curve(des, ages, w.disc.values, a40, grid, kn),
                       "disc_op": logit_curve(des, ages, w.disc_op.values, a40 & (ages >= 45), grid[grid >= 45], np.percentile(ages[a40 & (ages >= 45)], [10, 50, 90])),
                       "indet": logit_curve(des, ages, w.status.isin(INDET).astype(float).values, a40, grid, kn),
                       "menstruating": logit_curve(des, ages, (w.status == "Menstruating").astype(float).values, a40, grid, kn),
                       "knots": list(kn)}
    R["age1_disc_op_note"] = "disc_op spline fitted from 45 y because no osteoporosis-label discordance occurred before 50 y"

    # 停經狀態結構上無法判定 × 不一致 → 人數
    bands = [(40, 44), (45, 49), (50, 54), (55, 59)]
    R["affected"] = {f"{lo}–{hi}": {s: total(des, w.disc.values, (ages >= lo) & (ages <= hi) & (w.status == s).values) for s in INDET}
                     for lo, hi in bands}
    for lab, m in (("40–59", a40), ("45–55", prim)):
        R["affected"][f"total {lab}"] = total(des, w.disc.values, m & w.status.isin(INDET).values)
        R["affected"][f"total {lab} op"] = total(des, w.disc_op.values, m & w.status.isin(INDET).values)
        R["affected"][f"indet share {lab}"] = des.prop(w.status.isin(INDET).astype(float).values, m)
        R["affected"][f"disc among indet {lab}"] = des.prop(w.disc.values, m & w.status.isin(INDET).values)

    # 異質性（45–55）
    H = {}
    def grp(section, lab, m):
        H.setdefault(section, {})[lab] = {"any": des.prop(w.disc.values, prim & m), "op": des.prop(w.disc_op.values, prim & m),
                                          "mil": total(des, w.disc.values, prim & m), "share": des.prop(m.astype(float), prim)}
    for s in ["Menstruating", "No periods, menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy", "No periods, other reason", "Not ascertained"]:
        grp("Menstrual and surgical status", s, (w.status == s).values)
    for r in ["Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other, including multiracial"]:
        grp("Race and ethnicity", r, (w.race4 == r).values)
    bmi = w.BMXBMI.values
    for lab, m in (("<25", bmi < 25), ("25 to <30", (bmi >= 25) & (bmi < 30)), ("≥30", bmi >= 30)):
        grp("Body mass index, kg/m²", lab, m)
    grp("Female hormone use", "Ever", (w.hrt == 1).values); grp("Female hormone use", "Never", (w.hrt == 0).values)
    grp("Current smoking", "Yes", (w.smoke == 1).values); grp("Current smoking", "No", (w.smoke == 0).values)
    R["hetero"] = H

    # Table 1 新增列
    T1 = {}
    for s in ["All", "Menstruating", "No periods, menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy"]:
        m = prim if s == "All" else prim & (w.status == s).values
        T1[s] = {"smoke": des.prop(w.smoke.values, m), "frac": des.prop(w.frac.values, m), "par_hip": des.prop(w.par_hip.values, m),
                 "race4": {r: des.prop((w.race4 == r).astype(float).values, m) for r in ["Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other, including multiracial"]},
                 "ls": des.mean(w.LS.values, m), "n_ls": int((m & w.LS.notna().values).sum())}
    R["table1_extra"] = T1

    # 骨折：spline 對 T（年齡邊際標準化）＋ 盛行率比
    T = w["T"].values; Tg = np.linspace(*np.percentile(T[prim], [1, 99]), 81)  # 第五輪審閱：與圖說的 1st–99th 一致
    kT = np.percentile(T[prim], [5, 35, 65, 95])
    X = np.column_stack([rcs(T, kT), ages])
    b, V, dom = svyglm(des, X, w.frac.values, prim)
    pred, lo, hi = [], [], []
    ww = des.w * dom
    for t in Tg:
        Xt = np.column_stack([np.ones(len(T)), rcs(np.full(len(T), t), kT), ages])
        mu = 1 / (1 + np.exp(-(Xt @ b)))
        p = (ww * mu).sum() / ww.sum()
        g = (Xt * (ww * mu * (1 - mu))[:, None]).sum(axis=0) / ww.sum()
        se = np.sqrt(g @ V @ g); pred.append(100 * p); lo.append(100 * (p - q * se)); hi.append(100 * (p + q * se))
    R["frac_curve"] = {"T": list(Tg), "p": pred, "lo": lo, "hi": hi, "knots": list(kT),
                       "T_range_obs": [float(np.percentile(T[prim], 1)), float(np.percentile(T[prim], 99))]}
    # 盛行率比：不一致 vs 一致正常（排除 Z 低於同齡範圍的 18 人）
    expo = np.where((w.tcat == 0) & (w.Z > -2), 0.0, np.where((w.tcat >= 1) & (w.Z > -2), 1.0, np.nan))
    race_d = pd.get_dummies(w.race4).drop(columns="Non-Hispanic White").astype(float).values
    models = {"Crude": np.column_stack([expo]),
              "Age-adjusted": np.column_stack([expo, ages]),
              "Adjusted for age, BMI, race and ethnicity, hormone use": np.column_stack([expo, ages, bmi, race_d, w.hrt.values])}
    PR = {}
    for o in ("frac", "told"):
        PR[o] = {}
        for k, Xm in models.items():
            b, V, dm = svyglm(des, Xm, w[o].values, prim & ~np.isnan(expo), family="poisson")
            PR[o][k] = {"pr": float(np.exp(b[1])), "lo": float(np.exp(b[1] - q * np.sqrt(V[1, 1]))), "hi": float(np.exp(b[1] + q * np.sqrt(V[1, 1]))),
                        "n": int(dm.sum())}
    R["pr"] = PR

    # 規格曲線
    specs = []
    sites = {"Femoral neck": ["FN"], "Total hip": ["TH"], "Lumbar spine": ["LS"], "Lower of femoral neck and total hip": ["FN", "TH"],
             "Lowest of femoral neck, total hip and lumbar spine": ["FN", "TH", "LS"]}
    zrefs = {"Race-specific, 5-year age bands": ("RIDRETH1", "age5"), "All races, 5-year age bands": ("age5",), "Race-specific, 10-year age bands": ("RIDRETH1", "age10")}
    cache = {}
    du = Design(w.assign(w=1.0))
    seen = set()
    for (sn, ss), tr, (zn, zb), wt, hx in itertools.product(sites.items(), ("published", "internal"), zrefs.items(), ("Weighted", "Unweighted"), (False, True)):
        if ss == ["LS"] and tr == "internal":
            continue          # 腰椎只有本資料參考，T 參考這一維不適用
        parts = []
        for s in ss:
            key = (s, tr, zb)
            if key not in cache:
                cache[key] = site_scores(w, s, tref=tr, zby=zb)
            parts.append(cache[key])
        Z, tc = parts[0] if len(parts) == 1 else combine(parts)
        fl = flips(Z, tc)
        dm = prim & ((w.hrt == 0).values if hx else True)
        d = des if wt == "Weighted" else du
        a, o = d.prop(fl["any"], dm), d.prop(fl["op"], dm)
        specs.append({"site": sn, "tref": "NHANES III cut-points" if tr == "published" and ss != ["LS"] else ("2005–2010 reference" if tr == "internal" else "2005–2010 reference (spine)"),
                      "zref": zn, "weighting": wt, "hormone": "Never-users only" if hx else "All women",
                      "n": a["n"], "any": a["p"], "any_lo": a["lo"], "any_hi": a["hi"], "op": o["p"], "op_lo": o["lo"], "op_hi": o["hi"],
                      "primary": sn == "Femoral neck" and tr == "published" and zn.startswith("Race-specific, 5") and wt == "Weighted" and not hx})
    R["specs"] = specs

    # 補充圖：散佈圖依參考組上色、逐歲點
    w.loc[prim, ["SEQN", "RIDAGEYR", "RIDRETH1", "race4", "status", "DXXNKBMD", "T", "Z", "tcat", "w", "frac"]].to_csv(OUT / "women45_55_v2.csv", index=False)
    w.loc[a40, ["SEQN", "RIDAGEYR", "RIDRETH1", "DXXNKBMD", "w", "disc"]].to_csv(OUT / "women40_59_v2.csv", index=False)
    ref = zref(w, "DXXNKBMD")
    R["zcut"] = {f"{int(r)}|{int(a)}": float(v.m - 2 * v.sd) for (r, a), v in ref.iterrows()}
    R["zmean"] = {f"{int(r)}|{int(a)}": float(v.m) for (r, a), v in ref.iterrows()}
    R["zsd"] = {f"{int(r)}|{int(a)}": float(v.sd) for (r, a), v in ref.iterrows()}
    json.dump(R, open(OUT / "results2.json", "w"), indent=1, default=float, ensure_ascii=False)

    # 摘要
    x = X2["FN"]; y3 = X2["LOW3"]
    print("FN primary %.1f%% (%.1f–%.1f) = %.2f M (%.2f–%.2f)" % (x["primary"]["p"], x["primary"]["lo"], x["primary"]["hi"], x["primary"]["mil"]["t"], x["primary"]["mil"]["lo"], x["primary"]["mil"]["hi"]))
    print("LOW3 primary %.1f%% (%.1f–%.1f) = %.2f M; op %.2f%%; n=%d" % (y3["primary"]["p"], y3["primary"]["lo"], y3["primary"]["hi"], y3["primary"]["mil"]["t"], y3["op"]["p"], y3["primary"]["n"]))
    print("n with LS 45–55:", R["n_ls_45_55"], "internal T refs", R["tref_internal"])
    af = R["affected"]
    print("indet share 45–55 %.1f%%; disc among indet %.1f%%; affected 45–55 %.2f M (%.2f–%.2f); 40–59 %.2f M; op 45–55 %.3f M" % (
        af["indet share 45–55"]["p"], af["disc among indet 45–55"]["p"], af["total 45–55"]["t"], af["total 45–55"]["lo"], af["total 45–55"]["hi"], af["total 40–59"]["t"], af["total 45–55 op"]["t"]))
    for k, v in PR.items():
        print(k, {m: (round(r["pr"], 2), round(r["lo"], 2), round(r["hi"], 2), r["n"]) for m, r in v.items()})
    sp = pd.DataFrame(specs)
    print("specs", len(sp), "any range %.1f–%.1f median %.1f; op range %.2f–%.2f" % (sp["any"].min(), sp["any"].max(), sp["any"].median(), sp.op.min(), sp.op.max()))
    print(sp.groupby("site")["any"].agg(["min", "median", "max"]).round(1))


if __name__ == "__main__":
    main()
