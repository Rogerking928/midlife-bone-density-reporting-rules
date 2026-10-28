"""第五版分析（2026-10-28，正文外部審閱的選擇性建議）：全部為事後分析，登記在 Table S1。

1. Table 1 補兩欄：無月經（其他原因）、狀態未確定。各列定義與 analysis.py／analysis2.py／analysis3.py 相同；
   先用同一函式重算「Menstruating」欄，與既有結果逐值比對，確認定義一致。
2. 骨折限定 40 歲以後：OSD030 系列（髖、腕、脊椎，每次骨折的年齡；80 = 80 歲以上；77777／99999 = 拒答／不知道）。
   有任一次骨折發生在 ≥ 40 歲 → 1；從未骨折，或每一次骨折的年齡都已知且 < 40 歲 → 0；
   其餘 → 缺失。「每一次」依 OSQ020（該部位骨折次數）逐次檢查 OSD030 第 1…k 次的年齡；次數不明（如 9999）也算年齡不明
   （2026-10-28 第二輪審閱：原本只檢查第一次骨折的年齡）。
3. 缺失人數（Table 1 註腳，STROBE 14b）。
   不一致 vs 一致正常的盛行率比（粗、年齡校正、完整校正），與 analysis2.py 的骨折史模型相同。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

import analysis as A
import analysis2 as B
from analysis import Design, OUT, D

NEW = ["No periods, other reason", "Not ascertained"]
RACE5 = {3: "Non-Hispanic White", 4: "Non-Hispanic Black", 1: "Mexican American", 2: "Other Hispanic", 5: "Other, including multiracial"}


def frac40(w):
    parts = []
    for c in A.CYC:
        o = pd.read_sas(D / f"OSQ_{c}.xpt")
        cols = [x for x in o.columns if x.startswith(("OSD030", "OSQ020"))]
        parts.append(o[["SEQN"] + cols])
    o = pd.concat(parts, ignore_index=True)
    o = w[["SEQN", "OSQ010A", "OSQ010B", "OSQ010C"]].merge(o, on="SEQN", how="left")
    assert (o.SEQN.values == w.SEQN.values).all()
    out = np.full(len(o), np.nan)
    letters = "ABCDEFGHIJ"
    for i, rd in enumerate(o.to_dict("records")):
        sites = {s: rd["OSQ010" + s] for s in "ABC"}
        if all(v == 2 for v in sites.values()):
            out[i] = 0.0; continue
        if not any(v == 1 for v in sites.values()):
            continue                                   # 骨折題本身缺失
        known, unknown = [], False
        valid = lambda a: a is not None and a == a and 1 <= a <= 80
        for s, v in sites.items():
            if v != 1:
                continue
            k = rd.get("OSQ020" + s)
            k = int(k) if (k is not None and k == k and 1 <= k <= len(letters)) else None
            if k is None:                              # 骨折次數不明：收集已知年齡，但不能判為 0
                unknown = True
                known += [rd.get(f"OSD030{s}{L}") for L in letters if valid(rd.get(f"OSD030{s}{L}"))]
                continue
            for L in letters[:k]:                      # 第 1…k 次骨折，每一次的年齡都要已知
                a = rd.get(f"OSD030{s}{L}")
                if valid(a):
                    known.append(a)
                else:
                    unknown = True
        if any(a >= 40 for a in known):
            out[i] = 1.0
        elif not unknown:
            out[i] = 0.0
    return out


def main():
    w, _ = B.load_all()
    des = Design(w); q = stats.t.ppf(0.975, des.df)
    ages = w.RIDAGEYR.values.astype(float); prim = (ages >= 45) & (ages <= 55)
    w["Z"], w["tcat"] = B.site_scores(w, "FN")
    R = {}

    # ---- 1. Table 1 兩欄
    def col(m):
        return {"n": int(m.sum()), "age": des.mean(w.RIDAGEYR.values, m), "bmi": des.mean(w.BMXBMI.values, m),
                "race": {lab: des.prop((w.RIDRETH1.values == code).astype(float), m) for code, lab in RACE5.items()},
                "smoke": des.prop(w.smoke.values, m), "frac": des.prop(w.frac.values, m), "par_hip": des.prop(w.par_hip.values, m),
                "hrt": des.prop((w.RHQ540 == 1).astype(float).values, m & w.RHQ540.isin([1, 2]).values),
                "fn": des.mean(w.DXXNKBMD.values, m), "ls": des.mean(w.LS.values, m),
                "n_ls": int((m & w.LS.notna().values).sum()),
                "tcat": {str(c): des.prop((w.tcat.values == c).astype(float), m) for c in (0, 1, 2)},
                "zlow": des.prop((w.Z.values <= -2.0).astype(float), m)}
    # 定義一致性：重算 Menstruating 欄，與既有結果比對
    R1 = json.load(open(OUT / "results.json")); R2 = json.load(open(OUT / "results2.json")); R3 = json.load(open(OUT / "results3.json"))
    chk = col(prim & (w.status.values == "Menstruating")); s = "Menstruating"
    pairs = [(chk["n"], R1["table1"][s]["n"]), (chk["age"]["m"], R1["table1"][s]["age"]["m"]), (chk["bmi"]["m"], R1["table1"][s]["bmi"]["m"]),
             (chk["hrt"]["p"], R1["table1"][s]["hrt"]["p"]), (chk["fn"]["m"], R1["table1"][s]["fn"]["m"]),
             (chk["zlow"]["p"], R1["table1"][s]["zlow"]["p"]), (chk["tcat"]["1"]["p"], R1["table1"][s]["tcat"]["1"]["p"]),
             (chk["smoke"]["p"], R2["table1_extra"][s]["smoke"]["p"]), (chk["frac"]["p"], R2["table1_extra"][s]["frac"]["p"]),
             (chk["par_hip"]["p"], R2["table1_extra"][s]["par_hip"]["p"]), (chk["ls"]["m"], R2["table1_extra"][s]["ls"]["m"]),
             (chk["race"]["Mexican American"]["p"], R3["race5_table1"][s]["Mexican American"]["p"])]
    for a, b in pairs:
        assert abs(a - b) < 1e-9, (a, b)
    R["table1_more"] = {s: col(prim & (w.status.values == s)) for s in NEW}

    # ---- 2. 40 歲以後的骨折
    w["frac40"] = frac40(w)
    y = w.frac40.values
    R["frac40_counts"] = {"n": int((prim & ~np.isnan(y)).sum()), "cases": int(np.nansum(y[prim])),
                          "frac_any_cases": int(np.nansum(w.frac.values[prim])),
                          "missing_age": int((prim & (w.frac.values == 1) & np.isnan(y)).sum())}
    expo = np.where((w.tcat == 0) & (w.Z > -2), 0.0, np.where((w.tcat >= 1) & (w.Z > -2), 1.0, np.nan))
    race_d = pd.get_dummies(w.race4).drop(columns="Non-Hispanic White").astype(float).values
    bmi = w.BMXBMI.values
    models = {"Crude": np.column_stack([expo]),
              "Age-adjusted": np.column_stack([expo, ages]),
              "Adjusted for age, BMI, race and ethnicity, hormone use": np.column_stack([expo, ages, bmi, race_d, w.hrt.values])}
    # 與 analysis2.py 同一模型：先用骨折史（任何年齡）重算一次，確認和既有 PR 一致
    b, V, dm = B.svyglm(des, models["Age-adjusted"], w.frac.values, prim & ~np.isnan(expo), family="poisson")
    assert abs(np.exp(b[1]) - R2["pr"]["frac"]["Age-adjusted"]["pr"]) < 1e-9
    R["pr_frac40"] = {}
    for k, Xm in models.items():
        b, V, dm = B.svyglm(des, Xm, y, prim & ~np.isnan(expo), family="poisson")
        R["pr_frac40"][k] = {"pr": float(np.exp(b[1])), "lo": float(np.exp(b[1] - q * np.sqrt(V[1, 1]))),
                             "hi": float(np.exp(b[1] + q * np.sqrt(V[1, 1]))), "n": int(dm.sum())}
    P = w[prim]
    R["missing"] = {"bmi": int(P.BMXBMI.isna().sum()), "hormone": int(P.hrt.isna().sum()),
                    "hormone_not_ascertained": int((P.hrt.isna() & (P.status == "Not ascertained")).sum()),
                    "smoking": int(P.smoke.isna().sum()), "fracture": int(P.frac.isna().sum()),
                    "parent_hip": int(P.par_hip.isna().sum()), "told": int(P.told.isna().sum()), "lumbar_spine": int(P.LS.isna().sum())}
    # 完整校正模型（骨折史）排除的人：其中多少是狀態未確定（荷爾蒙題沒問到）
    dom = prim & ~np.isnan(expo) & w.frac.notna().values
    drop = dom & (w.hrt.isna() | w.BMXBMI.isna()).values
    R["full_model_excluded"] = {"n_domain": int(dom.sum()), "n": int(drop.sum()),
                                "not_ascertained": int((drop & (w.status.values == "Not ascertained")).sum())}
    json.dump(R, open(OUT / "results5.json", "w"), indent=1, default=float, ensure_ascii=False)
    print("missing:", R["missing"], "full model excluded:", R["full_model_excluded"])
    print("Table 1 new columns:", {s: R["table1_more"][s]["n"] for s in NEW})
    print("fracture ≥40:", R["frac40_counts"])
    for k, v in R["pr_frac40"].items():
        print(f"  {k}: PR {v['pr']:.2f} ({v['lo']:.2f}–{v['hi']:.2f}), n = {v['n']}")


if __name__ == "__main__":
    main()
