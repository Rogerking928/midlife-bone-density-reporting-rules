"""正式分析（2026-10-27）：股骨頸 DXA 依「停經後規則（T 分數／WHO 分類）」與「停經前規則（Z 分數，ISCD）」判讀，兩者不一致的比例。

事前（2026-10-23，讀 BMD 之前）定的：主要量＝45–55 歲女性中 T 路徑判為異常（骨質缺乏或骨質疏鬆）而 Z 路徑判為同齡範圍內的比例；
停損 < 2%；依年齡層與問卷停經狀態分層。其餘（調查設計 CI、骨折結局、敏感度分析、子宮切除分組）為事後，列在 Table S1。

T 路徑：Looker 2010 (JBMR 25:64) 依 NHANES III 20–29 歲非西裔白人女性換算的 BMD 切點（股骨頸 ≤ 0.56 骨質疏鬆、≤ 0.74 低骨量；
全髖 ≤ 0.64、≤ 0.82）。ISCD 2023：T 分數一律用白人女性參考、不分族裔。
Z 路徑：同族裔（RIDRETH1）、同 5 歲年齡層的加權均值與 SD（2005–2018 合併）；Z ≤ −2.0 =「低於同齡預期」（ISCD 2023）。
合併權重：各年齡有股骨 DXA 的週期數不同（2013–14 只做 ≥ 40 歲、2017–18 只做 ≥ 50 歲），權重除以該年齡實際貢獻的週期數。
CI：Taylor 線性化（strata × PSU）＋ Korn–Graubard 比例區間，df ＝ PSU 數 − strata 數。
"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D, OUT = ROOT / "data", ROOT / "out"
OUT.mkdir(exist_ok=True)
CYC = {"D": "2005–2006", "E": "2007–2008", "F": "2009–2010", "H": "2013–2014", "J": "2017–2018"}
FN_OP, FN_LOW = 0.56, 0.74          # Looker 2010 股骨頸
TH_OP, TH_LOW = 0.64, 0.82          # Looker 2010 全髖
T_MEAN, T_SD = 0.86, 0.12           # 由上面兩個切點反推（0.74 = μ − σ、0.56 = μ − 2.5σ），只用來畫連續 T 值


# ------------------------------------------------------------------ 資料
def load():
    flow, rows = [], []
    for c in CYC:
        d = pd.read_sas(D / f"DEMO_{c}.xpt")
        f = pd.read_sas(D / f"DXXFEM_{c}.xpt")[["SEQN", "DXAFMRST", "DXXNKBMD", "DXXOFBMD"]]
        r = pd.read_sas(D / f"RHQ_{c}.xpt")
        o = pd.read_sas(D / f"OSQ_{c}.xpt")
        b = pd.read_sas(D / f"BMX_{c}.xpt")[["SEQN", "BMXBMI"]]
        r = r.rename(columns={"RHD042": "reason", "RHD043": "reason"})
        if c == "D":
            r["bso"] = np.where((r.RHQ300 == 1) & (r.RHQ310 == 1), 1, np.where(r.RHQ300.isin([1, 2]), 2, np.nan))
        else:
            r["bso"] = r.RHQ305
        r = r[["SEQN", "RHQ031", "reason", "RHD280", "bso", "RHQ540"]]
        o = o[["SEQN", "OSQ010A", "OSQ010B", "OSQ010C", "OSQ060"]]
        x = d[["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "WTMEC2YR", "SDMVSTRA", "SDMVPSU"]]
        x = x.merge(f, on="SEQN", how="left").merge(r, on="SEQN", how="left").merge(o, on="SEQN", how="left").merge(b, on="SEQN", how="left")
        x["cycle"] = c
        rows.append(x)
    x = pd.concat(rows, ignore_index=True)
    # 流程（女性 40–59）
    s = x[x.RIAGENDR == 2]; flow.append(("Women, 5 cycles", len(s)))
    s = s[s.RIDAGEYR.between(40, 59)]; flow.append(("Aged 40–59 y", len(s)))
    noscan_age = s[(s.cycle == "J") & (s.RIDAGEYR < 50)]
    flow.append(("Excluded: 2017–2018 femur DXA limited to ≥50 y", len(noscan_age)))
    s = s.drop(noscan_age.index)
    for code, lab in ((3, "not scanned, pregnant"), (4, "not scanned, weight over DXA table limit"), (5, "not scanned, other reason"),
                      (2, "scan invalid")):
        flow.append((f"Excluded: {lab}", int((s.DXAFMRST == code).sum())))
    miss = s.DXAFMRST.isna() | ~s.DXAFMRST.isin([1, 2, 3, 4, 5])
    flow.append(("Excluded: no femur scan record", int(miss.sum())))
    s = s[s.DXAFMRST == 1]
    flow.append(("Excluded: femoral neck BMD missing", int(s.DXXNKBMD.isna().sum())))
    s = s[s.DXXNKBMD.notna()]
    flow.append(("Analysed, aged 40–59 y", len(s)))
    flow.append(("Primary sample, aged 45–55 y", int(s.RIDAGEYR.between(45, 55).sum())))
    return x, flow


def prepare(x):
    w = x[(x.RIAGENDR == 2) & (x.DXAFMRST == 1) & x.DXXNKBMD.notna() & (x.WTMEC2YR > 0)].copy()
    w = w[~((w.cycle == "J") & (w.RIDAGEYR < 50))]
    k = np.where(w.RIDAGEYR >= 50, 5, np.where(w.RIDAGEYR >= 40, 4, 3))      # 該年齡有股骨 DXA 的週期數
    w["w"] = w.WTMEC2YR / k
    w["age5"] = (w.RIDAGEYR // 5 * 5).astype(int)
    return w


def zref(w, col, by=("RIDRETH1", "age5"), lo=20, hi=79):
    g = w[w.RIDAGEYR.between(lo, hi) & w[col].notna()].groupby(list(by))
    m = g.apply(lambda t: np.average(t[col], weights=t.w))
    sd = g.apply(lambda t: np.sqrt(np.average((t[col] - np.average(t[col], weights=t.w)) ** 2, weights=t.w)))
    return pd.DataFrame({"m": m, "sd": sd})


def add_scores(w, site="FN", zby=("RIDRETH1", "age5"), tcut=None):
    col, (op, low) = ("DXXNKBMD", (FN_OP, FN_LOW)) if site == "FN" else ("DXXOFBMD", (TH_OP, TH_LOW))
    if tcut: op, low = tcut
    if zby == ("age10",):
        w = w.assign(age10=(w.RIDAGEYR // 10 * 10).astype(int))
    if zby == ("age5",) or zby == ("age10",):
        ref = zref(w, col, by=zby)
    else:
        ref = zref(w, col, by=zby)
    j = w.join(ref, on=list(zby))
    Z = (j[col] - j.m) / j.sd
    tcat = np.where(j[col] <= op, 2, np.where(j[col] <= low, 1, 0))       # 0 正常、1 低骨量、2 骨質疏鬆
    tcat = np.where(j[col].isna(), np.nan, tcat)
    return Z.values, tcat


def status(w):
    hyst = w.RHD280 == 1
    bso = w.bso == 1
    return np.select(
        [bso, hyst & ~bso, w.RHQ031 == 1, (w.RHQ031 == 2) & (w.reason == 7), w.RHQ031 == 2],
        ["Bilateral oophorectomy", "Hysterectomy, ovaries conserved", "Menstruating", "No periods, menopause", "No periods, other reason"],
        "Not ascertained")


# ------------------------------------------------------------------ 調查設計估計
class Design:
    def __init__(self, w):
        self.w = w.w.values
        self.h = w.SDMVSTRA.values
        self.psu = (w.SDMVSTRA.astype(int).astype(str) + "_" + w.SDMVPSU.astype(int).astype(str)).values
        ps = pd.DataFrame({"h": self.h, "p": self.psu}).drop_duplicates()
        self.df = ps.p.nunique() - ps.h.nunique()

    def _var(self, z):
        t = pd.DataFrame({"h": self.h, "p": self.psu, "z": z}).groupby(["h", "p"]).z.sum().reset_index()
        v = 0.0
        for _, g in t.groupby("h"):
            n = len(g)
            if n > 1:
                v += n / (n - 1) * ((g.z - g.z.mean()) ** 2).sum()
        return v

    def lin(self, y, dom):
        """比例（或均值）及其線性化變數。y 可為 0/1 或連續；dom 為布林遮罩。"""
        dom = np.asarray(dom, bool) & ~np.isnan(np.asarray(y, float))
        y = np.nan_to_num(np.asarray(y, float))
        N = (self.w * dom).sum()
        p = (self.w * dom * y).sum() / N
        z = self.w * dom * (y - p) / N
        return p, z, dom

    def prop(self, y, dom):
        p, z, dom = self.lin(y, dom)
        v = self._var(z)
        n = int(dom.sum())
        lo, hi = korn_graubard(p, v, n, self.df)
        return {"p": 100 * p, "lo": 100 * lo, "hi": 100 * hi, "se": 100 * np.sqrt(v), "n": n,
                "n_events": int(np.round((np.asarray(y, float)[dom] > 0).sum())) if set(np.unique(np.nan_to_num(y))) <= {0, 1} else None,
                "N_pop": float((self.w * dom).sum()), "N_pop_events": float((self.w * dom * np.nan_to_num(y)).sum())}

    def mean(self, y, dom):
        p, z, dom = self.lin(y, dom)
        se = np.sqrt(self._var(z)); t = stats.t.ppf(0.975, self.df)
        return {"m": p, "lo": p - t * se, "hi": p + t * se, "n": int(dom.sum())}

    def diff(self, y, d1, d2):
        p1, z1, _ = self.lin(y, d1); p2, z2, _ = self.lin(y, d2)
        se = np.sqrt(self._var(z1 - z2)); t = stats.t.ppf(0.975, self.df)
        d = p1 - p2
        return {"d": 100 * d, "lo": 100 * (d - t * se), "hi": 100 * (d + t * se), "p": float(2 * stats.t.sf(abs(d / se), self.df)) if se > 0 else np.nan}


def korn_graubard(p, v, n, df):
    """Korn & Graubard (1998) 區間；p = 0 時有效樣本數取 n。"""
    neff = p * (1 - p) / v if v > 0 and 0 < p < 1 else n
    neff = min(n, neff * (stats.t.ppf(0.975, n - 1) / stats.t.ppf(0.975, df)) ** 2) if n > 1 else n
    x = p * neff
    lo = stats.beta.ppf(0.025, x, neff - x + 1) if x > 0 else 0.0
    hi = stats.beta.ppf(0.975, x + 1, neff - x) if x < neff else 1.0
    return float(lo), float(hi)


# ------------------------------------------------------------------ 主程式
def flips(w, Z, tcat):
    within = Z > -2.0
    return {"any": ((tcat >= 1) & within).astype(float), "op": ((tcat == 2) & within).astype(float),
            "rev": ((tcat == 0) & (Z <= -2.0)).astype(float)}


def main():
    x, flow = load()
    w = prepare(x)
    w["status"] = status(w)
    Z, tcat = add_scores(w)
    w["Z"], w["tcat"], w["T"] = Z, tcat, (w.DXXNKBMD - T_MEAN) / T_SD
    F = flips(w, Z, tcat)
    for k, v in F.items():
        w["flip_" + k] = v
    des = Design(w)
    a40 = w.RIDAGEYR.between(40, 59); prim = w.RIDAGEYR.between(45, 55)
    R = {"flow": flow, "df": des.df, "n_psu_strata": None}

    # 主要量與交叉表（45–55）
    R["primary"] = {k: des.prop(F[k], prim) for k in F}
    cells = {}
    for tc, tl in ((0, "Normal"), (1, "Low bone mass"), (2, "Osteoporosis")):
        for zc, zl, m in (("within", "Within expected range for age (Z > −2.0)", w.Z > -2.0), ("below", "Below expected range for age (Z ≤ −2.0)", w.Z <= -2.0)):
            cells[f"{tl}|{zc}"] = des.prop(((w.tcat == tc) & m).astype(float), prim)
    R["xtab"] = cells

    # 年齡層
    R["age"] = {}
    for lo, hi in ((40, 44), (45, 49), (50, 54), (55, 59)):
        m = w.RIDAGEYR.between(lo, hi)
        R["age"][f"{lo}–{hi}"] = {k: des.prop(F[k], m) for k in ("any", "op")}
    R["age1"] = {int(a): {k: des.prop(F[k], w.RIDAGEYR == a) for k in ("any", "op")} for a in range(40, 60)}

    # 停經狀態（45–55）
    ORDER = ["Menstruating", "No periods, menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy", "No periods, other reason", "Not ascertained"]
    R["status"] = {s: {k: des.prop(F[k], prim & (w.status == s)) for k in ("any", "op")} for s in ORDER}
    R["status_share"] = {s: des.prop((w.status == s).astype(float), prim) for s in ORDER}

    # Table 1：45–55 依狀態
    T1 = {}
    for s in ["All"] + ORDER[:4]:
        m = prim if s == "All" else prim & (w.status == s)
        T1[s] = {"n": int(m.sum()), "age": des.mean(w.RIDAGEYR, m), "bmi": des.mean(w.BMXBMI, m),
                 "fn": des.mean(w.DXXNKBMD, m), "Z": des.mean(w.Z, m),
                 "race": {r: des.prop((w.RIDRETH1 == r).astype(float), m) for r in (3, 4, 1, 2, 5)},
                 "hrt": des.prop((w.RHQ540 == 1).astype(float), m & w.RHQ540.isin([1, 2])),
                 "tcat": {c: des.prop((w.tcat == c).astype(float), m) for c in (0, 1, 2)},
                 "zlow": des.prop((w.Z <= -2.0).astype(float), m)}
    R["table1"] = T1

    # 結局：骨折史、曾被告知骨質疏鬆（45–55），依一致性分組
    grp = np.select([(w.tcat == 0) & (w.Z > -2), (w.tcat >= 1) & (w.Z > -2), w.Z <= -2],
                    ["T normal", "T abnormal, Z within range", "Z below expected range"], "")
    w["grp"] = grp
    frac = np.where(w[["OSQ010A", "OSQ010B", "OSQ010C"]].isin([1]).any(axis=1), 1.0,
                    np.where(w[["OSQ010A", "OSQ010B", "OSQ010C"]].isin([2]).all(axis=1), 0.0, np.nan))
    told = np.where(w.OSQ060 == 1, 1.0, np.where(w.OSQ060 == 2, 0.0, np.nan))
    w["frac"], w["told"] = frac, told
    G = ["T normal", "T abnormal, Z within range", "Z below expected range"]
    R["outcome"] = {o: {g: des.prop(w[o].values, prim & (w.grp == g)) for g in G} for o in ("frac", "told")}
    R["outcome_diff"] = {o: {g: des.diff(w[o].values, prim & (w.grp == g), prim & (w.grp == "T normal")) for g in G[1:]} for o in ("frac", "told")}
    R["outcome_age"] = {g: des.mean(w.RIDAGEYR, prim & (w.grp == g)) for g in G}
    # 年齡標準化（直接法，45–49／50–55 兩層，標準人口＝45–55 全體加權）：骨折史是累積量，分組的平均年齡不同
    bands = [w.RIDAGEYR.between(45, 49), w.RIDAGEYR.between(50, 55)]
    sh = np.array([(des.w * b.values).sum() for b in bands]); sh = sh / sh.sum()
    t975 = stats.t.ppf(0.975, des.df)
    def std(o, g):
        ps, zs = zip(*[des.lin(w[o].values, prim & b & (w.grp == g))[:2] for b in bands])
        return sum(s * p for s, p in zip(sh, ps)), sum(s * z for s, z in zip(sh, zs))
    R["outcome_std"] = {}
    for o in ("frac", "told"):
        R["outcome_std"][o] = {}
        p0, z0 = std(o, "T normal")
        for g in G[:2]:
            p, z = std(o, g)
            se = np.sqrt(des._var(z)); sed = np.sqrt(des._var(z - z0))
            R["outcome_std"][o][g] = {"p": 100 * p, "lo": 100 * (p - t975 * se), "hi": 100 * (p + t975 * se),
                                     "d": 100 * (p - p0), "dlo": 100 * (p - p0 - t975 * sed), "dhi": 100 * (p - p0 + t975 * sed),
                                     "pval": float(2 * stats.t.sf(abs((p - p0) / sed), des.df))}
    R["std_shares"] = list(sh)

    # 敏感度分析（全部針對 45–55 主要量與骨質疏鬆翻轉）
    S = {}
    def put(name, Zs, tc, dom=prim, d=des):
        f = flips(w, Zs, tc)
        S[name] = {"any": d.prop(f["any"], dom), "op": d.prop(f["op"], dom)}
    put("Primary analysis", Z, tcat)
    Zh, th = add_scores(w, site="TH"); put("Total hip instead of femoral neck", Zh, th)
    lowest_t = np.fmax(tcat, th); lowest_z = np.fmin(Z, Zh)
    put("Lower of femoral neck and total hip", lowest_z, lowest_t)
    y = w[(w.RIDRETH1 == 3) & w.RIDAGEYR.between(20, 29)]
    mu = np.average(y.DXXNKBMD, weights=y.w); sd = np.sqrt(np.average((y.DXXNKBMD - mu) ** 2, weights=y.w))
    R["internal_tref"] = {"mean": mu, "sd": sd, "n": len(y)}
    _, t_int = add_scores(w, tcut=(mu - 2.5 * sd, mu - sd)); put("T reference from 2005–2010 white women aged 20–29 y", Z, t_int)
    Zr, _ = add_scores(w, zby=("age5",)); put("Z reference not race-specific", Zr, tcat)
    Z10, _ = add_scores(w, zby=("RIDRETH1", "age10")) if False else (None, None)
    w["age10"] = (w.RIDAGEYR // 10 * 10).astype(int)
    ref10 = zref(w, "DXXNKBMD", by=("RIDRETH1", "age10"))
    j = w.join(ref10, on=["RIDRETH1", "age10"]); put("Z reference by 10-year age band", ((j.DXXNKBMD - j.m) / j.sd).values, tcat)
    put("Excluding ever-users of female hormones", Z, tcat, dom=prim & (w.RHQ540 == 2))
    put("Cycles 2005–2014 only", Z, tcat, dom=prim & (w.cycle != "J"))
    class Unw(Design):
        pass
    du = Design(w.assign(w=1.0)); put("Unweighted", Z, tcat, d=du)
    R["sens"] = S

    # 族裔
    R["race"] = {r: {k: des.prop(F[k], prim & (w.RIDRETH1 == r)) for k in ("any", "op")} for r in (3, 4, 1, 2, 5)}

    # 圖用的個人層資料
    w[prim][["SEQN", "cycle", "RIDAGEYR", "RIDRETH1", "status", "DXXNKBMD", "T", "Z", "tcat", "w", "frac", "told"]].to_csv(OUT / "women45_55.csv", index=False)
    # Z 參考值（非西裔白人），給機制圖
    ref = zref(w, "DXXNKBMD")
    R["zref_nhw"] = {int(a): {"m": float(r.m), "sd": float(r.sd)} for (rr, a), r in ref.iterrows() if rr == 3}
    json.dump(R, open(OUT / "results.json", "w"), indent=1, default=float, ensure_ascii=False)
    pr = R["primary"]
    print("df", des.df, "| 45–55 n", pr["any"]["n"], "any %.1f (%.1f–%.1f)" % (pr["any"]["p"], pr["any"]["lo"], pr["any"]["hi"]),
          "op %.2f (%.2f–%.2f)" % (pr["op"]["p"], pr["op"]["lo"], pr["op"]["hi"]), "rev %.2f" % pr["rev"]["p"])
    for k, v in R["status"].items():
        print(f"  {k:35s} n={v['any']['n']:5d} any {v['any']['p']:.1f} ({v['any']['lo']:.1f}–{v['any']['hi']:.1f}) op {v['op']['p']:.2f}")
    for k, v in R["age"].items():
        print(f"  age {k}: n={v['any']['n']} any {v['any']['p']:.1f} op {v['op']['p']:.2f}")
    for k, v in S.items():
        print(f"  S {k:55s} n={v['any']['n']:5d} any {v['any']['p']:.1f} ({v['any']['lo']:.1f}–{v['any']['hi']:.1f}) op {v['op']['p']:.2f}")
    for o in ("frac", "told"):
        print(o, {g: (round(v["p"], 1), v["n"]) for g, v in R["outcome"][o].items()}, {g: round(v["d"], 1) for g, v in R["outcome_diff"][o].items()}, {g: round(v["p"], 3) for g, v in R["outcome_diff"][o].items()})
    print("flow", flow)
    print("internal T ref", R["internal_tref"])


if __name__ == "__main__":
    main()
