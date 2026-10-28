"""第三版圖表本（第二輪外部審閱）：Word ＋ 700 dpi PNG ＋ Source_data.xlsx，同一個交付資料夾（版式照 Bcc，共用 book_base）。

正文：Table 1–2、Figure 1–3。補充：Table S1–S3、Figure S1–S5。**只交 PNG，不產 PDF。**
用法：python3 code/book2.py [輸出資料夾]
"""
import os
import json, os, shutil, sys
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font as XFont

from book_base import OUT, build, numeric_table

DEST = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("DELIVERABLE_DIR", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "deliverable"))
R1 = json.load(open(OUT / "results.json"))
R = json.load(open(OUT / "results2.json"))
R3 = json.load(open(OUT / "results3.json"))
R4 = json.load(open(OUT / "results4.json"))
ci = lambda d, nd=1: f"{d['p']:.{nd}f} ({d['lo']:.{nd}f} to {d['hi']:.{nd}f})"
mci = lambda d, nd=1: f"{d['m']:.{nd}f} ({d['lo']:.{nd}f} to {d['hi']:.{nd}f})"
mil = lambda d: f"{d['t']:.2f} ({d['lo']:.2f} to {d['hi']:.2f})"
FLOW = dict(R1["flow"])
DESIGN = ("Estimates are weighted with NHANES examination weights, divided by the number of cycles that offered femur DXA at each age "
          "(2005–2010, all ages; 2013–2014, ≥40 y; 2017–2018, ≥50 y). 95% CIs use Taylor linearization with the masked strata and primary "
          "sampling units (Korn–Graubard for proportions).")
RULES = ("Postmenopausal rule: WHO category from the T-score (white female reference aged 20–29 years; low bone mass, −1.0 to −2.5; osteoporosis, "
         "−2.5 or lower). Premenopausal rule: Z-score against women of the same race and ethnicity and 5-year age band; −2.0 or lower is below "
         "the expected range for age (ISCD 2023 Official Positions).")
TREF = R["tref_internal"]
LSREF = (f"Lumbar spine (L1–L4) T-scores use the weighted mean ({TREF['LS']['mean']:.3f} g/cm²) and SD ({TREF['LS']['sd']:.3f}) of non-Hispanic White "
         f"women aged 20–29 years in NHANES 2005–2010 (n = {TREF['LS']['n']}), because the manufacturer reference values used clinically are not publicly documented.")

# ------------------------------------------------------------------ Table 1（精簡）
ST = ["All", "Menstruating", "No periods, menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy"]
HEAD = ["Characteristic", "All women aged 45–55 y", "Menstruating", "No periods, natural menopause", "Hysterectomy, ovaries conserved", "Bilateral oophorectomy"]
T1, TX, SH = R1["table1"], R["table1_extra"], R1["status_share"]
pc = lambda d: f"{d['p']:.1f}"
rows = [HEAD,
        ["Participants, n"] + [f"{T1[s]['n']:,}" for s in ST],
        ["Weighted share of women aged 45–55 y, %"] + ["100"] + [pc(SH[s]) for s in ST[1:]],
        ["Age, y, mean (95% CI)"] + [mci(T1[s]["age"]) for s in ST],
        ["Body mass index, kg/m², mean (95% CI)"] + [mci(T1[s]["bmi"]) for s in ST],
        ["Race and ethnicity, %", "", "", "", "", ""]]
for r in ("Non-Hispanic White", "Non-Hispanic Black", "Mexican American", "Other Hispanic", "Other, including multiracial"):
    rows.append([f"  {r}"] + [pc(R3["race5_table1"][s][r]) for s in ST])
rows += [["Current smoking, %"] + [pc(TX[s]["smoke"]) for s in ST],
         ["Prior hip, wrist or spine fracture, %"] + [pc(TX[s]["frac"]) for s in ST],
         ["Parent ever fractured a hip, %"] + [pc(TX[s]["par_hip"]) for s in ST],
         ["Ever used female hormones, %"] + [pc(T1[s]["hrt"]) for s in ST],
         ["Femoral neck BMD, g/cm², mean (95% CI)"] + [mci(T1[s]["fn"], 3) for s in ST],
         ["Lumbar spine BMD, g/cm², mean (95% CI)"] + [mci(TX[s]["ls"], 3) for s in ST],
         ["Postmenopausal rule, femoral neck, %", "", "", "", "", ""]]
for c, lab in (("0", "Normal"), ("1", "Low bone mass"), ("2", "Osteoporosis")):
    rows.append([f"  {lab}"] + [pc(T1[s]["tcat"][c]) for s in ST])
rows.append(["Premenopausal rule, femoral neck Z-score −2.0 or lower, %"] + [pc(T1[s]["zlow"]) for s in ST])
nO, nNA = R1["status"]["No periods, other reason"]["any"]["n"], R1["status"]["Not ascertained"]["any"]["n"]
TAB1 = {"kind": "table", "caption": "Table 1. Characteristics of women aged 45–55 years with femoral neck DXA, by menstrual and surgical status",
        "rows": rows, "first_w": 2.4,
        "note": (f"Unit of analysis: participant (one examination each). Status hierarchy: bilateral oophorectomy; hysterectomy with at least one ovary conserved; "
                 f"at least one menstrual period in the past 12 months; no periods because of menopause. The all-women column also includes {nO} women with "
                 f"no periods for another reason and {nNA} whose status was not ascertained (Figure 2D). Race and ethnicity groups are those used for the Z-score reference. "
                 f"Percentages exclude missing responses; lumbar spine BMD was available for {TX['All']['n_ls']:,} women. {DESIGN} BMD, bone mineral density.")}

# ------------------------------------------------------------------ Table 2（精簡：百萬人＋CI）
X, X3 = R3["xtab"]["FN"], R3["xtab"]["LOW3"]
rows = [["Postmenopausal rule → premenopausal rule", "Femoral neck: % (95% CI)", "Femoral neck: US women, millions (95% CI)",
         "Lowest of three sites: % (95% CI)", "Lowest of three sites: US women, millions (95% CI)"]]
for lab in ("Normal", "Low bone mass", "Osteoporosis"):
    for z, zl in (("within", "within expected range"), ("below", "below expected range")):
        a, b = X[f"{lab}|{z}"], X3[f"{lab}|{z}"]
        rows.append([f"{lab} → {zl}", ci(a), mil(a["mil"]), ci(b), mil(b["mil"])])
rows.append(["Discordant: low bone mass or osteoporosis → within expected range", ci(X["primary"]), mil(X["primary"]["mil"]), ci(X3["primary"]), mil(X3["primary"]["mil"])])
CN = R3["counts"]
rows.append(["Discordant with structurally indeterminate menopausal status: upper bound", "", mil(CN["45–55|upper"]), "", ""])
rows.append(["Discordant with structurally indeterminate menopausal status: expected premenopausal", "", mil(CN["45–55|expected"]), "", ""])
TAB2 = {"kind": "table", "caption": "Table 2. Cross-classification of DXA results under the postmenopausal and premenopausal reporting rules, women aged 45–55 years",
        "rows": rows, "first_w": 2.6,
        "note": (f"n = {X['primary']['n']:,} (femoral neck available for all; lumbar spine for {R['n_ls_45_55']:,}). Femoral neck is the pre-specified primary site. "
                 f"Lowest of three sites: worst T-score category and lowest Z-score across femoral neck, total hip and lumbar spine, as in clinical reporting; "
                 f"women without a valid spine scan are classified from the hip sites. Weighted total: {X['N']:.1f} million US women aged 45–55 years, averaged "
                 f"over the survey years; millions are the weighted proportion (Korn–Graubard 95% CI) multiplied by this total. "
                 "Structurally indeterminate: hysterectomy with ovaries conserved, no periods for a reason other than menopause, or status not ascertained. "
                 "Upper bound counts every such woman with a discordant category; the expected count weights each by her age-specific probability of still being "
                 "premenopausal, estimated among women with an intact uterus and ovaries (Figure 2C); 95% CIs for these two rows from 1,000 Rao–Wu bootstrap replicates "
                 f"that refit the probability model. The ISCD 2023 Official Positions prefer Z-scores before menopause and state that the WHO criteria may be applied "
                 "in the menopausal transition, so the choice of rule in that window is permitted rather than prescribed. Women in the transition who still "
                 f"menstruate are not counted here, which makes the expected count conservative; sensitivity in Table S3. {RULES} {LSREF} {DESIGN}")}

# ------------------------------------------------------------------ Figures
FIG1 = {"kind": "figure", "caption": "Figure 1. How the two reporting rules diverge in midlife", "path": OUT / "Figure_1.png",
        "note": ("(A) Worked example; Z-score from the non-Hispanic White reference for ages 50–54 years. "
                 f"(B) Weighted distribution of femoral neck BMD by age, women aged 40–59 years (n = {FLOW['Analysed, aged 40–59 y']:,}); kernel density within each "
                 "year of age, averaged with adjacent years; darker shading holds a larger share of women of that age. Red lines: Z-score −2.0 boundary for non-Hispanic "
                 "White (solid), non-Hispanic Black (dashed) and Mexican American (dotted) women. Amber band: discordant region for non-Hispanic White women "
                 "(T-score below −1.0, Z-score above −2.0); red band: from 50 years the White Z-score boundary lies below T-score −2.5, so a T-score in the osteoporosis range can still be within the expected range for age. Women aged 40–49 years come "
                 "from 2005–2014 and women aged 50–59 years from 2005–2018. (C) Cross-classification at the femoral neck, women aged 45–55 years; widths are "
                 "weighted US women in millions (Table 2); numbers may not sum because of rounding. (D) Osteoporosis by T-score with a Z-score within the expected range, femoral neck and lowest of three "
                 "sites (95% CI). Unit: participant.")}
FIG2 = {"kind": "figure", "caption": "Figure 2. Menopausal status, discordance and the number of women for whom the reporting rule is uncertain",
        "path": OUT / "Figure_2.png",
        "note": ("(A) Menstrual and surgical status by age, women aged 40–59 years; each status modelled by survey-weighted logistic regression with a "
                 "restricted cubic spline of age (4 knots) and rescaled to sum to 100%. Indeterminate: hysterectomy with ovaries conserved, no periods for a reason "
                 "other than menopause, or status not ascertained. (B) Discordant category (femoral neck), same model; points, single-year estimates (Figure S5). "
                 "(C) Weighted US women with a discordant category and structurally indeterminate status, 2-year age bands. Upper bound counts all such women; "
                 "expected weights each by the probability of still being premenopausal at her age, estimated among women with an intact uterus and ovaries. "
                 "Error bars, 95% CIs from 1,000 Rao–Wu bootstrap replicates; totals for ages 45–55 years in Table 2. "
                 "(D) Women aged 45–55 years; labels give participants and weighted share; teal marks the indeterminate groups. Unit: participant.")}
FIG3 = {"kind": "figure", "caption": "Figure 3. The postmenopausal and premenopausal rules compared against prior fracture, women aged 45–55 years",
        "path": OUT / "Figure_3.png",
        "note": (f"Self-reported doctor-diagnosed hip, wrist or spine fracture (n = {R3['n_frac']:,}; {R3['n_frac_cases']} with a fracture). "
                 "(A) Survey-weighted logistic regression on a restricted cubic spline of the femoral neck T-score (4 knots) and age, marginally standardized "
                 "to the sample age distribution; shown between the 1st and 99th T-score percentiles; few women have T-scores below −2, so the interval is wide there. "
                 "(B) Survey-weighted Poisson regression adjusted for age. (C) Weighted C-statistics for each score used alone; 95% CIs and differences from 1,000 "
                 f"Rao–Wu bootstrap replicates. Within this narrow age range T- and Z-scores are nearly collinear (r = {R4['corr_TZ_45_55']:.2f}), so similar "
                 f"continuous C-statistics are expected. Same positivity: the Z-score cut-point ({R4['z_matched_cut']:.2f}) that classifies the same weighted "
                 f"share of women as T-score −1.0 or lower ({R4['t_positivity']:.1f}% vs {R4['z_matched_positivity']:.1f}%). At equal positivity the two rules "
                 "discriminate equally: as a label for fracture risk, the −2.0 cut-point carries almost no information in women of this age. The Z-score −2.0 "
                 "threshold is intended to flag BMD below that expected for age, not to predict fracture, and self-reported prior fracture is a weak outcome "
                 "(C-statistics near 0.56). "
                 "Unit: participant.").replace("(-", "(−")}

# ------------------------------------------------------------------ Supplement
REG = [["Analysis", "Status", "Date fixed", "Exhibit"],
       ["Primary outcome: femoral neck T-score category abnormal while Z-score within expected range, women aged 45–55 y", "Pre-specified", "2026-10-23", "Table 2"],
       ["Stopping rule: abandon if primary proportion < 2%", "Pre-specified", "2026-10-23", f"Not met ({X['primary']['p']:.1f}%)"],
       ["Strata by age and questionnaire menstrual status", "Pre-specified", "2026-10-23", "Figure 2"],
       ["Osteoporosis label discordance (T-score −2.5 or lower, the treatment threshold, with Z-score within range)", "Pre-specified", "2026-10-23", "Table 2, Figures 1D, 2D"],
       ["Survey-design 95% CIs; counts as weighted proportion × total", "Post hoc", "2026-10-27", "All"],
       ["Status groups from hysterectomy and oophorectomy items", "Post hoc", "2026-10-27", "Table 1, Figure 2"],
       ["Lumbar spine and lowest-of-three-sites classification", "Post hoc", "2026-10-27", "Table 2, Figures 1D, S2"],
       ["Age splines; status composition; upper-bound and expected counts with indeterminate status", "Post hoc", "2026-10-27", "Figure 2, Table 2"],
       ["T- versus Z-score comparison against prior fracture (prevalence ratios, C-statistics)", "Post hoc", "2026-10-27", "Figure 3"],
       ["Prior fracture and diagnosed osteoporosis by concordance group", "Post hoc", "2026-10-27", "Table S4"],
       ["Specification curve and variance decomposition", "Post hoc", "2026-10-27", "Figure S2, Table S2"],
       ["Heterogeneity by participant characteristics", "Post hoc", "2026-10-27", "Figure S3"],
       ["Z-score dichotomized at the cut-point giving the same positivity as T-score −1.0, and at −1.0", "Post hoc", "2026-10-27", "Figure 3C"],
       ["Expected count with the premenopausal probability shifted 1 and 2 years earlier; late menopausal transition (2005–2006)", "Post hoc", "2026-10-27", "Table S3"],
       ["ISCD 2023 Official Positions checked against the original text: WHO criteria may be applied in the menopausal transition (permitted, not required)", "Post hoc", "2026-10-27", "Frames the choice of rule as discretionary"],
       ["T-score cut-points taken from published NHANES III femur thresholds", "Correction", "2026-10-27", "Earlier draft used an unverified reference mean and SD"],
       ["Pooled weights divided by the number of cycles offering femur DXA at each age", "Correction", "2026-10-27", "Earlier draft divided all weights by 5"]]
TABS1 = {"kind": "table", "caption": "Table S1. Analysis register", "rows": REG, "first_w": 3.6,
         "note": ("Pre-specified analyses were recorded in the project log before femoral neck BMD values were examined. "
                  "The clinical consequence is quantified with the T-score −2.5 threshold, which is published and reproducible; FRAX-based eligibility was "
                  "not modelled because the FRAX equations are not published.")}

sp = pd.DataFrame(R["specs"]).sort_values("any").reset_index(drop=True)
rows = [["Site", "T reference", "Z reference", "Weighting", "Hormone users", "n", "T abnormal, Z within range: % (95% CI)", "T osteoporosis, Z within range: % (95% CI)"]]
for _, r in sp.iterrows():
    rows.append([r.site + (" (primary)" if r.primary else ""), r.tref, r.zref, r.weighting, r.hormone, f"{r.n:,}",
                 f"{r['any']:.1f} ({r.any_lo:.1f} to {r.any_hi:.1f})", f"{r.op:.2f} ({r.op_lo:.2f} to {r.op_hi:.2f})"])
AN = R3["spec_anova"]
TABS2 = {"kind": "table", "caption": "Table S2. Specification curve: all analytic combinations, women aged 45–55 years", "rows": rows, "first_w": 1.5,
         "widths": [1.5, 1.05, 1.15, 0.75, 0.85, 0.5, 1.2, 1.2],
         "note": (f"{len(sp)} combinations, ordered as in Figure S2; the primary analysis lies at the {R3['spec_primary_pct_rank']:.0f}th percentile "
                  f"({round(R3['spec_primary_pct_rank'] * len(sp) / 100)} of {len(sp)} combinations give a lower estimate). "
                  f"Share of variance across combinations (type II sums of squares): site {100 * AN['site']:.1f}%, T reference {100 * AN['tref']:.1f}%, "
                  f"Z reference {100 * AN['zref']:.1f}%, weighting {100 * AN['weighting']:.1f}%, hormone users {100 * AN['hormone']:.1f}%, residual {100 * AN['Residual']:.1f}%. "
                  f"NHANES III cut-points apply to femoral neck and total hip; the 2005–2010 reference uses the weighted mean and SD of non-Hispanic White women "
                  f"aged 20–29 years in NHANES 2005–2010 (femoral neck {TREF['FN']['mean']:.3f}/{TREF['FN']['sd']:.3f}; total hip {TREF['TH']['mean']:.3f}/"
                  f"{TREF['TH']['sd']:.3f} g/cm²). {LSREF} For multi-site combinations the spine always uses its 2005–2010 reference. "
                  "Unweighted rows use equal weights with design-based CIs. The design is not fully crossed (lumbar spine has a single T reference), so the "
                  "variance shares for site and T reference are approximate.")}

G = [("T normal", "Concordant normal"), ("T abnormal, Z within range", "Discordant (T abnormal, Z within range)"), ("Z below expected range", "Z below expected range")]
rows = [["Outcome and group", "n", "Mean age, y", "Crude: % (95% CI)", "Age-standardized: % (95% CI)", "Age-standardized difference vs concordant normal: points (95% CI)"]]
for o, lab in (("frac", "Prior hip, wrist or spine fracture"), ("told", "Ever told by a doctor of osteoporosis")):
    rows.append([lab, "", "", "", "", ""])
    for g, gl in G:
        c = R1["outcome"][o][g]; s_ = R1["outcome_std"][o].get(g)
        std = f"{s_['p']:.1f} ({s_['lo']:.1f} to {s_['hi']:.1f})" if s_ else "—"
        dif = (f"{s_['d']:.1f} ({s_['dlo']:.1f} to {s_['dhi']:.1f})" if g != "T normal" else "Reference") if s_ else "—"
        rows.append([f"  {gl}", f"{c['n']:,}", f"{R1['outcome_age'][g]['m']:.1f}", ci(c), std, dif])
rows.append(["Prevalence ratio, discordant vs concordant normal", "", "", "", "", ""])
for o, lab in (("frac", "Prior fracture"), ("told", "Told of osteoporosis")):
    for k, v in R["pr"][o].items():
        rows.append([f"  {lab}: {k[0].lower() + k[1:] if k != 'Crude' else 'crude'}", f"{v['n']:,}", "", "", "", f"PR {v['pr']:.2f} ({v['lo']:.2f} to {v['hi']:.2f})"])
TABS3 = {"kind": "table", "caption": "Table S4. Prior fracture and diagnosed osteoporosis by concordance group, women aged 45–55 years", "rows": rows, "first_w": 2.9,
         "note": ("Self-reported, doctor-diagnosed fracture at any age. n excludes missing responses ("
                  + "; ".join(f"{(lambda s: s if s.startswith('Z ') else s.lower())(gl.split(' (')[0])}: {R1['outcome_age'][g]['n'] - R1['outcome']['frac'][g]['n']} missing fracture, "
                              f"{R1['outcome_age'][g]['n'] - R1['outcome']['told'][g]['n']} missing osteoporosis item" for g, gl in G)
                  + "). Age-standardized by the direct method to all women aged 45–55 years "
                  f"(45–49 y, {R1['std_shares'][0] * 100:.1f}%; 50–55 y, {R1['std_shares'][1] * 100:.1f}%); the Z-below group (n = 18) is shown crude only. "
                  "Prevalence ratios from survey-weighted Poisson regression; the fully adjusted model (age, body mass index, race and ethnicity, female hormone use) "
                  "excludes women with missing covariates. These comparisons contrast women with lower and normal BMD; the comparison of the two rules themselves is "
                  "in Figure 3. A prior osteoporosis diagnosis may itself reflect an earlier DXA, so that outcome is not independent of BMD. "
                  f"{DESIGN} PR, prevalence ratio.")}

SH4, LT = R4["shift"], R4["late"]
_P = pd.read_csv(OUT / "women45_55_v3.csv")
n_indet_disc = int(((_P.tcat >= 1) & (_P.Z > -2) & _P.status.isin(["Hysterectomy, ovaries conserved", "No periods, other reason", "Not ascertained"])).sum())
rows = [["Scenario", "Unweighted n with a discordant category", "US women aged 45–55 y with a discordant category, millions (95% CI)"],
        ["Structurally indeterminate status, 2005–2018", "", ""],
        ["  Upper bound: all such women", f"{n_indet_disc}", mil(CN["45–55|upper"])],
        ["  Expected: probability of still being premenopausal at own age (primary)", f"{n_indet_disc}", mil(CN["45–55|expected"])],
        ["  Expected: probability taken 1 year older (earlier ovarian decline)", f"{n_indet_disc}", mil(SH4["shift1"])],
        ["  Expected: probability taken 2 years older", f"{n_indet_disc}", mil(SH4["shift2"])],
        ["Exploratory, 2005–2006 cycle only (relative SE > 30%; direction only)", "", ""],
        [f"  Late menopausal transition: still menstruating, last period 3 months ago or more (share of menstruating women {LT['share_late']['p']:.1f}%, "
         f"95% CI {LT['share_late']['lo']:.1f} to {LT['share_late']['hi']:.1f}; n = {LT['n_late']})", f"{LT['n_late_disc']}", mil(LT["count_late_disc"]) + "*"],
        ["  For comparison: structurally indeterminate status, same cycle", f"{LT['n_indet_disc_D']}", mil(LT["count_indet_disc_D"]) + "*"]]
TABS4 = {"kind": "table", "caption": "Table S3. Sensitivity of the number of women for whom the choice of reporting rule is open, women aged 45–55 years",
         "rows": rows, "first_w": 4.6,
         "note": ("Ovaries conserved at hysterectomy may decline earlier than intact ovaries, which would lower the expected count; shifting the probability "
                  "curve 1 or 2 years older models this. Women still menstruating but in the late menopausal transition may be reported under either rule "
                  "(ISCD 2023), which would raise the count; time since the last period was asked only in 2005–2006 (n = "
                  f"{LT['n_menstruating']} menstruating women aged 45–55 years, {LT['n_late']} with a last period 3 months ago or more), so this row uses the "
                  f"2005–2006 examination weights alone (weighted total {LT['N_45_55_D']:.1f} million women aged 45–55 years) and is not added to the pooled estimate. "
                  "Three months is the closest questionnaire boundary to the 60 days of amenorrhea that define the late transition. 95% CIs: 1,000 Rao–Wu "
                  "bootstrap replicates for the expected counts; Taylor linearization for 2005–2006 rows. "
                  f"*Estimate may be unreliable: relative standard error {LT['rse_late_disc']:.0f}% and {LT['rse_indet_disc_D']:.0f}%, above the 30% threshold "
                  "commonly applied to NHANES estimates; these rows indicate only that the primary expected count is conservative.")}

LM = R3["ls_missing_compare"]
FIGS1 = {"kind": "figure", "caption": "Figure S1. Selection of participants", "path": OUT / "Figure_S1.png",
         "note": ("NHANES femur and spine DXA were offered to participants aged 8 years or older in 2005–2010, 40 years or older in 2013–2014 and 50 years or "
                  "older in 2017–2018. Lumbar spine results are missing mainly because one or more vertebrae were invalid; women without a valid spine result were "
                  f"older (mean {LM['age_missing']['m']:.1f} vs {LM['age_valid']['m']:.1f} years) and had higher body mass index ({LM['bmi_missing']['m']:.1f} vs "
                  f"{LM['bmi_valid']['m']:.1f} kg/m²), and femoral neck discordance was {LM['disc_missing']['p']:.1f}% vs {LM['disc_valid']['p']:.1f}%.")}
FIGS2 = {"kind": "figure", "caption": "Figure S2. Specification curve for the proportion with a discordant category, women aged 45–55 years", "path": OUT / "Figure_S2.png",
         "note": (f"(A) Estimates and 95% CIs for all {len(sp)} combinations, ordered by size and coloured by site; circled, pre-specified primary analysis "
                  f"({R3['spec_primary_pct_rank']:.0f}th percentile). (B) Settings of each combination. (C) Share of variance across combinations explained by each "
                  "setting (type II sums of squares). Values in Table S2.")}
FIGS3 = {"kind": "figure", "caption": "Figure S3. Discordance by participant characteristics, women aged 45–55 years", "path": OUT / "Figure_S3.png",
         "note": ("Discordant: T-score category low bone mass or osteoporosis with a Z-score within the expected range; dashed line, overall estimate. "
                  "Millions: weighted proportion of all women aged 45–55 years in that row with a discordant category, multiplied by the weighted total. "
                  "T-scores use a single white reference for all women whereas Z-scores are race-specific, so discordance differs by race and ethnicity by "
                  "construction; women with higher body mass index have higher BMD and fall less often in the discordant band. Status groups are in Figure 2D. "
                  "Unit: participant.")}
FIGS4 = {"kind": "figure", "caption": "Figure S4. Femoral neck T-scores and Z-scores by Z-score reference group, women aged 45–55 years", "path": OUT / "Figure_S4.png",
         "note": (f"n = {X['primary']['n']:,}. Each race and 5-year age band has its own reference mean and SD, so women in the same band lie on one straight line; "
                  "the several parallel lines within one colour are the 45–49, 50–54 and 55–59 year bands. Unit: participant.")}
FIGS5 = {"kind": "figure", "caption": "Figure S5. Discordance by single year of age compared with the age spline, women aged 40–59 years", "path": OUT / "Figure_S5.png",
         "note": ("Points, weighted single-year estimates (95% CI); line, restricted cubic spline from Figure 2B. Single-year estimates above 55 years rest on fewer "
                  "participants per year and fluctuate accordingly. Unit: participant.")}

ITEMS = [TAB1, TAB2, FIG1, FIG2, FIG3]
SUPP = [TABS1, TABS2, TABS4, TABS3,  # 2026-10-28 依正文引用順序對調 S3／S4
        FIGS1, FIGS2, FIGS3, FIGS4, FIGS5]
PNGS = ["Figure_1", "Figure_2", "Figure_3", "Figure_S1", "Figure_S2", "Figure_S3", "Figure_S4", "Figure_S5"]


def source_data(path):
    P40 = pd.read_csv(OUT / "women40_59_v2.csv")
    sheets = []
    for it in ITEMS + SUPP:
        if it["kind"] == "table":
            sheets.append((it["caption"].split(".")[0].replace(" ", "_"), numeric_table(it["rows"]), it["caption"]))
    P = pd.read_csv(OUT / "women45_55_v3.csv")
    curve = lambda d, pre: pd.DataFrame({"age": d["x"], f"{pre}_pct": d["p"], f"{pre}_lo": d["lo"], f"{pre}_hi": d["hi"]})
    hrow = lambda sec, g, v: {"section": sec, "group": g, "n": v["any"]["n"], "discordant_pct": v["any"]["p"], "lo": v["any"]["lo"], "hi": v["any"]["hi"],
                              "osteoporosis_flip_pct": v["op"]["p"], "op_lo": v["op"]["lo"], "op_hi": v["op"]["hi"],
                              "millions": v["mil"]["t"], "mil_lo": v["mil"]["lo"], "mil_hi": v["mil"]["hi"], "share_pct": v["share"]["p"]}
    H = R3["hetero"]; CN = R3["counts"]
    sheets += [
        ("Figure_1A", pd.DataFrame([{"age": 50, "race": "Non-Hispanic White", "femoral_neck_BMD": 0.70,
                                     "T_reference_mean": 0.86, "T_reference_SD": 0.12, "T_score": (0.70 - 0.86) / 0.12,
                                     "Z_reference_mean_age50_54": R["zmean"]["3|50"], "Z_reference_SD_age50_54": R["zsd"]["3|50"],
                                     "Z_score": (0.70 - R["zmean"]["3|50"]) / R["zsd"]["3|50"]}]), "Figure 1A (worked example)"),
        ("Figure_1C", pd.DataFrame([{"T_category": c, "Z_category": z, "millions": R3["xtab"]["FN"][f"{c}|{z}"]["mil"]["t"],
                                     "lo": R3["xtab"]["FN"][f"{c}|{z}"]["mil"]["lo"], "hi": R3["xtab"]["FN"][f"{c}|{z}"]["mil"]["hi"],
                                     "pct": R3["xtab"]["FN"][f"{c}|{z}"]["p"]}
                                    for c in ("Normal", "Low bone mass", "Osteoporosis") for z in ("within", "below")]), "Figure 1C (Sankey flows)"),
        ("Figure_1B", P40[["SEQN", "RIDAGEYR", "RIDRETH1", "DXXNKBMD", "w"]].rename(columns={"RIDAGEYR": "age", "RIDRETH1": "race_code", "DXXNKBMD": "femoral_neck_BMD", "w": "pooled_weight"}), "Figure 1B (participants)"),
        ("Figure_1B_Zcut", pd.DataFrame([{"race_code": k.split("|")[0], "age_band_start": k.split("|")[1], "Z_minus2_BMD": v} for k, v in R["zcut"].items()]), "Figure 1B Z −2.0 boundaries"),
        ("Figure_1D", pd.DataFrame([{"site": n_, "millions": R3["xtab"][k]["op"]["mil"]["t"], "lo": R3["xtab"][k]["op"]["mil"]["lo"], "hi": R3["xtab"][k]["op"]["mil"]["hi"]}
                                    for k, n_ in (("FN", "Femoral neck"), ("LOW3", "Lowest of three sites"))]), "Figure 1D"),
        ("Figure_2A", pd.DataFrame(R3["composition"]), "Figure 2A (percent)"),
        ("Figure_2B", curve(R3["spline_disc"], "discordant"), "Figure 2B"),
        ("Figure_2C", pd.DataFrame([{"age_band": f"{a}–{a + 1}", "layer": l, **CN[f"{a}–{a + 1}|{l}"]} for a in range(40, 60, 2) for l in ("upper", "expected")]), "Figure 2C (millions)"),
        ("Figure_2C_p", curve(R3["p_premenopausal"], "premenopausal"), "Figure 2C probability of still being premenopausal"),
        ("Figure_2D", pd.DataFrame([hrow("Menstrual and surgical status", g, v) for g, v in H["Menstrual and surgical status"].items()]), "Figure 2D"),
        ("Figure_3A", pd.DataFrame({"T_score": R["frac_curve"]["T"], "fracture_pct": R["frac_curve"]["p"], "lo": R["frac_curve"]["lo"], "hi": R["frac_curve"]["hi"]}), "Figure 3A"),
        ("Figure_3B", pd.DataFrame([{"comparison": k, **v} for k, v in R3["pr_tz"].items()]), "Figure 3B"),
        ("Figure_3C", pd.DataFrame([{"score": k, **v} for k, v in R4["auc"].items()]), "Figure 3C"),
        ("Figure_S1", pd.DataFrame(R1["flow"] + [["Not in primary age range", R3["flow_extra"]["excluded_age"]], ["Valid lumbar spine, 45–55 y", R3["flow_extra"]["ls_valid"]]]
                                   + [[f"No valid spine: {k}", v] for k, v in R3["flow_extra"]["ls_missing"].items()], columns=["step", "n"]), "Figure S1"),
        ("Figure_S2", pd.DataFrame(R["specs"]), "Figure S2A–B"),
        ("Figure_S2C", pd.DataFrame([{"setting": k, "share_of_variance": v} for k, v in R3["spec_anova"].items()]), "Figure S2C"),
        ("Figure_S3", pd.DataFrame([hrow(s_, g, v) for s_, d in H.items() if s_ != "Menstrual and surgical status" for g, v in d.items()]), "Figure S3"),
        ("Figure_S4", P[["SEQN", "RIDRETH1", "T", "Z"]].rename(columns={"RIDRETH1": "race_code"}), "Figure S4"),
        ("Figure_S5", pd.DataFrame([{"age": int(a), "pct": v["any"]["p"], "lo": v["any"]["lo"], "hi": v["any"]["hi"], "n": v["any"]["n"]} for a, v in R1["age1"].items()]), "Figure S5 (points)"),
    ]
    wb = Workbook(); toc = wb.active; toc.title = "Contents"
    toc.append(["Sheet", "Exhibit"]); toc["A1"].font = toc["B1"].font = XFont(bold=True)
    for name, df, cap in sheets:
        ws = wb.create_sheet(name[:31]); toc.append([name[:31], cap])
        ws.append(list(df.columns))
        for c in ws[1]: c.font = XFont(bold=True)
        for r in df.itertuples(index=False):
            ws.append([None if (isinstance(v, float) and v != v) else (v.item() if hasattr(v, "item") else v) for v in r])
    wb.save(path)


def main():
    os.makedirs(DEST, exist_ok=True)
    doc = os.path.join(DEST, "Tables_and_Figures.docx")
    nt, nf = build(ITEMS, doc, SUPP)
    for f in PNGS:
        shutil.copy(OUT / f"{f}.png", os.path.join(DEST, f"{f}.png"))
    sd = os.path.join(DEST, "Source_data.xlsx"); source_data(sd)
    print(f"{nt} tables, {nf} figures ->", DEST)


if __name__ == "__main__":
    main()
