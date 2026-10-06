
"""
================================================================================
 Meta-analysis of sequential treatment after denosumab
================================================================================
 Usage: modify INPUT_FILE and OUTPUT_DIR below, then run
        python run_meta_analysis.py
 Dependencies: pandas / numpy / scipy / matplotlib / openpyxl

 Outputs:
   OUTPUT_DIR/01_Cleaned_Dataset.xlsx      complete cleaning and conversion record
   OUTPUT_DIR/02_Meta_Analysis_Results.xlsx    primary, heterogeneity, sensitivity, and bias results
   OUTPUT_DIR/figures/*.png / *.svg     publication-ready static figures and animated SVG figures
================================================================================
"""
import os, re, math, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")


INPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meta_data.xlsx")

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meta_analysis_output")

FIG_DIR = os.path.join(OUTPUT_DIR, "figures")

FIG_DIR    = os.path.join(OUTPUT_DIR, "figures")


SHEET_MAP = {
    "DMAb-ZOL": ("Denosumab→Zoledronate",   "DMAb→ZOL",  "Zoledronate"),
    "DMAb-ROMO": ("Denosumab→Romosozumab", "DMAb→ROMO", "Romosozumab"),
    "DMAb-ALN": ("Denosumab→Alendronate", "DMAb→ALN",  "Alendronate"),
    "DMAb-RAL": ("Denosumab→Raloxifene",   "DMAb→RAL",  "Raloxifene"),
}
SITES = [("LS", "Lumbar spine"), ("FN", "Femoral neck"), ("TH", "Total hip")]
Z975  = 1.959963984540054

DISPERSION_AS_SD = ("Ha 2022",)

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)





#   SEM → SD :  SD = SEM × √n

#   SD  → SE :  SE = SD / √n


def _to_pct(x):
    """Parse a mean cell as a percentage; supports fractional numeric values and percent strings."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    if isinstance(x, str):
        s = x.strip().replace("−", "-").replace("％", "%")
        if s in ("", "/", "\\", "-", "—"):
            return None
        if "%" in s:
            return float(re.sub(r"[^0-9.\-+]", "", s))
        try:
            v = float(s)
        except ValueError:
            return None
        return v * 100.0 if abs(v) < 1 else v
    v = float(x)

    return v * 100.0 if abs(v) < 1 else v


def _parse_ci(s):
    """Parse a confidence-interval string into (lower, upper)."""
    t = (s.replace("−", "-").replace(" to ", ",")
          .replace("％", "%").replace("（", "(").replace("）", ")"))
    t = re.sub(r"\(\s*95\s*%?\s*C[IiLl]\s*\)|95\s*%?\s*C[IiLl]", "", t)
    t = re.sub(r"[()\[\]%\s]", "", t)
    if "," in t:
        parts = [p for p in t.split(",") if p not in ("", "-")]
    else:

        m = re.search(r"(?<=\d)-", t)
        parts = [t[:m.start()], t[m.end():]] if m else [t]
    vals = []
    for p in parts:
        p = re.sub(r"[^0-9.\-+]", "", p)
        if p not in ("", "-", "+"):
            vals.append(float(p))
    if len(vals) != 2:
        return None
    return (min(vals), max(vals))


def _parse_dispersion(raw, n, z=Z975):
    """Parse a dispersion cell into (SD%, SE%, original type, conversion note)."""
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return (None,) * 4
    if isinstance(raw, str):
        s = raw.strip()
        if s in ("", "/", "\\", "-", "—"):
            return (None,) * 4
        if re.search(r"C[IlL]|95", s):
            ci = _parse_ci(s)
            if ci is None:
                return (None,) * 4
            lo, hi = ci
            se = (hi - lo) / (2 * z)
            sd = se * math.sqrt(n)
            return sd, se, "95%CI", f"SE=({hi:.4g}−{lo:.4g})/(2×{z:.3f})={se:.4f}; SD=SE×√{n}={sd:.4f}"
        if "SD" in s.upper():
            sd = float(re.sub(r"[^0-9.\-+]", "", s.replace("−", "-")))
            se = sd / math.sqrt(n)
            return sd, se, "SD", f"Reported SD={sd:.4f}; SE=SD/√{n}={se:.4f}"
        try:
            v = float(s.replace("−", "-"))
        except ValueError:
            return (None,) * 4
        raw = v
    sem = float(raw)
    sem = sem * 100.0 if abs(sem) < 1 else sem
    sd = sem * math.sqrt(n)
    return sd, sem, "SEM", f"Reported SEM={sem:.4f}; SD=SEM×√{n}={sd:.4f}"


def load_and_clean(path, ci_crit="z", sd_override=DISPERSION_AS_SD):
    """Read all relevant XLSX worksheets into long-format analysis-unit data."""
    xl = pd.ExcelFile(path)
    rows = []
    for sheet in xl.sheet_names:
        if sheet not in SHEET_MAP:
            continue
        cn, en, drug = SHEET_MAP[sheet]
        raw = pd.read_excel(xl, sheet, header=None)

        hdr_idx = next(i for i in range(len(raw))
                       if raw.iloc[i].astype(str).str.contains("LS sample size").any())
        hdr = {str(v).strip(): j for j, v in enumerate(raw.iloc[hdr_idx]) if pd.notna(v)}
        for i in range(hdr_idx + 1, len(raw)):
            study = raw.iloc[i, 0]
            if pd.isna(study) or str(study).strip() == "":
                continue
            study = re.sub(r"（.*?）", "", str(study)).strip()
            for site, site_cn in SITES:
                cn_n, cn_m, cn_d = f"{site} sample size", f"{site} mean", f"{site}-SEM"
                if cn_n not in hdr or cn_m not in hdr:
                    continue
                n_raw = raw.iloc[i, hdr[cn_n]]
                mean  = _to_pct(raw.iloc[i, hdr[cn_m]])
                if pd.isna(n_raw) or str(n_raw).strip() in ("/", "\\", "") or mean is None:
                    continue
                n = int(float(n_raw))
                d_raw = raw.iloc[i, hdr[cn_d]] if cn_d in hdr else None

                if (not isinstance(d_raw, str)) and pd.notna(d_raw) and study.startswith(sd_override):
                    sd = float(d_raw)
                    sd = sd * 100.0 if abs(sd) < 1 else sd
                    sd_, se_, typ = sd, sd / math.sqrt(n), "SD (source verified)"
                    note = f"Unlabeled in source table; verified as SD={sd:.4f}; SE=SD/√{n}={sd_ / math.sqrt(n):.4f}"
                else:

                    zc = Z975 if ci_crit == "z" else stats.t.ppf(0.975, max(n - 1, 1))
                    sd_, se_, typ, note = _parse_dispersion(d_raw, n, zc)
                if sd_ is None:
                    continue
                rows.append(dict(plan=cn, plan_en=en, subsequent_drug=drug, study=study,
                                 site=site, site_name=site_cn, n=n, mean_change=mean,
                                 SD=sd_, SE=se_, dispersion_type=typ,
                                 raw_dispersion=str(d_raw), conversion_note=note,
                                 CI_lower=mean - Z975 * se_, CI_upper=mean + Z975 * se_))
    df = pd.DataFrame(rows)
    order = {k: i for i, k in enumerate([v[0] for v in SHEET_MAP.values()])}
    df["_o"] = df["plan"].map(order)
    df = df.sort_values(["_o", "site", "study"]).drop(columns="_o").reset_index(drop=True)
    return df







def _tau2_dl(y, v):
    w = 1.0 / v
    theta_fe = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - theta_fe) ** 2)
    k = len(y)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0, Q


def _tau2_reml(y, v, tol=1e-10, itmax=200):
    t2 = _tau2_dl(y, v)[0]
    for _ in range(itmax):
        w = 1.0 / (v + t2)
        mu = np.sum(w * y) / np.sum(w)

        num = np.sum(w ** 2 * ((y - mu) ** 2 - v + 1.0 / np.sum(w)))
        den = np.sum(w ** 2)
        new = max(0.0, num / den)
        if abs(new - t2) < tol:
            t2 = new
            break
        t2 = new
    return t2


def _tau2_pm(y, v, tol=1e-10, itmax=200):
    """Paule-Mandel moment estimator: solve sum[w(y-mu)^2] = k-1."""
    k = len(y)
    lo, hi = 0.0, max(1e-8, np.var(y) * 10 + np.max(v) * 10)
    f = lambda t2: np.sum((y - np.sum(y / (v + t2)) / np.sum(1 / (v + t2))) ** 2 / (v + t2)) - (k - 1)
    if f(0.0) <= 0:
        return 0.0
    for _ in range(itmax):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2


def meta_analyze(y, se, n=None, model="DL", hksj=False, alpha=0.05):
    """
    y: study effects (BMD change, %)
    se: study standard errors (%)
    model: 'DL' | 'REML' | 'PM' | 'FE'
    hksj: apply Hartung-Knapp-Sidik-Jonkman variance adjustment
    Returns pooled effect, CI, heterogeneity statistics, and prediction interval.
    """
    y = np.asarray(y, float); se = np.asarray(se, float)
    v = se ** 2; k = len(y)
    if k == 0:
        return None
    tau2_dl, Q = _tau2_dl(y, v) if k > 1 else (0.0, 0.0)
    if model == "FE":
        tau2 = 0.0
    elif model == "REML":
        tau2 = _tau2_reml(y, v) if k > 1 else 0.0
    elif model == "PM":
        tau2 = _tau2_pm(y, v) if k > 1 else 0.0
    else:
        tau2 = tau2_dl
    w = 1.0 / (v + tau2)
    theta = np.sum(w * y) / np.sum(w)
    se_theta = math.sqrt(1.0 / np.sum(w))
    if hksj and k > 1:
        se_theta = math.sqrt(np.sum(w * (y - theta) ** 2) / ((k - 1) * np.sum(w)))
        crit = stats.t.ppf(1 - alpha / 2, k - 1)
        pval = 2 * (1 - stats.t.cdf(abs(theta / se_theta), k - 1))
    else:
        crit = stats.norm.ppf(1 - alpha / 2)
        pval = 2 * (1 - stats.norm.cdf(abs(theta / se_theta)))
    df_q = max(k - 1, 1)
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if (k > 1 and Q > 0) else 0.0
    p_q = 1 - stats.chi2.cdf(Q, df_q) if k > 1 else np.nan

    if k >= 3:
        tcrit = stats.t.ppf(1 - alpha / 2, k - 2)
        pi_half = tcrit * math.sqrt(tau2 + se_theta ** 2)
        pi = (theta - pi_half, theta + pi_half)
    else:
        pi = (np.nan, np.nan)
    return dict(k=k, N=int(np.sum(n)) if n is not None else np.nan,
                theta=theta, se=se_theta, lo=theta - crit * se_theta,
                hi=theta + crit * se_theta, p=pval, Q=Q, df=df_q, p_Q=p_q,
                I2=I2, tau2=tau2, tau=math.sqrt(tau2), H2=Q / df_q if k > 1 else 1.0,
                pi_lo=pi[0], pi_hi=pi[1], model=model, hksj=hksj,
                w_pct=(w / np.sum(w) * 100))


def egger_test(y, se):
    """Egger linear regression test (exploratory): standardized normal deviate versus precision."""
    k = len(y)
    if k < 3:
        return dict(k=k, intercept=np.nan, se=np.nan, t=np.nan, p=np.nan)
    snd = y / se; prec = 1.0 / se
    res = stats.linregress(prec, snd)
    tval = res.intercept / res.intercept_stderr
    p = 2 * (1 - stats.t.cdf(abs(tval), k - 2))
    return dict(k=k, intercept=res.intercept, se=res.intercept_stderr, t=tval, p=p)




def run_primary(df, model="DL", hksj=False):
    """Pool by sequential regimen and measurement site; with k=1, report without pooling."""
    out = []
    for site, site_cn in SITES:
        for plan in [v[0] for v in SHEET_MAP.values()]:
            s = df[(df.plan == plan) & (df.site == site)]
            if len(s) == 0:
                continue
            r = meta_analyze(s.mean_change.values, s.SE.values, s.n.values,
                             model=model, hksj=hksj)
            out.append(dict(site=site, site_name=site_cn, plan=plan,
                            plan_en=s.plan_en.iloc[0], k_studies=r["k"], N_total=r["N"],
                            pooled_change=r["theta"], CI_lower=r["lo"], CI_upper=r["hi"],
                            SE=r["se"], P_value=r["p"], Q=r["Q"], Q_df=r["df"],
                            Q_p=r["p_Q"], I2=r["I2"], tau2=r["tau2"], tau=r["tau"],
                            prediction_lower=r["pi_lo"], prediction_upper=r["pi_hi"]))
    return pd.DataFrame(out)


def subgroup_test(res, site):
    """Between-regimen test (Q_between) using pooled effects and SEs."""
    s = res[(res.site == site) & (res.k_studies >= 2)]
    if len(s) < 2:
        return None
    w = 1.0 / s.SE.values ** 2
    theta = np.sum(w * s.pooled_change.values) / np.sum(w)
    Q = np.sum(w * (s.pooled_change.values - theta) ** 2)
    dfree = len(s) - 1
    return dict(site=site, plans_compared=len(s), Q_between=Q, df=dfree,
                P_value=1 - stats.chi2.cdf(Q, dfree),
                I2_between=max(0.0, (Q - dfree) / Q * 100) if Q > 0 else 0.0)




def leave_one_out(df):
    """S1 leave-one-out sensitivity analysis for strata with k>=3."""
    out = []
    for (plan, site), s in df.groupby(["plan", "site"], sort=False):
        if len(s) < 3:
            continue
        full = meta_analyze(s.mean_change.values, s.SE.values, s.n.values)
        for i in range(len(s)):
            sub = s.drop(s.index[i])
            r = meta_analyze(sub.mean_change.values, sub.SE.values, sub.n.values)
            out.append(dict(plan=plan, site=site, omitted_study=s.study.iloc[i], remaining_k=r["k"],
                            pooled_change=r["theta"], CI_lower=r["lo"], CI_upper=r["hi"],
                            I2=r["I2"], P_value=r["p"],
                            difference_from_primary=r["theta"] - full["theta"],
                            primary_pooled_change=full["theta"],
                            conclusion_changed="Yes" if (r["p"] < 0.05) != (full["p"] < 0.05) else "No"))
    return pd.DataFrame(out)


def _combine_two_arms(s, keys):
    """Combine two subgroup arms from the same study into one analysis unit."""
    a, b = s[s.study == keys[0]].iloc[0], s[s.study == keys[1]].iloc[0]
    n1, n2 = a.n, b.n
    m1, m2 = a.mean_change, b.mean_change
    s1, s2 = a.SD, b.SD
    n = n1 + n2
    m = (n1 * m1 + n2 * m2) / n
    sd = math.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2 + n1 * n2 / n * (m1 - m2) ** 2) / (n - 1))
    row = a.copy()
    row.study, row.n, row.mean_change, row.SD, row.SE = "Haeri 2025 (Combined)", n, m, sd, sd / math.sqrt(n)
    return row


def run_sensitivity(path, df_main):
    """S2-S10: rerun all regimen-by-site strata under nine scenarios."""
    scen = {}

    scen["S2 Fixed-effect model"]        = run_primary(df_main, model="FE")
    scen["S3 REML estimate of τ²"]         = run_primary(df_main, model="REML")
    scen["S4 Paule-Mandel estimate of τ²"] = run_primary(df_main, model="PM")
    scen["S5 HKSJ variance adjustment"]       = run_primary(df_main, hksj=True)


    scen["S6 Native SEM/SD only"] = run_primary(df_main[df_main.dispersion_type.str.startswith(("SEM", "SD"))])


    scen["S7 t critical value for CI conversion"] = run_primary(load_and_clean(path, ci_crit="t"))


    scen["S8 Interpret Ha 2022 as SEM"] = run_primary(load_and_clean(path, sd_override=()))


    d9 = []
    for (plan, site), s in df_main.groupby(["plan", "site"], sort=False):
        keys = [x for x in ["Haeri 2025 (Female)", "Haeri 2025 (Male)"] if x in set(s.study)]
        if len(keys) == 2:
            s = pd.concat([s[~s.study.isin(keys)],
                           pd.DataFrame([_combine_two_arms(s, keys)])], ignore_index=True)
        d9.append(s)
    scen["S9 Combine Haeri female/male"] = run_primary(pd.concat(d9, ignore_index=True))


    scen["S10 Exclude n<20"] = run_primary(df_main[df_main.n >= 20])

    rows = []
    for name, r in scen.items():
        for _, x in r.iterrows():
            base = df_main.pipe(lambda d: None)
            rows.append(dict(scenario=name, site=x.site, plan=x.plan, k=x.k_studies,
                             pooled_change=x.pooled_change, CI_lower=x.CI_lower, CI_upper=x.CI_upper,
                             I2=x.I2, P_value=x.P_value))
    sens = pd.DataFrame(rows)
    main = run_primary(df_main)[["site", "plan", "k_studies", "pooled_change", "CI_lower", "CI_upper", "I2", "P_value"]]
    main.columns = ["site", "plan", "k", "pooled_change", "CI_lower", "CI_upper", "I2", "P_value"]
    main.insert(0, "scenario", "S0 Primary analysis (DL random effects)")
    sens = pd.concat([main, sens], ignore_index=True)
    m = main.set_index(["site", "plan"])["pooled_change"]
    sens["difference_from_primary"] = sens.apply(lambda r: r.pooled_change - m.get((r.site, r.plan), np.nan), axis=1)
    sens["direction_consistent"] = sens.apply(
        lambda r: "—" if pd.isna(r.difference_from_primary) else
        ("Yes" if np.sign(r.pooled_change) == np.sign(m.get((r.site, r.plan), np.nan)) else "No"), axis=1)
    return sens




def run_bias(df):
    """Egger test; when k<10, report as exploratory only."""
    out = []
    for (plan, site), s in df.groupby(["plan", "site"], sort=False):
        if len(s) < 3:
            continue
        e = egger_test(s.mean_change.values, s.SE.values)
        out.append(dict(plan=plan, site=site, k=e["k"], Egger_intercept=e["intercept"],
                        intercept_SE=e["se"], t_value=e["t"], P_value=e["p"],
                        interpretation="k<10, exploratory only" if e["k"] < 10 else "Interpretable"))
    for site, site_cn in SITES:
        s = df[df.site == site]
        if len(s) >= 3:
            e = egger_test(s.mean_change.values, s.SE.values)
            out.append(dict(plan="All regimens combined (exploratory)", site=site, k=e["k"],
                            Egger_intercept=e["intercept"], intercept_SE=e["se"], t_value=e["t"],
                            P_value=e["p"], interpretation="Populations differ across regimens; funnel-plot symmetry reference only"))
    return pd.DataFrame(out)




def export_clean_workbook(df, path_out):
    """Export cleaned dataset with conversion notes and Excel verification formulas."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook(); ws = wb.active; ws.title = "Cleaned dataset"
    head = ["No.", "Sequential regimen", "plan_en", "subsequent_drug", "Study (analysis unit)", "site", "site_name",
            " sample sizen", "12-month BMD change (%)", "Raw dispersion entry", "Original dispersion type",
            "Converted SD (%)", "Converted SE (%)", "SD verification formula", "SE verification formula",
            "95%CI_lower(%)", "95%CI_upper(%)", "Conversion rationale"]
    ws.append(head)
    for i, r in df.reset_index(drop=True).iterrows():
        row = i + 2
        ws.append([i + 1, r.plan, r.plan_en, r.subsequent_drug, r.study, r.site, r.site_name,
                   r.n, round(r.mean_change, 4), r.raw_dispersion, r.dispersion_type,
                   round(r.SD, 4), round(r.SE, 4), None, None,
                   round(r.CI_lower, 4), round(r.CI_upper, 4), r.conversion_note])
        ws.cell(row, 14).value = f"=M{row}"
        ws.cell(row, 15).value = f"=L{row}/SQRT(H{row})"
        ws.cell(row, 14).value = f"=M{row}*SQRT(H{row})"   # SD = SE×√n
    hdr_fill = PatternFill("solid", fgColor="1F4E79")
    thin = Side(style="thin", color="BFBFBF")
    for c in range(1, len(head) + 1):
        cell = ws.cell(1, c); cell.font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
        cell.fill = hdr_fill; cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, max_col=len(head)):
        for cell in row:
            cell.font = Font(name="Arial", size=10)
            cell.border = Border(left=thin, right=thin, top=thin, bottom=thin)
            cell.alignment = Alignment(vertical="center")
    for c, w in zip(range(1, len(head) + 1),
                    [6, 20, 12, 13, 24, 7, 9, 9, 18, 26, 13, 13, 16, 16, 13, 13, 46]):
        ws.column_dimensions[get_column_letter(c)].width = w
    ws.freeze_panes = "E2"


    ws2 = wb.create_sheet("Conversion rules and notes")
    notes = [
        ["Item", "Description"],
        ["Effect measure", "Percentage BMD change from month 0 to month 12 after subsequent treatment following denosumab discontinuation (single-arm mean change, %)"],
        ["Default dispersion", "Entries without a specific label are interpreted as SEM, as stated in the source table."],
        ["SEM → SD", "SD = SEM × √n"],
        ["95%CI → SE", "SE = (upper - lower) / (2 x 1.96); then SD = SE x sqrt(n) (Cochrane Handbook 6.5.2.3)."],
        ["SD → SE", "SE = SD / √n"],
        ["Units", "All values are standardized to percentages (%); fractional entries such as -0.048 are multiplied by 100."],
        ["Analysis unit", "Independent subgroups from the same study (e.g., Haeri 2025 female/male) are included as separate analysis units and combined in sensitivity analysis S9."],
        ["Missing data", "Entries marked slash or backslash are treated as not reported for that site; no imputation is performed."],
        ["Source verification", "Ha 2022 was unlabeled in the source table but verified from the original report as SD; the primary analysis treats it as SD, while S8 checks the SEM interpretation."],
    ]
    for r in notes:
        ws2.append(r)
    for c in range(1, 3):
        ws2.cell(1, c).font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
        ws2.cell(1, c).fill = hdr_fill
    ws2.column_dimensions["A"].width = 16; ws2.column_dimensions["B"].width = 105
    for row in ws2.iter_rows(min_row=2, max_row=ws2.max_row, max_col=2):
        for cell in row:
            cell.font = Font(name="Arial", size=10)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
    wb.save(path_out)


def export_result_workbook(tables, path_out):
    """Write result tables to one formatted workbook."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter
    wb = Workbook(); wb.remove(wb.active)
    hdr_fill = PatternFill("solid", fgColor="1F4E79")
    thin = Side(style="thin", color="BFBFBF")
    for name, (df_t, note) in tables.items():
        ws = wb.create_sheet(name[:31])
        ws.append([note]); ws.append([])
        ws.append(list(df_t.columns))
        for _, r in df_t.iterrows():
            ws.append([round(v, 4) if isinstance(v, (int, float, np.floating)) and not pd.isna(v)
                       else ("" if pd.isna(v) else v) for v in r.tolist()])
        ws.cell(1, 1).font = Font(name="Arial", bold=True, size=11, color="1F4E79")
        for c in range(1, len(df_t.columns) + 1):
            cell = ws.cell(3, c)
            cell.font = Font(name="Arial", bold=True, color="FFFFFF", size=10)
            cell.fill = hdr_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            ws.column_dimensions[get_column_letter(c)].width = max(
                11, min(30, int(df_t.iloc[:, c - 1].astype(str).str.len().max() if len(df_t) else 11) + 4,
                        ))
        for row in ws.iter_rows(min_row=4, max_row=ws.max_row, max_col=len(df_t.columns)):
            for cell in row:
                cell.font = Font(name="Arial", size=10)
                cell.border = Border(left=thin, right=thin, top=thin, bottom=thin)
                cell.alignment = Alignment(vertical="center", horizontal="center")
        ws.freeze_panes = "A4"
    wb.save(path_out)


# ================================ 7. Create figures ====================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 9

PLAN_COLOR = {"Denosumab→Zoledronate": "#2C6E9B", "Denosumab→Romosozumab": "#2E8B57",
              "Denosumab→Alendronate": "#7B5AA6", "Denosumab→Raloxifene": "#C0504D"}
GRID, INK, MUTE = "#D8DEE4", "#1A1A1A", "#6B7280"
SVG_FONT = "'Noto Sans CJK SC','Microsoft YaHei',Arial,sans-serif"


def _fmt(v, d=2):
    return "—" if (v is None or (isinstance(v, float) and (np.isnan(v)))) else f"{v:+.{d}f}"


def build_forest_spec(df, res, site, site_cn):
    """Build forest-plot rows for each regimen, studies, and pooled estimate."""
    rows = []
    for plan in [v[0] for v in SHEET_MAP.values()]:
        s = df[(df.plan == plan) & (df.site == site)]
        if len(s) == 0:
            continue
        r = res[(res.plan == plan) & (res.site == site)]
        if len(r) == 0:
            continue
        r = r.iloc[0]
        rows.append(dict(kind="header", label=f"{plan}（{r.plan_en}）", color=PLAN_COLOR[plan]))
        m = meta_analyze(s.mean_change.values, s.SE.values, s.n.values)
        for i, (_, x) in enumerate(s.iterrows()):
            rows.append(dict(kind="study", label=f"　{x.study}", n=int(x.n),
                             est=x.mean_change, lo=x.CI_lower, hi=x.CI_upper,
                             w=m["w_pct"][i], color=PLAN_COLOR[plan]))
        het = (f"I²={r.I2:.1f}%，τ²={r.tau2:.2f}，P={r.Q_p:.3f}"
               if r.k_studies > 1 else "One study only; not pooled")
        rows.append(dict(kind="summary", label=f"  Pooled effect (k={int(r.k_studies)}，N={int(r.N_total)}）",
                         est=r.pooled_change, lo=r.CI_lower, hi=r.CI_upper, p=r.P_value,
                         het=het, color=PLAN_COLOR[plan]))
        rows.append(dict(kind="blank"))
    if rows and rows[-1]["kind"] == "blank":
        rows.pop()
    return rows


def draw_forest_png(rows, title, subtitle, out_png, xlabel="12-month BMD change after discontinuation (%)", show_title=False):
    n = len(rows)
    fig_h = max(3.2, 0.42 * n + 1.9)
    fig, ax = plt.subplots(figsize=(11.2, fig_h))
    vals = [r for r in rows if r["kind"] in ("study", "summary")]
    lo = min(min(r["lo"] for r in vals), 0) - 1.2
    hi = max(max(r["hi"] for r in vals), 0) + 1.2
    ax.set_xlim(lo, hi); ax.set_ylim(-0.5, n - 0.5); ax.invert_yaxis()
    ax.axvline(0, color=MUTE, lw=1.0, ls="--", zorder=1)
    ax.set_xlabel(xlabel, fontsize=10)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_yticks([]); ax.tick_params(axis="x", colors=INK, labelsize=9)
    ax.grid(axis="x", color=GRID, lw=0.6, alpha=0.6, zorder=0)

    xt = lo - (hi - lo) * 0.55
    xr = hi + (hi - lo) * 0.05
    xw = hi + (hi - lo) * 0.40
    ax.text(xt, -1.15, "study / Analysis unit", ha="left", va="center", fontweight="bold", fontsize=9.5, clip_on=False)
    ax.text(xr, -1.15, "Change (95% CI)", ha="left", va="center", fontweight="bold", fontsize=9.5, clip_on=False)
    ax.text(xw, -1.15, "Weight", ha="left", va="center", fontweight="bold", fontsize=9.5, clip_on=False)

    for i, r in enumerate(rows):
        if r["kind"] == "header":
            ax.text(xt, i, r["label"], ha="left", va="center", fontweight="bold",
                    fontsize=9.5, color=r["color"], clip_on=False)
        elif r["kind"] == "study":
            ax.plot([r["lo"], r["hi"]], [i, i], color=r["color"], lw=1.4, solid_capstyle="round", zorder=3)
            for e in (r["lo"], r["hi"]):
                ax.plot([e, e], [i - 0.13, i + 0.13], color=r["color"], lw=1.4, zorder=3)
            sz = 40 + 220 * (r["w"] / 100.0)
            ax.scatter([r["est"]], [i], s=sz, marker="s", color=r["color"], zorder=4, edgecolors="white", lw=0.6)
            ax.text(xt, i, f'{r["label"]}（n={r["n"]}）', ha="left", va="center", fontsize=9, clip_on=False)
            ax.text(xr, i, f'{_fmt(r["est"])} ({_fmt(r["lo"])}, {_fmt(r["hi"])})', ha="left",
                    va="center", fontsize=8.8, clip_on=False)
            ax.text(xw, i, f'{r["w"]:.1f}%', ha="left", va="center", fontsize=8.8, color=MUTE, clip_on=False)
        elif r["kind"] == "summary":
            c, e, l, h = r["color"], r["est"], r["lo"], r["hi"]
            ax.add_patch(Polygon([[l, i], [e, i - 0.26], [h, i], [e, i + 0.26]],
                                 closed=True, facecolor=c, edgecolor=c, alpha=0.95, zorder=5))
            ax.text(xt, i, r["label"], ha="left", va="center", fontsize=9, fontweight="bold", clip_on=False)
            star = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
            ax.text(xr, i, f'{_fmt(e)} ({_fmt(l)}, {_fmt(h)}){star}', ha="left", va="center",
                    fontsize=9, fontweight="bold", clip_on=False)
            ax.text(xw, i, "100%", ha="left", va="center", fontsize=8.8, color=MUTE, clip_on=False)
            ax.text(xt, i + 0.52, "　　" + r["het"], ha="left", va="center", fontsize=7.8,
                    color=MUTE, style="italic", clip_on=False)
    if show_title:
        fig.suptitle(title, x=0.02, ha="left", fontsize=13, fontweight="bold", y=0.985)
        fig.text(0.02, 0.945, subtitle, ha="left", fontsize=8.8, color=MUTE)
    fig.subplots_adjust(left=0.30, right=0.72, top=0.90 if show_title else 0.97, bottom=0.13)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _nice_ticks(lo, hi, target=8):
    """Generate clean ticks using 1/2/2.5/5 step sizes."""
    raw = (hi - lo) / target
    mag = 10 ** math.floor(math.log10(raw))
    step = min([m for m in (1, 2, 2.5, 5, 10) if m * mag >= raw] or [10]) * mag
    t0 = math.ceil(lo / step) * step
    return [round(t0 + i * step, 6) for i in range(int((hi - t0) / step) + 1)]


def write_forest_svg(rows, title, subtitle, out_svg, xlabel="12-month BMD change after discontinuation (%)"):
    """Write SVG with SMIL entry animation."""
    vals = [r for r in rows if r["kind"] in ("study", "summary")]
    lo = min(min(r["lo"] for r in vals), 0) - 1.2
    hi = max(max(r["hi"] for r in vals), 0) + 1.2
    W, RH, TOP = 1180, 26, 92
    PX0, PX1 = 372, 852
    H = TOP + RH * len(rows) + 74
    sx = lambda v: PX0 + (v - lo) / (hi - lo) * (PX1 - PX0)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="{SVG_FONT}">',
         f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>',
         f'<text x="24" y="34" font-size="19" font-weight="700" fill="{INK}">{title}</text>',
         f'<text x="24" y="55" font-size="12" fill="{MUTE}">{subtitle}</text>',
         f'<text x="24" y="{TOP-12}" font-size="12.5" font-weight="700" fill="{INK}">study / Analysis unit</text>',
         f'<text x="{PX1+30}" y="{TOP-12}" font-size="12.5" font-weight="700" fill="{INK}">Change (95% CI)</text>',
         f'<text x="{W-70}" y="{TOP-12}" font-size="12.5" font-weight="700" fill="{INK}">Weight</text>',
         f'<line x1="24" y1="{TOP-5}" x2="{W-24}" y2="{TOP-5}" stroke="{GRID}" stroke-width="1"/>']

    for t in _nice_ticks(lo, hi):
        o.append(f'<line x1="{sx(t):.1f}" y1="{TOP-2}" x2="{sx(t):.1f}" y2="{TOP+RH*len(rows)+6}" '
                 f'stroke="{GRID}" stroke-width="0.8"/>')
        o.append(f'<text x="{sx(t):.1f}" y="{TOP+RH*len(rows)+24}" font-size="11" fill="{MUTE}" '
                 f'text-anchor="middle">{t:+.6g}</text>')
    o.append(f'<line x1="{sx(0):.1f}" y1="{TOP-2}" x2="{sx(0):.1f}" y2="{TOP+RH*len(rows)+6}" '
             f'stroke="{MUTE}" stroke-width="1.2" stroke-dasharray="5 4"/>')
    o.append(f'<text x="{(PX0+PX1)/2:.0f}" y="{TOP+RH*len(rows)+48}" font-size="12" fill="{INK}" '
             f'text-anchor="middle">{xlabel}</text>')

    for i, r in enumerate(rows):
        y = TOP + RH * i + RH / 2
        d = 0.10 * i
        if r["kind"] == "header":
            o.append(f'<text x="24" y="{y+4:.0f}" font-size="13.5" font-weight="700" '
                     f'fill="{r["color"]}" opacity="1">{r["label"]}'
                     f'<animate attributeName="opacity" from="0" to="1" begin="{d:.2f}s" dur="0.4s" fill="freeze"/>'
                     f'</text>')
        elif r["kind"] == "study":
            c = r["color"]; x0, x1, xe = sx(r["lo"]), sx(r["hi"]), sx(r["est"])
            side = 7 + 9 * (r["w"] / 100.0)
            o.append(f'<g><title>{r["label"].strip()}: {_fmt(r["est"])}% ({_fmt(r["lo"])}, {_fmt(r["hi"])})</title>'
                     f'<line x1="{x0:.1f}" y1="{y:.0f}" x2="{x1:.1f}" y2="{y:.0f}" stroke="{c}" stroke-width="2" '
                     f'stroke-linecap="round">'
                     f'<animate attributeName="x1" from="{xe:.1f}" to="{x0:.1f}" begin="{d:.2f}s" dur="0.65s" fill="freeze"/>'
                     f'<animate attributeName="x2" from="{xe:.1f}" to="{x1:.1f}" begin="{d:.2f}s" dur="0.65s" fill="freeze"/>'
                     f'</line>')
            for xx in (x0, x1):
                o.append(f'<line x1="{xx:.1f}" y1="{y-5:.0f}" x2="{xx:.1f}" y2="{y+5:.0f}" stroke="{c}" '
                         f'stroke-width="2" opacity="1">'
                         f'<animate attributeName="opacity" from="0" to="1" begin="{d+0.55:.2f}s" dur="0.3s" fill="freeze"/></line>')
            o.append(f'<rect x="{xe-side/2:.1f}" y="{y-side/2:.1f}" width="{side:.1f}" height="{side:.1f}" '
                     f'fill="{c}" stroke="#FFFFFF" stroke-width="1.2" opacity="1">'
                     f'<animate attributeName="opacity" from="0" to="1" begin="{d:.2f}s" dur="0.35s" fill="freeze"/></rect></g>')
            o.append(f'<text x="24" y="{y+4:.0f}" font-size="12.5" fill="{INK}">{r["label"]}（n={r["n"]}）</text>')
            o.append(f'<text x="{PX1+30}" y="{y+4:.0f}" font-size="12" fill="{INK}">'
                     f'{_fmt(r["est"])} ({_fmt(r["lo"])}, {_fmt(r["hi"])})</text>')
            o.append(f'<text x="{W-70}" y="{y+4:.0f}" font-size="12" fill="{MUTE}">{r["w"]:.1f}%</text>')
        elif r["kind"] == "summary":
            c = r["color"]; x0, x1, xe = sx(r["lo"]), sx(r["hi"]), sx(r["est"])
            pts = f'{x0:.1f},{y:.0f} {xe:.1f},{y-8:.0f} {x1:.1f},{y:.0f} {xe:.1f},{y+8:.0f}'
            o.append(f'<g transform="translate({xe:.1f},{y:.0f})">'
                     f'<g transform="scale(1,1)">'
                     f'<polygon points="{x0-xe:.1f},0 0,-8 {x1-xe:.1f},0 0,8" fill="{c}" stroke="{c}"/>'
                     f'<animateTransform attributeName="transform" type="scale" from="0 0" to="1 1" '
                     f'begin="{d+0.2:.2f}s" dur="0.5s" fill="freeze" calcMode="spline" '
                     f'keySplines="0.2 0.9 0.3 1" keyTimes="0;1"/></g></g>')
            star = "***" if r["p"] < 0.001 else "**" if r["p"] < 0.01 else "*" if r["p"] < 0.05 else ""
            o.append(f'<text x="24" y="{y+4:.0f}" font-size="12.5" font-weight="700" fill="{INK}">{r["label"]}</text>')
            o.append(f'<text x="{PX1+30}" y="{y+4:.0f}" font-size="12.5" font-weight="700" fill="{INK}">'
                     f'{_fmt(r["est"])} ({_fmt(r["lo"])}, {_fmt(r["hi"])}){star}</text>')
            o.append(f'<text x="{W-70}" y="{y+4:.0f}" font-size="12" fill="{MUTE}">100%</text>')
            o.append(f'<text x="40" y="{y+18:.0f}" font-size="10.5" font-style="italic" fill="{MUTE}">{r["het"]}</text>')
    o.append(f'<text x="24" y="{H-10}" font-size="10.5" fill="{MUTE}">'
             f'Note: square area is proportional to random-effects weight; diamonds show pooled effects and 95% CIs; *P&lt;0.05, **P&lt;0.01, ***P&lt;0.001.</text>')
    o.append("</svg>")
    open(out_svg, "w", encoding="utf-8").write("\n".join(o))


def save_fig_both(fig, base, titles):
    """Keep titles in SVG; omit titles from PNG for manuscript captions."""
    fig.savefig(base + ".svg", facecolor="white")
    for a in titles:
        a.remove()
    fig.savefig(base + ".png", dpi=300, facecolor="white")
    plt.close(fig)
    animate_mpl_svg(base + ".svg")


def animate_mpl_svg(path, stagger=0.045, dur=0.65):
    """Inject staggered CSS entry animation into matplotlib SVG."""
    css = ["<style><![CDATA[",
           "@keyframes mfade{from{opacity:0;transform:translateY(8px)}to{opacity:1;transform:translateY(0)}}",
           'g[id^="line2d"],g[id^="PathCollection"],g[id^="patch"],g[id^="LineCollection"]'
           "{animation:mfade %.2fs ease-out both;transform-box:fill-box;transform-origin:center}" % dur]
    for i in range(1, 61):
        css.append('g[id^="line2d"]:nth-of-type(%d),g[id^="PathCollection"]:nth-of-type(%d),'
                   'g[id^="patch"]:nth-of-type(%d){animation-delay:%.2fs}' % (i, i, i, i * stagger))
    css.append("]]></style>")
    s = open(path, encoding="utf-8").read()
    i = s.find(">", s.find("<svg")) + 1
    open(path, "w", encoding="utf-8").write(s[:i] + "\n".join(css) + s[i:])


def fig_summary_matrix(res, out_base):
    """Figure 4: pooled effects by regimen and measurement site."""
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4), sharex=True)
    plans = [v[0] for v in SHEET_MAP.values()]
    for ax, (site, site_cn) in zip(axes, SITES):
        sub = res[res.site == site]
        yy, seen = [], []
        for j, plan in enumerate(plans):
            r = sub[sub.plan == plan]
            if len(r) == 0:
                continue
            r = r.iloc[0]; c = PLAN_COLOR[plan]; y = len(seen)
            ax.plot([r.CI_lower, r.CI_upper], [y, y], color=c, lw=2.6, solid_capstyle="round", zorder=3)
            ax.scatter([r.pooled_change], [y], s=110, color=c, zorder=4, edgecolors="white", lw=1.2,
                       marker="D" if r.k_studies > 1 else "o")
            lbl = f"{r.plan_en}\n(k={int(r.k_studies)}, N={int(r.N_total)})"
            seen.append(lbl); yy.append(y)
            ax.text(r.pooled_change, y - 0.30, f"{r.pooled_change:+.2f}%", ha="center", fontsize=9,
                    fontweight="bold", color=c)
        ax.axvline(0, color=MUTE, ls="--", lw=1)
        ax.set_yticks(yy); ax.set_yticklabels(seen, fontsize=9)
        ax.set_title(site_cn + f"（{site}）", fontsize=11.5, fontweight="bold")
        ax.set_ylim(-0.8, len(seen) - 0.2); ax.invert_yaxis()
        ax.grid(axis="x", color=GRID, lw=0.6)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
    gl = min(res.CI_lower.min(), 0) - 1.2; gh = max(res.CI_upper.max(), 0) + 1.2
    for ax in axes:
        ax.set_xlim(gl, gh)
    axes[1].set_xlabel("12-month BMD change after denosumab discontinuation (%, random-effects pooled estimate and 95% CI)", fontsize=10)
    _t = [fig.suptitle("Figure 4  Pooled effects across sequential regimens and measurement sites", x=0.01, ha="left",
                       fontsize=13, fontweight="bold")]
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig_both(fig, out_base, _t)


def fig_loo(loo, res, site, out_base):
    """Figure 5: leave-one-out sensitivity analysis."""
    d = loo[loo.site == site]
    plans = [p for p in [v[0] for v in SHEET_MAP.values()] if p in set(d.plan)]
    fig, axes = plt.subplots(1, len(plans), figsize=(4.4 * len(plans), 3.9), squeeze=False)
    for ax, plan in zip(axes[0], plans):
        s = d[d.plan == plan].reset_index(drop=True); c = PLAN_COLOR[plan]
        main = res[(res.plan == plan) & (res.site == site)].iloc[0]
        ax.axvspan(main.CI_lower, main.CI_upper, color=c, alpha=0.12, zorder=0)
        ax.axvline(main.pooled_change, color=c, ls="--", lw=1.4, zorder=1)
        ax.axvline(0, color=MUTE, ls=":", lw=1)
        for i, r in s.iterrows():
            ax.plot([r.CI_lower, r.CI_upper], [i, i], color=c, lw=2, solid_capstyle="round")
            ax.scatter([r.pooled_change], [i], s=55, color=c, zorder=3, edgecolors="white", lw=0.9)
        ax.set_yticks(range(len(s)))
        ax.set_yticklabels(["Omit " + x.replace(" (injection after 6 months)", "") for x in s.omitted_study], fontsize=8.5)
        ax.invert_yaxis(); ax.set_title(plan, fontsize=10.5, fontweight="bold", color=c)
        ax.grid(axis="x", color=GRID, lw=0.5)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
        ax.set_xlabel("pooled_change（%）", fontsize=9)
    _t = [fig.suptitle(f"Figure 5  Leave-one-out sensitivity analysis ({dict(SITES)[site]})  Shaded band indicates the primary-analysis 95% CI",
                 x=0.01, ha="left", fontsize=12.5, fontweight="bold")]
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig_both(fig, out_base, _t)


def fig_scenarios(sens, site, out_base):
    """Figure 6: scenario sensitivity analysis summary."""
    d = sens[sens.site == site]
    scen = list(dict.fromkeys(d.scenario))
    plans = [p for p in [v[0] for v in SHEET_MAP.values()] if p in set(d.plan)]
    fig, ax = plt.subplots(figsize=(10.6, 0.42 * len(scen) + 2.4))
    off = np.linspace(-0.26, 0.26, len(plans))
    for j, plan in enumerate(plans):
        c = PLAN_COLOR[plan]
        for i, sc in enumerate(scen):
            r = d[(d.scenario == sc) & (d.plan == plan)]
            if len(r) == 0:
                continue
            r = r.iloc[0]; y = i + off[j]
            ax.plot([r.CI_lower, r.CI_upper], [y, y], color=c, lw=2, solid_capstyle="round")
            ax.scatter([r.pooled_change], [y], s=42, color=c, zorder=3, edgecolors="white", lw=0.8,
                       label=plan if i == 0 else None)
    ax.axvline(0, color=MUTE, ls="--", lw=1)
    ax.set_yticks(range(len(scen))); ax.set_yticklabels(scen, fontsize=9)
    ax.invert_yaxis(); ax.grid(axis="x", color=GRID, lw=0.5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("Pooled change (%, with 95% CI)", fontsize=10)
    ax.legend(frameon=False, fontsize=9, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.30))
    _t = [fig.suptitle(f"Figure 6  Scenario sensitivity analysis ({dict(SITES)[site]}）", x=0.01, ha="left",
                 fontsize=12.5, fontweight="bold")]
    fig.tight_layout(rect=[0, 0.03, 1, 0.94])
    save_fig_both(fig, out_base, _t)


def fig_funnel(df, res, site, out_base, egger_row=None):
    """Figure 7: exploratory funnel plot with inverted SE axis and pseudo-95% limits."""
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    d = df[df.site == site]
    semax = d.SE.max() * 1.15
    plans = [p for p in [v[0] for v in SHEET_MAP.values()] if p in set(d.plan)]
    for plan in plans:
        s = d[d.plan == plan]
        r = res[(res.plan == plan) & (res.site == site)].iloc[0]
        ax.scatter(s.mean_change, s.SE, s=62, color=PLAN_COLOR[plan], alpha=0.9,
                   edgecolors="white", lw=1, label=f"{plan}（k={len(s)}）", zorder=3)
        ax.plot([r.pooled_change, r.pooled_change], [0, semax], color=PLAN_COLOR[plan], lw=1.2, ls="--", alpha=0.8)
        se = np.linspace(0.001, semax, 60)
        ax.plot(r.pooled_change - Z975 * se, se, color=PLAN_COLOR[plan], lw=0.9, alpha=0.45)
        ax.plot(r.pooled_change + Z975 * se, se, color=PLAN_COLOR[plan], lw=0.9, alpha=0.45)
    ax.set_ylim(semax, 0); ax.set_xlabel("BMD change (%)", fontsize=10)
    ax.set_ylabel("Standard error, SE (%)", fontsize=10)
    ax.grid(color=GRID, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8.8, loc="lower left")
    sub = "Pseudo-95% limits are centered on each regimen's pooled effect; k<10 in each stratum, so symmetry is visual/exploratory only."
    if egger_row is not None:
        sub += f" | Cross-regimen Egger intercept P={egger_row:.3f}"
    _t = [fig.suptitle(f"Figure 7  Funnel plot ({dict(SITES)[site]}, exploratory)", x=0.01, ha="left",
                 fontsize=12.5, fontweight="bold")]
    _t.append(fig.text(0.01, 0.925, sub, fontsize=8.6, color=MUTE))
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    save_fig_both(fig, out_base, _t)




def main():
    print("1. Read and clean data …")
    df = load_and_clean(INPUT_FILE)
    print(f"   Analysis unit {len(df)} / analysis arms {df.study.nunique()} / regimens {df.plan.nunique()}")

    print("2. Primary analysis (DerSimonian-Laird random effects)…")
    res = run_primary(df)

    print("3. Between-regimen difference test …")
    sg = pd.DataFrame([x for x in (subgroup_test(res, s) for s, _ in SITES) if x])
    sg["site_name"] = sg.site.map(dict(SITES))

    print("4. Sensitivity analyses (leave-one-out + 10 scenarios)…")
    loo = leave_one_out(df)
    sens = run_sensitivity(INPUT_FILE, df)

    print("5. Small-study effect test …")
    bias = run_bias(df)

    print("6. Export tables …")
    export_clean_workbook(df, os.path.join(OUTPUT_DIR, "01_Cleaned_Dataset.xlsx"))

    t1 = df.pivot_table(index=["plan", "study"], columns="site",
                        values=["n", "mean_change", "SD"], aggfunc="first")
    t1.columns = [f"{b}-{a}" for a, b in t1.columns]
    t1 = t1.reset_index()
    tables = {
        "Table1_Included_Studies": (t1, "Table 1  Characteristics of included analysis units and converted effect measures (%, converted SD)"),
        "Table2_Primary_Analysis": (res.round(4), "Table 2  Random-effects pooled results by sequential regimen and site (DerSimonian-Laird)"),
        "Table3_Between_Regimen_Test": (sg.round(4), "Table 3  Between-regimen difference test within each site (Q_between)"),
        "Table4_Leave_One_Out": (loo.round(4), "Table 4  Leave-one-out sensitivity analysis (strata with k>=3 only)"),
        "Table5_Scenario_Sensitivity": (sens.round(4), "Table 5  Scenario sensitivity analysis (S0 primary + S2-S11)"),
        "Table6_Small_Study_Effects": (bias.round(4), "Table 6  Egger test (k<10 in each stratum; exploratory only)"),
        "Table7_Cleaned_Unit_Data": (df.drop(columns=["plan_en", "subsequent_drug"]).round(4),
                                "Table 7  Cleaned and standardized analysis-unit data"),
    }
    export_result_workbook(tables, os.path.join(OUTPUT_DIR, "02_Meta_Analysis_Results.xlsx"))

    print("7. Create figures …")
    for idx, (site, site_cn) in enumerate(SITES, start=1):
        spec = build_forest_spec(df, res, site, site_cn)
        title = f"Figure {idx}  Forest plot of BMD change after sequential treatment following denosumab: {site_cn} BMD change"
        subtitle = "Effect: percentage BMD change from month 0 to month 12 after subsequent treatment following denosumab discontinuation; DerSimonian-Laird random-effects model"
        base = os.path.join(FIG_DIR, f"Fig{idx}_Forest_Plot_{site_cn}")
        draw_forest_png(spec, title, subtitle, base + ".png")
        write_forest_svg(spec, title, subtitle, base + ".svg")
    fig_summary_matrix(res, os.path.join(FIG_DIR, "Fig4_Regimen_Site_Pooled_Effects"))
    fig_loo(loo, res, "LS", os.path.join(FIG_DIR, "Fig5_Leave-one-out_Lumbar spine"))
    fig_scenarios(sens, "LS", os.path.join(FIG_DIR, "Fig6_Scenario sensitivity_Lumbar spine"))
    er = bias[(bias.plan.str.startswith("All")) & (bias.site == "LS")]
    fig_funnel(df, res, "LS", os.path.join(FIG_DIR, "Fig7_Funnel_Plot_Lumbar_Spine"),
               float(er.P_value.iloc[0]) if len(er) else None)

    print("Done. Output directory:", OUTPUT_DIR)
    return df, res, sg, loo, sens, bias


if __name__ == "__main__":
    main()
