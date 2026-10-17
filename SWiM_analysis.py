# -*- coding: utf-8 -*-
"""SWiM minimal submission analysis.

Outputs only:
- Figure4_effect_direction.png
- Table3_direction_counts_with_exactP_and_WilsonCI.csv
- Table9_RoB_summary.csv
"""

from pathlib import Path
from math import sqrt

import pandas as pd
from scipy.stats import binomtest
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
INPUT = HERE / "SWiM_data.xlsx"
OUT = HERE / "SWiM_reanalysis_output"
OUT.mkdir(exist_ok=True)


# Remove the obsolete output from the previous version, if present.
legacy_summary = OUT / "SWiM_summary_counts.csv"
if legacy_summary.exists():
    legacy_summary.unlink()


def wilson_ci(k, n, z=1.959963984540054):
    p = k / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * sqrt(
        p * (1 - p) / n + z * z / (4 * n * n)
    ) / den
    return center - half, center + half


def time_code(s):
    s = str(s)
    for w in ("W1", "W2", "W3", "W4"):
        if s.startswith(w):
            return w
    return "TTE"


def strategy_abbrev(strategy):
    return {
        "Zoledronic acid": "ZOL",
        "Bisphosphonate, ibandronate": "IBN",
        "IV BP, mostly ZOL": "IV-BP",
        "Oral BP": "Oral BP",
        "Romosozumab": "ROMO",
        "Teriparatide": "TPTD",
        "SERM": "SERM",
        "No treatment / discontinuation": "No Tx",
        "Bridge combination": "Bridge",
        "Complex sequence": "Complex",
    }.get(strategy, strategy)


# ============================================================
# Read only the two retained sheets
# ============================================================

study = pd.read_excel(
    INPUT,
    sheet_name="StudyLevel_All_175"
)

rob = pd.read_excel(
    INPUT,
    sheet_name="Report_RoB_43"
)

bmd = study[
    study["Synthesis_Class"].eq("Primary BMD trajectory")
].copy()


# ============================================================
# 1) Table 3
# Exact two-sided sign test + Wilson 95% CI
# ============================================================

rows = []

for (strategy, outcome, tw), g in bmd.groupby(
    ["Strategy_Family", "Standard_Outcome", "Time_Window"],
    sort=False
):
    gb = g[
        g["Study_Binary_Code"].isin([1, -1])
    ]

    if len(gb) >= 2:
        n = len(gb)
        fav = int(
            (gb["Study_Binary_Code"] == 1).sum()
        )
        unf = n - fav

        p = binomtest(
            fav,
            n,
            p=0.5,
            alternative="two-sided"
        ).pvalue

        lo, hi = wilson_ci(fav, n)

        rows.append([
            strategy,
            outcome,
            tw,
            n,
            fav,
            unf,
            p,
            lo,
            hi,
            "; ".join(
                gb["Study_Unit_ID"].astype(str)
            ),
        ])


t8 = pd.DataFrame(
    rows,
    columns=[
        "Strategy_Family",
        "Outcome",
        "Time_Window",
        "Binary_n",
        "Favorable",
        "Unfavorable",
        "Exact_two_sided_P",
        "Wilson95_low",
        "Wilson95_high",
        "Study_IDs",
    ]
)


# ============================================================
# 2) Table 9
# Report-level risk-of-bias summary
# ============================================================

def rs(name, sub):

    c = lambda tool, overall: int(
        (
            (sub["Tool"] == tool)
            & (sub["Overall"] == overall)
        ).sum()
    )

    return [
        name,
        sub["Census_ID"].nunique(),
        sub["Study_Unit_ID"].nunique(),
        c("RoB 2", "Low risk"),
        c("RoB 2", "Some concerns"),
        c("RoB 2", "High risk"),
        c("ROBINS-I", "Moderate"),
        c("ROBINS-I", "Serious"),
    ]


rr = [
    rs("SWiM all", rob)
]

for name, token in [
    ("BMD", "BMD"),
    ("Fracture", "Fracture"),
    ("BTM/TBS", "BTM/TBS"),
    ("Structural", "Structural"),
]:
    rr.append(
        rs(
            name,
            rob[
                rob["Included_outcomes"]
                .astype(str)
                .str.contains(
                    token,
                    regex=False,
                    na=False
                )
            ],
        )
    )


rdf = pd.DataFrame(
    rr,
    columns=[
        "Evidence",
        "Reports",
        "Studies",
        "RoB2_Low",
        "RoB2_SomeConcerns",
        "RoB2_High",
        "ROBINSI_Moderate",
        "ROBINSI_Serious",
    ]
)


# ============================================================
# 3) Figure 4
# Effect-direction plot
# ============================================================

site_x = {
    "LS-BMD": 0,
    "TH-BMD": 1,
    "FN-BMD": 2
}

groups = {}
order = []
reports_by_key = {}

for _, r in bmd.iterrows():

    key = (
        r["SWiM_Main_Section"],
        r["Strategy_Family"],
        r["Study_Unit_ID"],
        r["Time_Window"],
    )

    if key not in groups:
        groups[key] = {}
        order.append(key)
        reports_by_key[key] = []

    for report in str(r["Reports"]).split(";"):
        report = report.strip()
        if report and report not in reports_by_key[key]:
            reports_by_key[key].append(report)

    groups[key][
        r["Standard_Outcome"]
    ] = r["Study_Level_Direction"]


labels = [
    (
        f"{strategy_abbrev(k[1])} | "
        f"{'; '.join(reports_by_key[k])} | "
        f"{time_code(k[3])}"
    )
    for k in order
]


fig, ax = plt.subplots(
    figsize=(10.23, 18.26),
    dpi=100
)


for direction, marker, label in [
    (
        "Favorable trajectory",
        "^",
        "Favorable"
    ),
    (
        "Unfavorable trajectory",
        "v",
        "Unfavorable"
    ),
    (
        "Mixed/no clear trajectory",
        "D",
        "Mixed / no clear direction"
    ),
]:

    xs = []
    ys = []

    for yi, k in enumerate(order):

        for site, d in groups[k].items():

            if (
                d == direction
                and site in site_x
            ):
                xs.append(
                    site_x[site]
                )
                ys.append(
                    yi
                )

    ax.scatter(
        xs,
        ys,
        marker=marker,
        s=55,
        label=label
    )


ax.set_xticks(
    [0, 1, 2],
    [
        "LS-BMD",
        "TH-BMD",
        "FN-BMD"
    ]
)

ax.set_yticks(
    range(len(labels)),
    labels,
    fontsize=8.5
)

ax.invert_yaxis()

ax.set_xlim(
    -0.45,
    2.45
)

ax.set_title(
    "Effect-direction plot of post-denosumab BMD trajectories",
    fontsize=15
)

ax.set_xlabel(
    "BMD site"
)

ax.set_ylabel(
    "Strategy | Reports | post-DMAb time window"
)

ax.legend(
    loc="upper right",
    frameon=False
)


prev = None

for yi, k in enumerate(order):

    if (
        prev is not None
        and k[0] != prev
    ):
        ax.axhline(
            yi - 0.5,
            linewidth=0.8
        )

    prev = k[0]


# ============================================================
# 4) Hard checks
# Stop if the retained submission data have changed
# ============================================================

# Study-level BMD data
assert len(bmd) == 103, (
    f"Expected 103 primary BMD rows, "
    f"found {len(bmd)}"
)

assert (
    bmd["Study_Binary_Code"]
    .isin([1, -1])
    .sum()
    == 99
), "Expected 99 binary BMD directions"

assert (
    (bmd["Study_Binary_Code"] == 1).sum()
    == 41
), "Expected 41 favorable BMD directions"

assert (
    (bmd["Study_Binary_Code"] == -1).sum()
    == 58
), "Expected 58 unfavorable BMD directions"

assert (
    (~bmd["Study_Binary_Code"].isin([1, -1])).sum()
    == 4
), "Expected 4 mixed/no-clear BMD directions"


# RoB data
assert (
    rob["Census_ID"].nunique()
    == 43
), "Expected 43 RoB reports"

assert (
    rob["Study_Unit_ID"].nunique()
    == 24
), "Expected 24 RoB study units"


# Table 3
assert len(t8) == 17, (
    f"Expected 17 Table 3 rows, "
    f"found {len(t8)}"
)


# Table 9
assert len(rdf) == 5, (
    f"Expected 5 Table 9 rows, "
    f"found {len(rdf)}"
)


# Figure 4
assert bmd["Reports"].notna().all(), (
    "Expected non-missing Reports for all primary BMD rows"
)

assert bmd["Reports"].astype(str).str.strip().ne("").all(), (
    "Expected non-empty Reports for all primary BMD rows"
)

assert all(reports_by_key[k] for k in order), (
    "Expected at least one report label for every Figure 4 row"
)

assert len(order) == 46, (
    f"Expected 46 Figure 4 rows, "
    f"found {len(order)}"
)

assert (
    sum(
        len(v)
        for v in groups.values()
    )
    == 103
), "Expected 103 Figure 4 site-direction points"


# ============================================================
# 5) Write the three retained outputs
# Only reached after all hard checks pass
# ============================================================

t8.to_csv(
    OUT / "Table3_direction_counts_with_exactP_and_WilsonCI.csv",
    index=False,
    encoding="utf-8-sig"
)

rdf.to_csv(
    OUT / "Table9_RoB_summary.csv",
    index=False,
    encoding="utf-8-sig"
)

fig.tight_layout()

fig.savefig(
    OUT / "Figure4_effect_direction.png",
    dpi=160
)

plt.close(fig)


print(
    "Verification complete. "
    "Results were written to:",
    OUT
)
