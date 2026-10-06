# -*- coding: utf-8 -*-
"""SWiM final reanalysis after excluding two conference-abstract-only reports.
Put this script and SWiM_review_source_data_excluding_2_conference_abstracts_English.xlsx in the same folder,
then run:  python SWiM_reanalysis.py
"""
from pathlib import Path
from math import sqrt
import pandas as pd
import numpy as np
from scipy.stats import binomtest
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
INPUT = HERE / "SWiM_data.xlsx"
OUT = HERE / "SWiM_reanalysis_output"
OUT.mkdir(exist_ok=True)

def wilson_ci(k,n,z=1.959963984540054):
    p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half, center+half

def time_code(s):
    s=str(s)
    for w in ("W1","W2","W3","W4"):
        if s.startswith(w): return w
    return "TTE"

def strategy_abbrev(strategy):
    return {
        "Zoledronic acid":"ZOL",
        "Bisphosphonate, ibandronate":"IBN",
        "IV BP, mostly ZOL":"IV-BP",
        "Oral BP":"Oral BP",
        "Romosozumab":"ROMO",
        "Teriparatide":"TPTD",
        "SERM":"SERM",
        "No treatment / discontinuation":"No Tx",
        "Bridge combination":"Bridge",
        "Complex sequence":"Complex",
    }.get(strategy,strategy)

main=pd.read_excel(INPUT,sheet_name="MainEffect_Selected_434")
study=pd.read_excel(INPUT,sheet_name="StudyLevel_All_175")
rob=pd.read_excel(INPUT,sheet_name="Report_RoB_43")
bmd=study[study["Synthesis_Class"].eq("Primary BMD trajectory")].copy()

# 1) Overall counts
fract=main[main["Outcome_Domain"].eq("Fracture")]
btm=main[(main["Outcome_Domain"].eq("Bone turnover marker")) | (main["Standard_Outcome"].eq("TBS"))]
struct=rob[rob["Included_outcomes"].astype(str).str.contains("Structural",na=False)]
summary={
    "SWiM contributing reports":rob["Census_ID"].nunique(),
    "SWiM underlying study units (report/RoB census)":rob["Study_Unit_ID"].nunique(),
    "Primary BMD study-level rows":len(bmd),
    "Primary BMD binary directions":bmd["Study_Binary_Code"].isin([1,-1]).sum(),
    "Primary BMD favorable":(bmd["Study_Binary_Code"]==1).sum(),
    "Primary BMD unfavorable":(bmd["Study_Binary_Code"]==-1).sum(),
    "Primary BMD mixed/no clear":(~bmd["Study_Binary_Code"].isin([1,-1])).sum(),
    "Fracture underlying study units (outcome-level MainEffect)":fract["Study_Unit_ID"].nunique(),
    "BTM/TBS underlying study units (outcome-level MainEffect)":btm["Study_Unit_ID"].nunique(),
    "Structural/microstructural underlying study units (report/RoB census)":struct["Study_Unit_ID"].nunique(),
}
pd.DataFrame(list(summary.items()),columns=["Metric","Value"]).to_csv(OUT/"SWiM_summary_counts.csv",index=False,encoding="utf-8-sig")

# 2) Table 8: exact two-sided sign test and Wilson 95% CI
rows=[]
for (strategy,outcome,tw),g in bmd.groupby(["Strategy_Family","Standard_Outcome","Time_Window"],sort=False):
    gb=g[g["Study_Binary_Code"].isin([1,-1])]
    if len(gb)>=2:
        n=len(gb); fav=int((gb["Study_Binary_Code"]==1).sum()); unf=n-fav
        p=binomtest(fav,n,p=0.5,alternative="two-sided").pvalue
        lo,hi=wilson_ci(fav,n)
        rows.append([strategy,outcome,tw,n,fav,unf,p,lo,hi,"; ".join(gb["Study_Unit_ID"].astype(str))])
t8=pd.DataFrame(rows,columns=["Strategy_Family","Outcome","Time_Window","Binary_n","Favorable","Unfavorable","Exact_two_sided_P","Wilson95_low","Wilson95_high","Study_IDs"])
t8.to_csv(OUT/"Table8_direction_counts_with_exactP_and_WilsonCI.csv",index=False,encoding="utf-8-sig")

# 3) Table 9: report-level RoB summary
def rs(name,sub):
    c=lambda tool,overall:int(((sub["Tool"]==tool)&(sub["Overall"]==overall)).sum())
    return [name,sub["Census_ID"].nunique(),sub["Study_Unit_ID"].nunique(),c("RoB 2","Low risk"),c("RoB 2","Some concerns"),c("RoB 2","High risk"),c("ROBINS-I","Moderate"),c("ROBINS-I","Serious")]
rr=[rs("SWiM all",rob)]
for name,token in [("BMD","BMD"),("Fracture","Fracture"),("BTM/TBS","BTM/TBS"),("Structural","Structural")]:
    rr.append(rs(name,rob[rob["Included_outcomes"].astype(str).str.contains(token,regex=False,na=False)]))
rdf=pd.DataFrame(rr,columns=["Evidence","Reports","Studies","RoB2_Low","RoB2_SomeConcerns","RoB2_High","ROBINSI_Moderate","ROBINSI_Serious"])
rdf.to_csv(OUT/"Table9_RoB_summary.csv",index=False,encoding="utf-8-sig")

# 4) Figure S1 effect-direction plot
site_x={"LS-BMD":0,"TH-BMD":1,"FN-BMD":2}
groups={}
order=[]
for _,r in bmd.iterrows():
    key=(r["SWiM_Main_Section"],r["Strategy_Family"],r["Study_Unit_ID"],r["Time_Window"])
    if key not in groups:
        groups[key]={}; order.append(key)
    groups[key][r["Standard_Outcome"]]=r["Study_Level_Direction"]
labels=[f"{strategy_abbrev(k[1])} | {k[2]} | {time_code(k[3])}" for k in order]
fig,ax=plt.subplots(figsize=(10.23,18.26),dpi=100)
for direction,marker,label in [("Favorable trajectory","^","Favorable"),("Unfavorable trajectory","v","Unfavorable"),("Mixed/no clear trajectory","D","Mixed / no clear direction")]:
    xs=[]; ys=[]
    for yi,k in enumerate(order):
        for site,d in groups[k].items():
            if d==direction and site in site_x:
                xs.append(site_x[site]); ys.append(yi)
    ax.scatter(xs,ys,marker=marker,s=55,label=label)
ax.set_xticks([0,1,2],["LS-BMD","TH-BMD","FN-BMD"])
ax.set_yticks(range(len(labels)),labels,fontsize=8.5)
ax.invert_yaxis(); ax.set_xlim(-0.45,2.45)
ax.set_title("Effect-direction plot of post-denosumab BMD trajectories",fontsize=15)
ax.set_xlabel("BMD site"); ax.set_ylabel("Strategy | independent study unit | post-DMAb time window")
ax.legend(loc="upper right",frameon=False)
prev=None
for yi,k in enumerate(order):
    if prev is not None and k[0]!=prev: ax.axhline(yi-0.5,linewidth=0.8)
    prev=k[0]
fig.tight_layout(); fig.savefig(OUT/"Figure_S1_effect_direction.png",dpi=160); plt.close(fig)

# 5) Hard checks: if these fail, data were changed.
assert summary["SWiM contributing reports"]==43
assert summary["SWiM underlying study units (report/RoB census)"]==24
assert summary["Primary BMD study-level rows"]==103
assert summary["Primary BMD binary directions"]==89
assert summary["Fracture underlying study units (outcome-level MainEffect)"]==17
assert summary["BTM/TBS underlying study units (outcome-level MainEffect)"]==20
assert summary["Structural/microstructural underlying study units (report/RoB census)"]==6
print("Verification complete. Results were written to:", OUT)
print(pd.DataFrame(list(summary.items()),columns=["Metric","Value"]).to_string(index=False))
