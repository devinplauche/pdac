"""Power analysis for IL-6 ELISA repeat: n per group for 80% power, alpha=0.05,
two-sample t-test on log2 concentrations, using observed triplicate variances.
Uses the larger of the DMSO and treatment SDs (conservative) within each line."""
import numpy as np, pandas as pd
from scipy import stats

XLSX = "/home/hatch/workspace/user/files/LILY_IL-6_ELISA_processed.xlsx"
OUT = "/home/hatch/workspace/pdac/analyses/il6_power_analysis/output"
df = pd.read_excel(XLSX, sheet_name="Results")
rep_cols = ["Y1 (Normalized Conc)", "Y2 (Normalized Conc)", "Y3 (Normalized Conc)"]

rows = []
for sample, g in df.groupby("Sample"):
    vals = g[rep_cols].iloc[0].values.astype(float)
    vals = vals[vals > 0]
    if len(vals) < 2: continue
    log2v = np.log2(vals)
    rows.append((sample, log2v.mean(), log2v.std(ddof=1), len(vals)))
base = pd.DataFrame(rows, columns=["sample","m","sd","n"]).set_index("sample")

def n_needed(delta, sd, alpha=0.05, power=0.80):
    d = abs(delta) / sd  # Cohen's d on log2 scale
    n = 2
    while n < 500:
        ncp = d * np.sqrt(n/2)
        df_ = 2*n - 2
        tcrit = stats.t.ppf(1-alpha/2, df_)
        pw = 1 - stats.nct.cdf(tcrit, df_, ncp) + stats.nct.cdf(-tcrit, df_, ncp)
        if pw >= power: return n
        n += 1
    return n

results = []
for line in ["AZ","BU","GOL 321","GOL 205"]:
    ctrl = base.loc[f"{line} DMSO"]
    for cond in ["10 nM RASi daily","100 nM RASi daily","10 nM RASi 5 Day","100 nM RASi 5 Day"]:
        key = f"{line} {cond}"
        if key not in base.index: continue
        tr = base.loc[key]
        sd = max(ctrl.sd, tr.sd)
        delta = tr.m - ctrl.m
        n = n_needed(delta, sd)
        pw_n3 = 1 - stats.nct.cdf(stats.t.ppf(0.975, 4), 4, abs(delta)/sd*np.sqrt(3/2)) + stats.nct.cdf(-stats.t.ppf(0.975,4), 4, abs(delta)/sd*np.sqrt(3/2))
        results.append({"line":line,"condition":cond,"fold":round(2**delta,2),
                        "delta_log2":round(delta,3),"pooled_sd_log2":round(sd,3),
                        "cohens_d":round(abs(delta)/sd,2),"n_for_80pct":n,
                        "power_at_n3":round(pw_n3,3)})

r = pd.DataFrame(results)
r.to_csv(f"{OUT}/power_table.csv", index=False)
print(r.to_string(index=False))
