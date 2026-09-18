# IL-6 ELISA power analysis — sizing the repeat

**Source:** `~/workspace/user/files/LILY_IL-6_ELISA_processed.xlsx` (same data as the IL-6 fold-change report, 2026-09-17).
**Method:** two-sample t-test on log2-transformed concentrations, α = 0.05, target 80% power.
Variance = max(DMSO SD, treatment SD) within each line (conservative). Full table: `output/power_table.csv`.

## Answer, per condition (n per group for 80% power)

| Line | Condition | Fold vs DMSO | n for 80% | Power at n=3 |
|---|---|---|---|---|
| AZ | 100 nM RASi 5 Day | 15.1× | 2 | 1.00 |
| GOL 205 | 100 nM RASi 5 Day | 20.1× | 2 | 1.00 |
| GOL 321 | 100 nM RASi 5 Day | 2.6× | 2 | 1.00 |
| GOL 321 | 100 nM RASi daily | 0.36× | 3 | 0.91 |
| GOL 321 | 10 nM RASi 5 Day | 0.25× | 3 | 0.94 |
| BU | 100 nM RASi 5 Day | 0.91× | 3 | 0.96 |
| AZ | 100 nM RASi daily | 3.5× | 4 | 0.72 |
| GOL 205 | 100 nM RASi daily | 2.6× | 4 | 0.60 |
| AZ | 10 nM RASi 5 Day | 0.29× | 5 | 0.59 |
| BU | 10 nM RASi 5 Day | 0.88× | 6 | 0.45 |
| AZ | 10 nM RASi daily | 3.1× | 11 | 0.23 |
| GOL 205 | 10 nM RASi 5 Day | 1.6× | 25 | 0.12 |

Excluded: GOL 205 10 nM daily (at assay floor — reads zero, no variance to power against);
BU 10/100 nM daily (essentially flat, n=2, nothing to detect).

## Takeaways for the repeat

- **The 100 nM 5-day hits don't need more replicates.** All three lines (15×, 20×, 2.6×) are already
  at ~100% power with triplicates; she can run those at n=3 (or even 2) and spend the plate space elsewhere.
- **The GOL 321 "weak" 2.6× effect is not weak for power purposes** — triplicates were tight (SD 0.11 on
  log2 scale), so n=2–3 already detects it. The effect to worry about is noise, not size.
- **AZ 10 nM daily is the expensive one: n=11.** That's the high-CV triplicate condition — a repeat with
  cleaner technicals (or dropping that condition) saves the most wells. Same story for GOL 205 10 nM 5-day
  (n=25), where the effect (1.6×) is small and variance moderate.
- **GOL 321 10 nM daily (n=24) is not a real target.** One replicate reads near-zero, which explodes the
  log2 SD; powering against that outlier is meaningless. Exclude or re-run cleanly.
- BU 5-day conditions are flat (0.88–0.91×) — powered to detect nothing; keep them only as controls.

## Caveats

- Power computed from observed (post-hoc) effect sizes — optimistic for the big hits, which is exactly
  why the recommendation is to *reduce* replicates there rather than claim a powered discovery.
- Log2 t-test assumes roughly log-normal data; fine for ELISA concentrations in range.
- The BU line's DMSO and all BU treatments sit near the assay floor; fold changes there are noisy.
