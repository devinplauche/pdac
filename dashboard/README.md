# PDAC CAF Research Dashboard

Next.js (App Router) + Tailwind v4 research dashboard for the PDAC repo —
computational analyses behind a Michigan State PDAC PhD project on
daraxonrasib-treated cancer-associated fibroblasts.

## Develop

```bash
npm install
npm run build-data   # regenerate src/data/*.json from ../analyses outputs
npm run dev
```

## Build / deploy

```bash
npm run build   # static export to ./out — deploy to Vercel as-is
```

`output: "export"` is set in `next.config.ts`, so the build is fully static
(no server required). On Vercel, use the default Next.js preset — it will
serve `out/`.

## Data pipeline

`scripts/build-data.mjs` converts repo analysis outputs into committed JSON:

- `analyses/il6_power_analysis/output/power_table.csv` → `src/data/power.json`
- `analyses/resistance_datasets_survey/output/catalog.csv` → `src/data/datasets.json`
- `analyses/signature_panel/signatures_v1.yaml` → `src/data/signatures.json`

`src/data/analyses.json` is hand-curated (per-analysis summaries, findings,
caveats) from each `report.md`. Figures are copied to `public/figures/`.
Missing inputs degrade gracefully (status note, not a crash).
