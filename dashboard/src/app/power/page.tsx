import powerData from "@/data/power.json";
import type { PowerRow } from "@/data/types";
import { EmptyNote } from "@/components/ui";

function fmt(n: number | null, digits = 2) {
  return n === null ? "—" : n.toFixed(digits).replace(/\.?0+$/, "");
}

export default function PowerPage() {
  if (!powerData.available) {
    return (
      <div>
        <h1 className="font-serif text-3xl font-semibold mb-4">IL-6 ELISA power analysis</h1>
        <EmptyNote>
          The power table hasn't been built yet. Run <code>npm run build-data</code> in the
          dashboard directory to generate it from <code>output/power_table.csv</code>.
        </EmptyNote>
      </div>
    );
  }

  const rows = [...(powerData.rows as PowerRow[])].sort(
    (a, b) => (a.n_for_80pct ?? 999) - (b.n_for_80pct ?? 999)
  );

  return (
    <div>
      <h1 className="font-serif text-3xl font-semibold tracking-tight mb-2">
        IL-6 ELISA: sizing the repeat
      </h1>
      <p className="text-stone-600 mb-6 max-w-3xl leading-relaxed">
        Two-sample t-test on log2 concentrations, α = 0.05, target 80% power,
        conservative variance (max of DMSO / treatment SD). Sorted by required
        replicates — cheap conditions first.
      </p>

      <div className="rounded-xl border border-stone-200 bg-white overflow-hidden mb-8">
        <div className="overflow-x-auto">
          <table className="w-full text-sm min-w-[640px]">
            <thead>
              <tr className="bg-stone-50 text-left text-xs uppercase tracking-wide text-stone-500">
                <th className="px-4 py-3 font-semibold">Line</th>
                <th className="px-4 py-3 font-semibold">Condition</th>
                <th className="px-4 py-3 font-semibold text-right">Fold vs DMSO</th>
                <th className="px-4 py-3 font-semibold text-right">Cohen's d</th>
                <th className="px-4 py-3 font-semibold text-right">n for 80%</th>
                <th className="px-4 py-3 font-semibold text-right">Power at n=3</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((r, i) => {
                const n = r.n_for_80pct ?? 999;
                const tone =
                  n <= 3
                    ? "bg-emerald-50/60"
                    : n <= 6
                      ? "bg-amber-50/60"
                      : "bg-red-50/60";
                return (
                  <tr key={i} className={`border-t border-stone-100 ${tone}`}>
                    <td className="px-4 py-2.5 font-medium">{r.line}</td>
                    <td className="px-4 py-2.5">{r.condition}</td>
                    <td className="px-4 py-2.5 text-right font-mono">{fmt(r.fold)}</td>
                    <td className="px-4 py-2.5 text-right font-mono">{fmt(r.cohens_d)}</td>
                    <td className="px-4 py-2.5 text-right font-mono font-semibold">
                      {r.n_for_80pct ?? "—"}
                    </td>
                    <td className="px-4 py-2.5 text-right font-mono">
                      {r.power_at_n3 === null ? "—" : r.power_at_n3.toFixed(2)}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      </div>

      <h2 className="font-serif text-xl font-semibold mb-3">Takeaways for the repeat</h2>
      <div className="prose-academic text-stone-700 max-w-3xl">
        <ul>
          <li>
            <strong>The 100 nM 5-day hits don't need more replicates.</strong> AZ
            (15.1×), GOL 205 (20.1×) and GOL 321 (2.6×) are at ~100% power with
            triplicates — run them at n=2–3 and spend the plate space elsewhere.
          </li>
          <li>
            <strong>The GOL 321 "weak" 2.6× effect is not weak for power purposes</strong> —
            its triplicates were tight (SD 0.11 on log2 scale), so n=2–3 already
            detects it. The effect to worry about is noise, not size.
          </li>
          <li>
            <strong>AZ 10 nM daily is the expensive one: n=11</strong> — the high-CV
            triplicate condition. Cleaner technicals (or dropping it) saves the most
            wells. Same story for GOL 205 10 nM 5-day (n=25): small 1.6× effect,
            moderate variance.
          </li>
          <li>
            <strong>GOL 321 10 nM daily ("n=24") is not a real target.</strong> One
            near-zero replicate explodes the log2 SD; powering against that outlier
            is meaningless. Re-run cleanly or exclude.
          </li>
          <li>
            BU 5-day conditions are flat (0.88–0.91×) — powered to detect nothing; keep
            only as controls. GOL 205 10 nM daily reads zero (assay floor, no variance
            to power against). BU daily conditions were excluded as essentially flat.
          </li>
          <li>
            Caveat: power is computed from observed (post-hoc) effect sizes, so it is
            optimistic for the big hits — which is exactly why the recommendation is to
            <em> reduce</em> replicates there rather than claim a powered discovery.
          </li>
        </ul>
      </div>
    </div>
  );
}
