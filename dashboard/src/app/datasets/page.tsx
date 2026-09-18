"use client";

import { useMemo, useState } from "react";
import datasetsData from "@/data/datasets.json";
import type { DatasetRow } from "@/data/types";
import { EmptyNote } from "@/components/ui";

type Row = DatasetRow;

const ROWS = datasetsData.rows as DatasetRow[];

function matches(r: Row, q: string) {
  if (!q) return true;
  const hay = [r.accession, r.publication, r.model_system, r.treatment, r.data_type, r.notes]
    .join(" ")
    .toLowerCase();
  return q
    .toLowerCase()
    .split(/\s+/)
    .every((t) => hay.includes(t));
}

export default function DatasetsPage() {
  const [query, setQuery] = useState("");
  const [organism, setOrganism] = useState("all");
  const [resistant, setResistant] = useState("all");
  const [open, setOpen] = useState<string | null>(null);

  const organisms = useMemo(
    () => Array.from(new Set(ROWS.map((r) => r.organism))).sort(),
    []
  );

  const filtered = useMemo(
    () =>
      ROWS.filter(
        (r) =>
          matches(r, query) &&
          (organism === "all" || r.organism === organism) &&
          (resistant === "all" ||
            (resistant === "yes" ? /^yes\b/i.test(r.resistant_model) : !/^yes\b/i.test(r.resistant_model)))
      ),
    [query, organism, resistant]
  );

  if (!datasetsData.available) {
    return (
      <div>
        <h1 className="font-serif text-3xl font-semibold mb-4">Dataset catalog</h1>
        <EmptyNote>
          The catalog hasn't been built yet. Run <code>npm run build-data</code> in the
          dashboard directory to generate it from the survey's <code>catalog.csv</code>.
        </EmptyNote>
      </div>
    );
  }

  return (
    <div>
      <h1 className="font-serif text-3xl font-semibold tracking-tight mb-2">Dataset catalog</h1>
      <p className="text-stone-600 mb-6 max-w-3xl leading-relaxed">
        {ROWS.length} curated public datasets around RMC-6236 /
        daraxonrasib and KRAS-inhibitor resistance, from the survey
        (2026-09-18, verification addendum included). Click a row for full details.
      </p>

      <div className="flex flex-col sm:flex-row gap-3 mb-6">
        <input
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Search accession, treatment, model…"
          className="flex-1 rounded-lg border border-stone-300 px-3 py-2 text-sm bg-white focus:outline-none focus:ring-2 focus:ring-accent/40"
        />
        <select
          value={organism}
          onChange={(e) => setOrganism(e.target.value)}
          className="rounded-lg border border-stone-300 px-3 py-2 text-sm bg-white"
        >
          <option value="all">All organisms</option>
          {organisms.map((o) => (
            <option key={o} value={o}>{o}</option>
          ))}
        </select>
        <select
          value={resistant}
          onChange={(e) => setResistant(e.target.value)}
          className="rounded-lg border border-stone-300 px-3 py-2 text-sm bg-white"
        >
          <option value="all">Resistant model: any</option>
          <option value="yes">Resistant model: yes</option>
          <option value="no">Resistant model: no / other</option>
        </select>
      </div>

      <p className="text-xs text-stone-500 mb-3">
        Showing {filtered.length} of {ROWS.length}
      </p>

      <div className="space-y-3">
        {filtered.map((r) => {
          const isOpen = open === r.accession;
          return (
            <article key={r.accession} className="rounded-xl border border-stone-200 bg-white overflow-hidden">
              <button
                onClick={() => setOpen(isOpen ? null : r.accession)}
                className="w-full text-left p-4 hover:bg-stone-50 transition"
              >
                <div className="flex flex-wrap items-center gap-x-3 gap-y-1">
                  <code className="text-sm font-semibold text-accent">{r.accession}</code>
                  <span className="text-xs px-2 py-0.5 rounded-full bg-stone-100 border border-stone-200 text-stone-600">
                    {r.organism}
                  </span>
                  <span className="text-xs px-2 py-0.5 rounded-full bg-stone-100 border border-stone-200 text-stone-600">
                    {r.data_type}
                  </span>
                  {/^yes\b/i.test(r.resistant_model) && (
                    <span className="text-xs px-2 py-0.5 rounded-full bg-emerald-50 border border-emerald-200 text-emerald-700">
                      resistant model
                    </span>
                  )}
                  <span className="ml-auto text-stone-400 text-sm">{isOpen ? "▾" : "▸"}</span>
                </div>
                <p className="text-sm text-stone-600 mt-1.5 leading-relaxed line-clamp-2">{r.model_system}</p>
              </button>
              {isOpen && (
                <div className="px-4 pb-4 pt-1 text-sm border-t border-stone-100">
                  <dl className="grid gap-x-6 gap-y-2.5 mt-3 sm:grid-cols-[140px_1fr]">
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Publication</dt>
                    <dd className="text-stone-700 leading-relaxed">{r.publication}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Model system</dt>
                    <dd className="text-stone-700 leading-relaxed">{r.model_system}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Treatment</dt>
                    <dd className="text-stone-700 leading-relaxed">{r.treatment}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Data type</dt>
                    <dd className="text-stone-700">{r.data_type}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Resistant model</dt>
                    <dd className="text-stone-700">{r.resistant_model}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Raw available</dt>
                    <dd className="text-stone-700">{r.raw_available}</dd>
                    <dt className="text-stone-400 text-xs uppercase tracking-wide pt-0.5">Source</dt>
                    <dd className="text-stone-700">
                      {r.source_db} ·{" "}
                      <a
                        href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${encodeURIComponent(r.accession.split(" ")[0])}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-accent hover:underline"
                      >
                        GEO record
                      </a>
                    </dd>
                  </dl>
                  <p className="mt-3 text-stone-600 leading-relaxed bg-stone-50 rounded-lg p-3 border border-stone-100">
                    {r.notes}
                  </p>
                </div>
              )}
            </article>
          );
        })}
        {filtered.length === 0 && (
          <EmptyNote>No datasets match the current filters.</EmptyNote>
        )}
      </div>
    </div>
  );
}
