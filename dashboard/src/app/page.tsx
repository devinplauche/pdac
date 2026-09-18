import Link from "next/link";
import analysesData from "@/data/analyses.json";
import { StatusBadge } from "@/components/ui";

export default function Home() {
  const analyses = analysesData.analyses;
  return (
    <div>
      <h1 className="font-serif text-3xl font-semibold tracking-tight mb-3">
        Daraxonrasib-treated CAFs in PDAC
      </h1>
      <div className="prose-academic text-stone-700 max-w-3xl">
        <p>
          This dashboard collects the computational analyses behind a Michigan
          State PDAC PhD project studying <strong>daraxonrasib (RMC-6236)</strong>-treated
          cancer-associated fibroblasts (CAFs) and ECM remodeling. Each card
          below links to a per-analysis summary with key findings, honest
          caveats, and the underlying figures.
        </p>
        <p>
          The analyses live in the repo's <code>analyses/</code> directory and
          are fully reproducible. The gene-signature compendium and dataset
          catalog are browsable in their own sections above.
        </p>
      </div>

      <h2 className="font-serif text-xl font-semibold mt-10 mb-4">Analyses</h2>
      <div className="grid gap-4 sm:grid-cols-2">
        {analyses.map((a) => (
          <Link
            key={a.slug}
            href={`/analyses/${a.slug}`}
            className="block rounded-xl border border-stone-200 bg-white p-5 hover:shadow-md hover:border-stone-300 transition"
          >
            <div className="flex items-start justify-between gap-2 mb-2">
              <h3 className="font-semibold leading-snug">{a.title}</h3>
              <StatusBadge status={a.status} />
            </div>
            <p className="text-sm text-stone-600 leading-relaxed">{a.summary}</p>
            <span className="inline-block mt-3 text-sm text-accent">Read summary →</span>
          </Link>
        ))}
      </div>

      <h2 className="font-serif text-xl font-semibold mt-10 mb-4">Quick links</h2>
      <div className="grid gap-4 sm:grid-cols-3">
        <Link href="/signatures" className="rounded-xl border border-stone-200 bg-white p-5 hover:shadow-md transition">
          <h3 className="font-semibold mb-1">Signature compendium</h3>
          <p className="text-sm text-stone-600">7 versioned CAF / RAS-inhibition signatures with gene lists and PubMed provenance.</p>
        </Link>
        <Link href="/datasets" className="rounded-xl border border-stone-200 bg-white p-5 hover:shadow-md transition">
          <h3 className="font-semibold mb-1">Dataset catalog</h3>
          <p className="text-sm text-stone-600">16 curated public datasets around RMC-6236 and KRAS-inhibitor resistance — searchable.</p>
        </Link>
        <Link href="/power" className="rounded-xl border border-stone-200 bg-white p-5 hover:shadow-md transition">
          <h3 className="font-semibold mb-1">Power analysis</h3>
          <p className="text-sm text-stone-600">Replicates needed per IL-6 ELISA condition for 80% power.</p>
        </Link>
      </div>
    </div>
  );
}
