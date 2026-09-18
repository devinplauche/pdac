import { notFound } from "next/navigation";
import Link from "next/link";
import analysesData from "@/data/analyses.json";
import type { AnalysisEntry } from "@/data/types";
import { StatusBadge, SectionTitle } from "@/components/ui";

const ANALYSES = analysesData.analyses as AnalysisEntry[];

export function generateStaticParams() {
  return ANALYSES.map((a) => ({ slug: a.slug }));
}

export default async function AnalysisPage({
  params,
}: {
  params: Promise<{ slug: string }>;
}) {
  const { slug } = await params;
  const a = ANALYSES.find((x) => x.slug === slug);
  if (!a) notFound();

  return (
    <div>
      <Link href="/" className="text-sm text-stone-500 hover:text-accent">
        ← All analyses
      </Link>
      <div className="flex flex-wrap items-center gap-3 mt-2 mb-3">
        <h1 className="font-serif text-3xl font-semibold tracking-tight">{a.title}</h1>
        <StatusBadge status={a.status} />
      </div>
      <p className="text-stone-600 leading-relaxed max-w-3xl mb-6">{a.summary}</p>

      <div className="rounded-xl border border-stone-200 bg-white p-5 mb-6">
        <h2 className="text-xs font-semibold uppercase tracking-wide text-stone-400 mb-2">
          Research question
        </h2>
        <p className="text-stone-700 leading-relaxed">{a.question}</p>
      </div>

      {a.figure && (
        <figure className="mb-6">
          {/* plain img: static export has no image optimizer */}
          <img
            src={a.figure}
            alt={a.figure_caption ?? a.title}
            className="rounded-xl border border-stone-200 bg-white w-full"
          />
          {a.figure_caption && (
            <figcaption className="text-xs text-stone-500 mt-2 leading-relaxed">
              {a.figure_caption}
            </figcaption>
          )}
        </figure>
      )}

      <SectionTitle>Key findings</SectionTitle>
      <div className="prose-academic text-stone-700 max-w-3xl">
        <ul>
          {a.findings.map((f, i) => (
            <li key={i}>{f}</li>
          ))}
        </ul>
      </div>

      {a.caveats.length > 0 && (
        <>
          <SectionTitle>Caveats</SectionTitle>
          <div className="prose-academic text-stone-600 max-w-3xl text-sm">
            <ul>
              {a.caveats.map((c, i) => (
                <li key={i}>{c}</li>
              ))}
            </ul>
          </div>
        </>
      )}
    </div>
  );
}
