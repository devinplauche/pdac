import signaturesData from "@/data/signatures.json";
import type { SignatureEntry } from "@/data/types";
import { PubMedLink, SectionTitle, EmptyNote } from "@/components/ui";

const SIGNATURES = signaturesData.signatures as SignatureEntry[];

function GeneChips({ genes, tone }: { genes: string[]; tone: "human" | "mouse" }) {
  return (
    <div className="flex flex-wrap gap-1.5">
      {genes.map((g) => (
        <span
          key={g}
          className={`text-xs font-mono px-2 py-1 rounded-md border ${
            tone === "human"
              ? "bg-blue-50 border-blue-200 text-blue-900"
              : "bg-violet-50 border-violet-200 text-violet-900"
          }`}
        >
          {g}
        </span>
      ))}
    </div>
  );
}

export default function SignaturesPage() {
  if (!signaturesData.available) {
    return (
      <div>
        <h1 className="font-serif text-3xl font-semibold mb-4">Signature compendium</h1>
        <EmptyNote>
          The signature data hasn't been built yet. Run <code>npm run build-data</code> in the
          dashboard directory to generate it from <code>analyses/signature_panel/signatures_v1.yaml</code>.
        </EmptyNote>
      </div>
    );
  }
  return (
    <div>
      <h1 className="font-serif text-3xl font-semibold tracking-tight mb-2">
        Signature compendium
      </h1>
      <p className="text-stone-600 mb-1 max-w-3xl leading-relaxed">
        Version {signaturesData.version} · curated {signaturesData.curated} ·
        {" "}{SIGNATURES.length} signatures with human and mouse
        gene lists. Every provenance verified against PubMed.
      </p>
      <p className="text-sm text-stone-500 max-w-3xl mb-6">{signaturesData.scoring_note}</p>

      <div className="space-y-6">
        {SIGNATURES.map((s) => (
          <article key={s.key} className="rounded-xl border border-stone-200 bg-white p-5 sm:p-6">
            <div className="flex flex-wrap items-baseline gap-x-3 gap-y-1 mb-2">
              <h2 className="font-serif text-xl font-semibold">{s.label}</h2>
              <code className="text-xs text-stone-400">{s.key}</code>
              {s.direction && (
                <span className="text-xs italic text-stone-500">{s.direction}</span>
              )}
            </div>
            <p className="text-sm text-stone-600 leading-relaxed mb-4 max-w-3xl">{s.description}</p>

            <SectionTitle>Genes</SectionTitle>
            <p className="text-xs font-semibold uppercase tracking-wide text-stone-400 mb-2">
              Human ({s.human.length})
            </p>
            <GeneChips genes={s.human} tone="human" />
            {s.human_note && <p className="text-xs text-stone-500 mt-2 italic">{s.human_note}</p>}
            <p className="text-xs font-semibold uppercase tracking-wide text-stone-400 mb-2 mt-4">
              Mouse ({s.mouse.length})
            </p>
            <GeneChips genes={s.mouse} tone="mouse" />
            {s.mouse_note && <p className="text-xs text-stone-500 mt-2 italic">{s.mouse_note}</p>}

            <SectionTitle>Provenance</SectionTitle>
            <p className="text-sm text-stone-600 leading-relaxed">
              {s.source.first_author} et al., <em>{s.source.journal}</em> {s.source.year}
              {" · "}
              <PubMedLink pmid={s.source.pmid} />
              {s.source.doi && (
                <>
                  {" · "}
                  <a
                    href={`https://doi.org/${s.source.doi}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-accent hover:underline"
                  >
                    DOI
                  </a>
                </>
              )}
            </p>
            {s.source.note && (
              <p className="text-xs text-stone-500 mt-2 leading-relaxed max-w-3xl">{s.source.note}</p>
            )}
          </article>
        ))}
      </div>
    </div>
  );
}
