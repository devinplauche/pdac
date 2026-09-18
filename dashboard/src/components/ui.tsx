export function StatusBadge({ status }: { status: string }) {
  const done = status === "complete";
  return (
    <span
      className={`inline-block text-xs font-medium px-2 py-0.5 rounded-full border ${
        done
          ? "bg-emerald-50 text-emerald-700 border-emerald-200"
          : "bg-amber-50 text-amber-700 border-amber-200"
      }`}
    >
      {done ? "Complete" : "In progress"}
    </span>
  );
}

export function PubMedLink({ pmid, label }: { pmid: number | string | null | undefined; label?: string }) {
  if (!pmid) return null;
  return (
    <a
      href={`https://pubmed.ncbi.nlm.nih.gov/${pmid}/`}
      target="_blank"
      rel="noopener noreferrer"
      className="text-accent hover:underline"
    >
      {label ?? `PMID ${pmid}`}
    </a>
  );
}

export function SectionTitle({ children }: { children: React.ReactNode }) {
  return <h2 className="font-serif text-xl font-semibold mt-8 mb-3">{children}</h2>;
}

export function EmptyNote({ children }: { children: React.ReactNode }) {
  return (
    <div className="rounded-lg border border-stone-200 bg-stone-50 p-4 text-sm text-stone-600">
      {children}
    </div>
  );
}
