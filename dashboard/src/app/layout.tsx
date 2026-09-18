import type { Metadata } from "next";
import Link from "next/link";
import "./globals.css";

export const metadata: Metadata = {
  title: "PDAC CAF Research Dashboard",
  description:
    "Computational analyses supporting a PDAC PhD thesis on daraxonrasib-treated cancer-associated fibroblasts.",
};

const NAV = [
  { href: "/", label: "Overview" },
  { href: "/signatures", label: "Signatures" },
  { href: "/datasets", label: "Datasets" },
  { href: "/power", label: "Power analysis" },
];

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="min-h-screen flex flex-col">
        <header className="border-b border-stone-200 bg-white/80 backdrop-blur sticky top-0 z-10">
          <div className="max-w-5xl mx-auto px-4 py-3 flex items-center justify-between gap-4">
            <Link href="/" className="font-serif text-lg font-semibold tracking-tight">
              PDAC · CAF Research
            </Link>
            <nav className="flex gap-1 sm:gap-2 text-sm">
              {NAV.map((n) => (
                <Link
                  key={n.href}
                  href={n.href}
                  className="px-2.5 py-1.5 rounded-md text-stone-600 hover:text-ink hover:bg-stone-100"
                >
                  {n.label}
                </Link>
              ))}
            </nav>
          </div>
        </header>
        <main className="flex-1 w-full max-w-5xl mx-auto px-4 py-8">{children}</main>
        <footer className="border-t border-stone-200 mt-12">
          <div className="max-w-5xl mx-auto px-4 py-6 text-xs text-stone-500 leading-relaxed">
            Research dashboard for a Michigan State PDAC PhD project on
            daraxonrasib (RMC-6236)-treated cancer-associated fibroblasts and
            ECM remodeling. All analyses are reproducible — see the repo's
            <code className="mx-1">analyses/</code> directory for code and full reports.
            Gene-symbol provenance verified against PubMed, 2026-09-18.
          </div>
        </footer>
      </body>
    </html>
  );
}
