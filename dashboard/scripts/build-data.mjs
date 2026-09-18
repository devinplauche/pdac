// Build script: converts the repo's analysis outputs into JSON consumed by the
// Next.js dashboard. Run manually: `npm run build-data`. Outputs are committed
// to src/data/*.json. Missing inputs degrade gracefully (empty data + note).
import { readFileSync, writeFileSync, mkdirSync, existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import YAML from "yaml";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const analyses = resolve(root, "..", "analyses");
const outDir = resolve(root, "src", "data");
mkdirSync(outDir, { recursive: true });

// Minimal RFC-4180 CSV parser (handles quoted fields with commas/quotes).
function parseCSV(text) {
  const rows = [];
  let row = [],
    field = "",
    inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          i++;
        } else inQuotes = false;
      } else field += c;
    } else if (c === '"') inQuotes = true;
    else if (c === ",") {
      row.push(field);
      field = "";
    } else if (c === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else if (c === "\r") {
      /* skip */
    } else field += c;
  }
  if (field !== "" || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [header, ...data] = rows.filter((r) => r.some((x) => x.trim() !== ""));
  return data.map((r) => Object.fromEntries(header.map((h, i) => [h.trim(), (r[i] ?? "").trim()])));
}

function readCSV(p) {
  if (!existsSync(p)) return null;
  try {
    return parseCSV(readFileSync(p, "utf8"));
  } catch (e) {
    console.warn(`build-data: failed to parse ${p}: ${e.message}`);
    return null;
  }
}

const num = (v) => {
  const n = Number(v);
  return Number.isFinite(n) ? n : null;
};

// --- 1. IL-6 power table ---
const powerRows = readCSV(resolve(analyses, "il6_power_analysis/output/power_table.csv"));
const power = {
  generated: new Date().toISOString().slice(0, 10),
  available: powerRows !== null,
  rows: (powerRows ?? []).map((r) => ({
    line: r.line,
    condition: r.condition,
    fold: num(r.fold),
    delta_log2: num(r.delta_log2),
    cohens_d: num(r.cohens_d),
    n_for_80pct: num(r.n_for_80pct),
    power_at_n3: num(r.power_at_n3),
  })),
};
writeFileSync(resolve(outDir, "power.json"), JSON.stringify(power, null, 2));
console.log(`power.json: ${power.rows.length} rows (available=${power.available})`);

// --- 2. Resistance dataset catalog ---
const catRows = readCSV(resolve(analyses, "resistance_datasets_survey/output/catalog.csv"));
const datasets = {
  generated: new Date().toISOString().slice(0, 10),
  available: catRows !== null,
  rows: (catRows ?? []).map((r) => ({
    accession: r.accession,
    source_db: r.source_db,
    publication: r.publication,
    model_system: r.model_system,
    organism: r.organism,
    treatment: r.treatment,
    data_type: r.data_type,
    resistant_model: r.resistant_model,
    raw_available: r.raw_available,
    notes: r.notes,
  })),
};
writeFileSync(resolve(outDir, "datasets.json"), JSON.stringify(datasets, null, 2));
console.log(`datasets.json: ${datasets.rows.length} rows (available=${datasets.available})`);

// --- 3. Signature compendium ---
const sigPath = resolve(analyses, "signature_panel/signatures_v1.yaml");
let signatures = { generated: new Date().toISOString().slice(0, 10), available: false, version: null, signatures: [] };
if (existsSync(sigPath)) {
  try {
    const doc = YAML.parse(readFileSync(sigPath, "utf8"));
    signatures = {
      generated: signatures.generated,
      available: true,
      version: doc.version ?? null,
      curated: doc.curated ?? null,
      scoring_note: doc.scoring_note ?? null,
      signatures: Object.entries(doc.signatures ?? {}).map(([key, s]) => ({
        key,
        label: s.label,
        description: s.description,
        direction: s.direction ?? null,
        human: s.human ?? [],
        mouse: s.mouse ?? [],
        mouse_note: s.mouse_note ?? null,
        human_note: s.human_note ?? null,
        source: {
          first_author: s.source?.first_author,
          year: s.source?.year,
          journal: s.source?.journal,
          pmid: s.source?.pmid ?? null,
          doi: s.source?.doi ?? null,
          pmcid: s.source?.pmcid ?? null,
          note: s.source?.note ?? null,
        },
      })),
    };
  } catch (e) {
    console.warn(`build-data: failed to parse signatures yaml: ${e.message}`);
  }
} else {
  console.warn(`build-data: missing ${sigPath}`);
}
writeFileSync(resolve(outDir, "signatures.json"), JSON.stringify(signatures, null, 2));
console.log(`signatures.json: ${signatures.signatures.length} signatures (available=${signatures.available})`);

console.log("build-data: done.");
