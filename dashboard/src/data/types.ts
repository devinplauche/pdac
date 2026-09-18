export interface AnalysisEntry {
  slug: string;
  title: string;
  status: string;
  summary: string;
  question: string;
  findings: string[];
  caveats: string[];
  figure: string | null;
  figure_caption: string | null;
}

export interface PowerRow {
  line: string;
  condition: string;
  fold: number | null;
  delta_log2: number | null;
  cohens_d: number | null;
  n_for_80pct: number | null;
  power_at_n3: number | null;
}

export interface DatasetRow {
  accession: string;
  source_db: string;
  publication: string;
  model_system: string;
  organism: string;
  treatment: string;
  data_type: string;
  resistant_model: string;
  raw_available: string;
  notes: string;
}

export interface SignatureSource {
  first_author?: string;
  year?: number;
  journal?: string;
  pmid?: number | null;
  doi?: string | null;
  pmcid?: string | null;
  note?: string | null;
}

export interface SignatureEntry {
  key: string;
  label: string;
  description: string;
  direction: string | null;
  human: string[];
  mouse: string[];
  mouse_note: string | null;
  human_note: string | null;
  source: SignatureSource;
}
