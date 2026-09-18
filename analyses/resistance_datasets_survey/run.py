"""Survey public omics datasets for daraxonrasib (RMC-6236) resistance.

Queries NCBI GEO (gds), ENA and EBI Search live, dumps raw hits to
output/raw_hits.csv. The curated verdict lives in report.md; the curated
table is output/catalog.csv (human-verified from the raw hits).
"""
import csv
import json
import os
import urllib.parse
import urllib.request

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

QUERIES = {
    # direct drug mentions
    "RMC-6236 OR daraxonrasib": "drug_direct",
    "RAS(ON) inhibitor": "drug_class",
    # resistance-focused
    "daraxonrasib resistance": "resistance_direct",
    "RMC-6236 resistance": "resistance_direct",
    # closest substitutes: PDAC KRAS-inhibitor resistance
    "adagrasib resistance pancreatic": "substitute",
    "sotorasib resistance pancreatic": "substitute",
    "MRTX1133 resistance pancreatic": "substitute",
    "KRAS inhibition resistance pancreatic": "substitute",
}

def get(url, timeout=30):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def geo_search(term, retmax=60):
    q = urllib.parse.urlencode({"db": "gds", "term": term,
                                "retmode": "json", "retmax": retmax})
    d = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
    res = d["esearchresult"]
    return int(res["count"]), res["idlist"]

def geo_summaries(idlist):
    if not idlist:
        return {}
    q = urllib.parse.urlencode({"db": "gds", "id": ",".join(idlist),
                                "retmode": "json"})
    d = json.loads(get(f"{EUTILS}/esummary.fcgi?{q}"))["result"]
    out = {}
    for i in d["uids"]:
        r = d[i]
        out[i] = {
            "accession": r.get("accession", ""),
            "entrytype": r.get("entrytype", ""),
            "title": r.get("title", ""),
        }
    return out

def ena_search(term):
    # ENA portal API: studies matching the term
    q = urllib.parse.urlencode({
        "result": "study",
        "query": f'"{term}"',
        "format": "json",
        "fields": "study_accession,study_title,tax_id",
        "limit": 20,
    })
    try:
        d = json.loads(get(f"https://www.ebi.ac.uk/ena/portal/api/search?{q}"))
    except Exception as e:  # noqa: BLE001 - network best effort
        return ("error", str(e))
    rows = d if isinstance(d, list) else []
    return ("ok", [(r.get("study_accession"), r.get("study_title")) for r in rows])

def ebi_search(term):
    # EBI Search API over ArrayExpress
    q = urllib.parse.urlencode({"query": term, "format": "json", "size": 10})
    try:
        d = json.loads(get(
            f"https://www.ebi.ac.uk/ebisearch/ws/rest/arrayexpress?{q}"))
    except Exception as e:  # noqa: BLE001 - network best effort
        return ("error", str(e))
    hits = [(e["fields"].get("acc", ["?"])[0],
             e["fields"].get("name", ["?"])[0])
            for e in d.get("entries", [])]
    return ("ok", hits)

def main():
    rows = []
    for term, bucket in QUERIES.items():
        try:
            count, ids = geo_search(term)
        except Exception as e:  # noqa: BLE001
            rows.append({"source": "GEO", "bucket": bucket, "query": term,
                         "count": "error", "accession": "", "entrytype": "",
                         "title": str(e)})
            continue
        summ = geo_summaries(ids)
        if not ids:
            rows.append({"source": "GEO", "bucket": bucket, "query": term,
                         "count": count, "accession": "", "entrytype": "",
                         "title": ""})
        for i in ids:
            s = summ.get(i, {})
            rows.append({"source": "GEO", "bucket": bucket, "query": term,
                         "count": count, "accession": s.get("accession"),
                         "entrytype": s.get("entrytype"),
                         "title": (s.get("title") or "")[:200]})

    for term in ["RMC-6236", "daraxonrasib"]:
        status, payload = ena_search(term)
        if status == "error":
            rows.append({"source": "ENA", "bucket": "drug_direct",
                         "query": term, "count": "error", "accession": "",
                         "entrytype": "", "title": payload})
        elif not payload:
            rows.append({"source": "ENA", "bucket": "drug_direct",
                         "query": term, "count": 0, "accession": "",
                         "entrytype": "", "title": ""})
        else:
            for acc, title in payload:
                rows.append({"source": "ENA", "bucket": "drug_direct",
                             "query": term, "count": len(payload),
                             "accession": acc, "entrytype": "study",
                             "title": (title or "")[:200]})
        status, payload = ebi_search(term)
        if status == "error":
            rows.append({"source": "ArrayExpress", "bucket": "drug_direct",
                         "query": term, "count": "error", "accession": "",
                         "entrytype": "", "title": payload})
        elif not payload:
            rows.append({"source": "ArrayExpress", "bucket": "drug_direct",
                         "query": term, "count": 0, "accession": "",
                         "entrytype": "", "title": ""})
        else:
            for acc, title in payload:
                rows.append({"source": "ArrayExpress", "bucket": "drug_direct",
                             "query": term, "count": len(payload),
                             "accession": acc, "entrytype": "experiment",
                             "title": (title or "")[:200]})

    path = os.path.join(OUT, "raw_hits.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["source", "bucket", "query",
                                          "count", "accession", "entrytype",
                                          "title"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
