#!/usr/bin/env bash
# GSE335219_RAW.tar download + inspect (Aronchik et al., Nature Medicine 2026)
# Step 1 of the gse335219_cnv_check analysis. Re-runnable; safe to run twice.
set -euo pipefail
D="$(cd "$(dirname "$0")/data" && pwd)"
cd "$D"
TAR=GSE335219_RAW.tar
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE335nnn/GSE335219/suppl/GSE335219_RAW.tar"
EXPECTED_BYTES=68229120

echo "== download =="
# NOTE 2026-09-22: plain curl to the GEO endpoint repeatedly died near 92%
# ("curl: (18) transfer closed") with no Range support, restarting from 0
# each retry. wget -c over HTTPS completed the 65 MB file in one go.
if [ -s "$TAR" ] && [ "$(stat -c%s "$TAR")" = "$EXPECTED_BYTES" ]; then
  echo "already present and complete: $TAR"
else
  wget -c -t 20 -w 2 --retry-connrefused "$URL" -O "$TAR"
fi

echo "== file type / size =="
file "$TAR"; stat -c '%s bytes' "$TAR"

echo "== tar listing =="
tar -tf "$TAR"

echo "== extract =="
tar -xf "$TAR"
ls -lh *.maf.gz
