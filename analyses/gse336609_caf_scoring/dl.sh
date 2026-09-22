#!/bin/bash
# Re-download GSE336609 processed matrices from GEO (NCBI HTTPS).
# NOTE: raw data / scRNA-seq were NOT deposited (patient privacy); these are
# the only processed files GEO holds: per-sample count matrices with human
# (tumor) and mouse (stroma) reads extracted separately.
base="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE336nnn/GSE336609/suppl"
cd "$(dirname "$0")/data" || exit 1
for f in \
  GSE336609_AAP01_KRAS-G12D_AA_Human_Reads_MRTX1133_vs_Ctrl.csv.gz \
  GSE336609_AAP01_KRAS-G12D_AA_Mouse_Reads_MRTX1133_vs_Ctrl.csv.gz \
  GSE336609_AAP16_KRAS-G12V_AA_Human_Reads_RMC-6236_vs_Ctrl.csv.gz \
  GSE336609_AAP16_KRAS-G12V_AA_Mouse_Reads_RMC-6236_vs_Ctrl.csv.gz ; do
  if [ -f "$f" ]; then echo "have $f"; continue; fi
  ok=0
  for i in 1 2 3 4 5; do
    if curl -sL --max-time 900 -o "$f" "$base/$f"; then ok=1; break; fi
    echo "retry $i $f"; sleep 5
  done
  [ "$ok" == "0" ] && echo "FAILED $f"
done
echo DONE
ls -la *.csv.gz | awk '{print $5, $9}'
