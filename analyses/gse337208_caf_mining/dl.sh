#!/bin/bash
# Download GSE337208 per-sample count matrices from GEO FTP.
base="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9850nnn"
declare -A f=( [GSM9850794]=343P [GSM9850795]=MRTX1133343P [GSM9850796]=KPC-T-CTLA4 [GSM9850797]=KPC-T-MRTX [GSM9850798]=KPC-T-1133-CTLA4 [GSM9850799]=Hy-RMC6236 [GSM9850800]=Hy-Veh [GSM9850801]=Hy-RMC6236-aCTLA4 [GSM9850802]=Hy-RMC6236-aPD1 )
declare -A sz=( [GSM9850794_343P_matrix.mtx.gz]=130453108 [GSM9850795_MRTX1133343P_matrix.mtx.gz]=84475067 [GSM9850796_KPC-T-CTLA4_matrix.mtx.gz]=67314221 [GSM9850797_KPC-T-MRTX_matrix.mtx.gz]=61602910 [GSM9850798_KPC-T-1133-CTLA4_matrix.mtx.gz]=51504014 [GSM9850799_Hy-RMC6236_matrix.mtx.gz]=69024361 [GSM9850800_Hy-Veh_matrix.mtx.gz]=68214748 [GSM9850801_Hy-RMC6236-aCTLA4_matrix.mtx.gz]=56610218 [GSM9850802_Hy-RMC6236-aPD1_matrix.mtx.gz]=58950968 )
cd "$(dirname "$0")/data" || exit 1
for gsm in "${!f[@]}"; do
  pfx="${f[$gsm]}"
  for suf in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
    out="${gsm}_${pfx}_${suf}"
    if [ -f "$out" ]; then
      exp="${sz[$out]}"
      if [ -z "$exp" ] || [ "$(stat -c%s "$out")" == "$exp" ]; then continue; fi
      echo "re-downloading partial $out"; rm -f "$out"
    fi
    ok=0
    for i in 1 2 3 4 5; do
      if curl -s --max-time 900 -o "$out" "$base/$gsm/suppl/${gsm}_${pfx}_${suf}"; then ok=1; break; fi
      echo "retry $i $out"; sleep 5
    done
    [ "$ok" == "0" ] && echo "FAILED $out"
  done
done
echo DONE
ls -la *_matrix.mtx.gz | awk '{print $5, $9}'
