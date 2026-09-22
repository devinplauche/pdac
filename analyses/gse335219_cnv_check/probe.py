#!/usr/bin/env python3
"""Probe MAF structure: barcodes, FILTER values, gnomAD_AF, classifications."""
import gzip, sys
from collections import Counter

def probe(path, name):
    print(f"== {name}")
    bc = Counter(); filt = Counter(); cls = Counter(); mstat = Counter()
    af_known = af_zero = af_high = 0
    matched_norm = Counter()
    n = 0
    with gzip.open(path, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if header is None:
                header = row
                print("  ncols:", len(header))
                idx = {c: i for i, c in enumerate(header)}
                continue
            n += 1
            bc[row[idx['Tumor_Sample_Barcode']]] += 1
            filt[row[idx['FILTER']]] += 1
            cls[row[idx['Variant_Classification']]] += 1
            mstat[row[idx['Mutation_Status']]] += 1
            matched_norm[row[idx['Matched_Norm_Sample_Barcode']]] += 1
            af = row[idx['gnomAD_AF']]
            if af and af != '.' and af != '':
                af_known += 1
                try:
                    if float(af) > 0.01: af_high += 1
                except ValueError: pass
            else:
                af_zero += 1
    print("  n variants:", n)
    print("  Tumor_Sample_Barcode:", dict(bc))
    print("  Matched_Norm:", dict(matched_norm))
    print("  Mutation_Status:", dict(mstat))
    print("  FILTER top:", filt.most_common(12))
    print("  Variant_Classification top:", cls.most_common(12))
    print(f"  gnomAD_AF: known={af_known} blank={af_zero} >1%={af_high}")

probe(sys.argv[1], "parental P18 (GSM9808344)")
probe(sys.argv[2], "resistant P28 (GSM9808345)")
