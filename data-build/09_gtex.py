#!/usr/bin/env python3
"""gtex/GTEx_v8_PER3_rs228701.tsv: the effect of rs228701 (chr1_7829744_C_G_b38) on PER3 RNA in each tissue of GTEx v8
that has eQTL data, from the GTEx Portal API (dyneqtl: the normalised effect size, NES, of the ALT allele G, its
p-value, and the number of donors of each genotype).

Run:  python3 09_gtex.py OUTDIR
"""
import csv, json, os, sys, time, urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
API = 'https://gtexportal.org/api/v2/'
get = lambda path: json.load(urllib.request.urlopen(API + path, timeout=120))
tissues = [t['tissueSiteDetailId'] for t in get('dataset/tissueSiteDetail?datasetId=gtex_v8&itemsPerPage=250')['data'] if t.get('hasEGenes')]
rows = []
for t in tissues:                     # in the order of the API
    try:
        r = get(f'association/dyneqtl?gencodeId=ENSG00000049246.14&variantId=chr1_7829744_C_G_b38&tissueSiteDetailId={t}&datasetId=gtex_v8')
    except Exception as e:
        print('no result for', t, e, file=sys.stderr)
        continue
    rows.append([t, r['homoRefCount'], r['hetCount'], r['homoAltCount'], '%.4f' % r['nes'], '%.3g' % r['pValue']])
    time.sleep(0.2)
with open(os.path.join(OUT, 'GTEx_v8_PER3_rs228701.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['tissue', 'n_CC', 'n_CG', 'n_GG', 'NES', 'p_value'])
    w.writerows(rows)
print(len(rows), 'tissues')
