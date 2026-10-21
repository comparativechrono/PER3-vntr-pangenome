#!/usr/bin/env python3
"""gtex/GTEx_v8_PER3_lineage_markers.tsv: the effect on PER3 RNA (GTEx v8, every tissue with eQTL data) of the
variants that mark the repeat's lineages and sub-lineages, from the GTEx Portal API v2 (association/dyneqtl: the
normalised effect size, NES, of the ALT allele, its p-value, and the number of donors of each genotype).

  rs1773135   chr1_7830204_A_G_b38   repeat length (A = five units; the second tag, with rs228701)
  rs228663    chr1_7832049_T_A_b38   4B (A)
  rs61773390  chr1_7824465_G_T_b38   4B.1 (T)
  rs228697    chr1_7827519_C_G_b38   4B.2 (G)

Run:  python3 11_gtex_lineage_markers.py OUTDIR      (made on 21 October 2026)
"""
import csv, json, os, sys, time, urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
API = 'https://gtexportal.org/api/v2/'
VARIANTS = [('rs1773135', 'chr1_7830204_A_G_b38'), ('rs228663', 'chr1_7832049_T_A_b38'),
            ('rs61773390', 'chr1_7824465_G_T_b38'), ('rs228697', 'chr1_7827519_C_G_b38')]


def get(path, tries=4):
    for k in range(tries):
        try:
            return json.load(urllib.request.urlopen(API + path, timeout=120))
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


tissues = [t['tissueSiteDetailId'] for t in get('dataset/tissueSiteDetail?datasetId=gtex_v8&itemsPerPage=250')['data']
           if t.get('hasEGenes')]
rows = []
for rs, vid in VARIANTS:
    for t in tissues:
        try:
            r = get(f'association/dyneqtl?gencodeId=ENSG00000049246.14&variantId={vid}&tissueSiteDetailId={t}&datasetId=gtex_v8')
        except Exception as e:
            print('no result for', rs, t, e, file=sys.stderr)
            continue
        rows.append([rs, vid, t, r['homoRefCount'], r['hetCount'], r['homoAltCount'], '%.4f' % r['nes'], '%.3g' % r['pValue']])
        time.sleep(0.15)
with open(os.path.join(OUT, 'GTEx_v8_PER3_lineage_markers.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['rsid', 'variant_id', 'tissue', 'n_hom_ref', 'n_het', 'n_hom_alt', 'NES', 'p_value'])
    w.writerows(rows)
print(len(rows), 'rows')
