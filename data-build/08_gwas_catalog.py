#!/usr/bin/env python3
"""gwas/gwas_catalog_PER3_associations.tsv: the associations that the GWAS Catalog maps to PER3 (REST API v2).

Run:  python3 08_gwas_catalog.py OUTDIR
The Catalog grows: run again, and the table may have more rows than the 28 of 21 October 2026.
"""
import csv, json, os, sys, urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
url = 'https://www.ebi.ac.uk/gwas/rest/api/v2/associations?mapped_gene=PER3&page=0&size=200'
d = json.load(urllib.request.urlopen(url, timeout=120))
rows = []
for a in d['_embedded']['associations']:
    allele = a['snp_allele'][0] if a.get('snp_allele') else {'rs_id': '', 'effect_allele': ''}
    pos = a['locations'][0].split(':')[1] if a.get('locations') else ''
    rows.append([allele['rs_id'], pos, allele.get('effect_allele') or '?', '%g' % a['p_value'],
                 '; '.join(t['efo_trait'] for t in a.get('efo_traits', [])), '; '.join(a.get('reported_trait', [])),
                 a['accession_id'], a.get('pubmed_id', ''), a.get('first_author', '')])
with open(os.path.join(OUT, 'gwas_catalog_PER3_associations.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['rsid', 'position_GRCh38', 'effect_allele', 'p_value', 'trait', 'reported_trait', 'study', 'pubmed_id', 'first_author'])
    w.writerows(rows)
print(len(rows), 'associations')
