#!/usr/bin/env python3
"""primates/assemblies.tsv: the genome assembly behind each primate DNA sequence of PER3_repeat_exon_DNA.fa, and how
it was made (sequencing technologies, assembler, level, contig N50) - to judge how far a tandem repeat in it can be
trusted. NCBI E-utilities (elink nuccore -> assembly) and the NCBI Datasets API v2 (genome dataset_report).

Run:  python3 12_primate_assemblies.py PER3_repeat_exon_DNA.fa OUTDIR      (made on 21 October 2026)
"""
import csv, json, os, sys, time, urllib.parse, urllib.request

fa, OUT = sys.argv[1], sys.argv[2]
EU = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
DS = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/dataset_report'


def get(url, tries=4):
    for k in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=120))
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


names = [l[1:].split()[0] for l in open(fa) if l.startswith('>')]
rows = []
for n in names:
    sp, acc, src = n.split('|')
    link = get(EU + 'elink.fcgi?' + urllib.parse.urlencode({'dbfrom': 'nuccore', 'db': 'assembly', 'id': acc, 'retmode': 'json'}))
    ids = [i for ls in link['linksets'] for db in ls.get('linksetdbs', []) for i in db.get('links', [])]
    summ = get(EU + 'esummary.fcgi?' + urllib.parse.urlencode({'db': 'assembly', 'id': ','.join(ids), 'retmode': 'json'})) if ids else {'result': {'uids': []}}
    cands = [summ['result'][u] for u in summ['result']['uids']]
    # prefer the RefSeq (GCF_) record
    cands.sort(key=lambda r: (not r.get('assemblyaccession', '').startswith('GCF_'), r.get('assemblyaccession', '')))
    a = cands[0] if cands else {}
    asm = a.get('assemblyaccession', '')
    rep = {}
    if asm:
        try:
            rep = (get(DS % asm).get('reports') or [{}])[0]
        except Exception:
            rep = {}
    ai, st = rep.get('assembly_info', {}), rep.get('assembly_stats', {})
    rows.append([n, sp.replace('_', ' '), acc, src, asm, ai.get('assembly_name', a.get('assemblyname', '')),
                 ai.get('assembly_level', a.get('assemblystatus', '')), ai.get('sequencing_tech', ''),
                 ai.get('assembly_method', ''), st.get('contig_n50', a.get('contign50', '')),
                 ai.get('release_date', a.get('seqreleasedate', '')[:10])])
    time.sleep(0.4)
with open(os.path.join(OUT, 'assemblies.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['sequence', 'species', 'sequence_accession', 'source', 'assembly_accession', 'assembly_name',
                'assembly_level', 'sequencing_technology', 'assembler', 'contig_n50', 'release_date'])
    w.writerows(rows)
print(len(rows), 'sequences')
