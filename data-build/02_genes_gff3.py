#!/usr/bin/env python3
"""reference/genes_GRCh38.gff3 and reference/PER3_exons.bed: the protein-coding genes of chr1:7.4-8.1 Mb (GRCh38)
with the exons and CDS of their MANE Select transcript, from the Ensembl REST API (GFF3 output, asked for in
100-kb pieces - one request for the whole stretch is refused).

Run:  python3 02_genes_gff3.py OUTDIR
"""
import os, sys, time, urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
URL = 'https://rest.ensembl.org/overlap/region/human/1:{a}-{b}?feature=gene;feature=transcript;feature=exon;feature=cds;content-type=text/x-gff3'


def get(url):
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode()
        except Exception as e:
            print('  retry', i, e, file=sys.stderr)
            time.sleep(2 + 3 * i)
    raise SystemExit('could not fetch ' + url)


rows, seen = [], set()
for a in range(7_400_000, 8_100_000, 100_000):
    for line in get(URL.format(a=a, b=a + 99_999)).splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        key = (f[2], f[3], f[4], f[8])
        if key not in seen:
            seen.add(key)
            rows.append(f)
    time.sleep(0.3)


def attr(f, k):
    for kv in f[8].split(';'):
        if kv.startswith(k + '='):
            return kv[len(k) + 1:]
    return ''


genes = {attr(f, 'ID'): f for f in rows if f[2] == 'gene' and attr(f, 'biotype') == 'protein_coding'}
mane = {attr(f, 'ID'): f for f in rows if attr(f, 'Parent') in genes and 'MANE_Select' in attr(f, 'tag')}
parts = [f for f in rows if f[2] in ('exon', 'CDS') and attr(f, 'Parent') in mane]

out = []
for gid, g in sorted(genes.items(), key=lambda x: int(x[1][3])):
    out.append(g)
    for tid, t in mane.items():
        if attr(t, 'Parent') == gid:
            out.append(t)
            out += sorted((p for p in parts if attr(p, 'Parent') == tid), key=lambda p: (int(p[3]), p[2] != 'exon'))
with open(os.path.join(OUT, 'genes_GRCh38.gff3'), 'w') as fh:
    fh.write('##gff-version 3\n##sequence-region chr1 1 8100000\n')
    fh.write('# Protein-coding genes of GRCh38 chr1:7,400,000-8,100,000 with their MANE Select transcript (Ensembl release via REST, %s)\n' % time.strftime('%B %Y'))
    for f in out:
        fh.write('\t'.join(['chr1'] + f[1:]) + '\n')

# the exons of PER3's MANE Select transcript as BED (0-based starts), numbered along the gene (+ strand)
per3 = [t for t in mane.values() if attr(genes[attr(t, 'Parent')], 'Name') == 'PER3'][0]
exons = sorted((p for p in parts if p[2] == 'exon' and attr(p, 'Parent') == attr(per3, 'ID')), key=lambda p: int(p[3]))
with open(os.path.join(OUT, 'PER3_exons.bed'), 'w') as fh:
    for i, e in enumerate(exons, 1):
        fh.write(f'chr1\t{int(e[3]) - 1}\t{e[4]}\tPER3_exon{i}\t0\t{e[6]}\n')
print(len(genes), 'genes,', len(mane), 'MANE Select transcripts,', len(exons), 'PER3 exons')
