"""03. The same haplotypes in a second pangenome graph: HPRC release 2 PGGB (GRCh38-based VCF, 2025-09-30 build).

In the PGGB VCF the four-unit allele is a 54-bp deletion at chr1:7,829,912 (REF 55 bp -> ALT 1 bp); in the
Minigraph-Cactus VCF it is a 54-bp deletion at chr1:7,829,992 (the graphs place the same deletion in different
units of the tandem repeat). Here the units of each haplotype (step 02, from Minigraph-Cactus) are compared with
the PGGB genotype. Also: the units of the CHM13 reference haplotype, and samples present in only one graph.

Output: outputs/tables/graph_concordance.tsv
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per3lib as L

vcf = os.path.join(L.DATA, 'hprc', 'hprc_pggb_per3_7.78-7.86Mb.vcf.gz')
L.run(['tabix', '-f', '-p', 'vcf', vcf]) if not os.path.exists(vcf + '.tbi') else None
samples = L.vcf_samples(vcf)
recs = L.bcftools_query(vcf, '%POS\t%REF\t%ALT[\t%GT]\n', region='GRCh38#0#chr1:7829800-7830300')
# indels whose length is a multiple of the unit change the number of units
len_changing = []
for line in recs:
    f = line.split('\t')
    pos, ref, alt = int(f[0]), f[1], f[2]
    if ',' in alt:
        raise SystemExit('03: multi-allelic record at %d, not expected' % pos)
    d = len(alt) - len(ref)
    if d != 0:
        len_changing.append((pos, ref, alt, d, f[3:]))
units_changing = [x for x in len_changing if x[3] % L.UNIT_LEN == 0]
other_indels = [x for x in len_changing if x[3] % L.UNIT_LEN != 0]
assert len(units_changing) == 1, units_changing
pos, ref, alt, d, gts = units_changing[0]
print('03: PGGB unit-changing record: chr1:%d, REF %d bp, ALT %d bp (%+d bp)' % (pos, len(ref), len(alt), d))
pggb = {}
for s, gt in zip(samples, gts):
    al = gt.replace('/', '|').split('|')
    pggb[s] = [5 + (d // L.UNIT_LEN) * (a == '1') if a in '01' else None for a in al]

mc = {r['haplotype']: int(r['units']) for r in L.read_tsv(os.path.join(L.TAB, 'hprc_haplotype_alleles.tsv'))}
rows, tab = [], collections.Counter()
for h, n in sorted(mc.items()):
    s, hp = h.rsplit('_hap', 1)
    p = pggb.get(s, [None, None])[int(hp) - 1] if s in pggb else None
    tab[(n, p)] += 1
    rows.append(dict(haplotype=h, units_minigraph_cactus=n, units_pggb='' if p is None else p,
                     agree='' if p is None else int(p == n)))
L.write_tsv(os.path.join(L.TAB, 'graph_concordance.tsv'), rows,
            ['haplotype', 'units_minigraph_cactus', 'units_pggb', 'agree'])
both = [r for r in rows if r['units_pggb'] != '']
print('03: haplotypes in both graphs: %d; identical number of units: %d' % (len(both), sum(r['agree'] for r in both)))
print('03: table (MC units, PGGB units): ', dict(tab))
only = sorted(set(pggb) - {h.rsplit('_hap', 1)[0] for h in mc})
print('03: in PGGB only:', {s: pggb[s] for s in only})
print('03: other (non-unit) indels near the repeat in PGGB:', [(x[0], x[3]) for x in other_indels])
