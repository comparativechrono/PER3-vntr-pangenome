"""02. Read the VNTR in each of the 462 HPRC haplotypes: number of units (two ways), the sequence of each unit,
the allele's structure and lineage (4A, 4B, 5A, 5B), and its protein.

Units are named after the five units of GRCh38 (R1-R5, in order); a unit that differs from all of them is named
after the closest one, with its differences (position within the unit, 1-based: reference>alternative).

Outputs (outputs/tables): hprc_haplotype_alleles.tsv, unit_catalogue.tsv, allele_structures.tsv
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per3lib as L

L.ensure_dirs()
haps = L.read_fasta(os.path.join(L.WORK, 'hprc_haplotypes.fa'))
ref = next(iter(L.read_fasta(os.path.join(L.WORK, 'grch38_region.fa')).values()))
assert len(ref) == L.REGION[1] - L.REGION[0] + 1 == 600

# the reference units and their coordinates
starts = L.unit_starts(ref)
assert len(starts) == 5, starts
ref_units = L.split_units(ref)
assert L.REGION[0] + starts[0] == L.REPEAT[0], (L.REGION[0] + starts[0], L.REPEAT[0])
assert L.REGION[0] + starts[-1] + L.UNIT_LEN - 1 == L.REPEAT[1]
flank = len(ref) - 5 * L.UNIT_LEN                    # 330 bases outside the units
names = {u: 'R%d' % (i + 1) for i, u in enumerate(ref_units)}


def name_unit(u):
    if u in names:
        return names[u]
    best = min(ref_units, key=lambda r: sum(a != b for a, b in zip(r, u)) + abs(len(r) - len(u)))
    diffs = ['%d%s>%s' % (i + 1, a, b) for i, (a, b) in enumerate(zip(best, u)) if a != b]
    if len(u) != len(best):
        diffs.append('len%d' % len(u))
    names[u] = names[best] + '*' + ','.join(diffs)
    return names[u]


LINEAGE = {'R1-R3-R4-R5': '4A', 'R1-R3-R4*8T>C-R5': '4B',
           'R1-R2-R3*9A>G-R4-R5': '5A', 'R1-R2-R3-R4-R5': '5B'}

samples = {r['sample']: r for r in L.read_tsv(os.path.join(L.DATA, 'hprc', 'hprc_samples.tsv'))}
rows, unit_count = [], collections.Counter()
for h, s in haps.items():
    sample, hp = h.rsplit('_hap', 1)
    units = L.split_units(s)
    nm = [name_unit(u) for u in units]
    unit_count.update(units)
    structure = '-'.join(nm)
    prot = ''.join(L.translate(u) for u in units)
    n_motif = s.count(L.MOTIF)
    rows.append(dict(haplotype=h, sample=sample, hap=hp,
                     population=samples[sample]['population'], superpopulation=samples[sample]['superpopulation'],
                     length=len(s), units=len(units), units_motif=n_motif,
                     units_from_length=round((len(s) - flank) / L.UNIT_LEN, 3),
                     structure=structure, lineage=LINEAGE.get(structure, 'other'),
                     protein=prot, residues=len(prot), ser_thr=sum(prot.count(x) for x in 'ST')))

# checks: the two ways of counting units
disagree = [r for r in rows if r['units'] != r['units_from_length']]
assert all(r['units'] == r['units_motif'] for r in rows)
cols = ['haplotype', 'sample', 'hap', 'population', 'superpopulation', 'length', 'units', 'units_motif',
        'units_from_length', 'structure', 'lineage', 'protein', 'residues', 'ser_thr']
L.write_tsv(os.path.join(L.TAB, 'hprc_haplotype_alleles.tsv'), sorted(rows, key=lambda r: r['haplotype']), cols)

# the 1000 Genomes haplotypes with a deletion call in the repeat (step 01), as the short-read panel describes them
kg_path = os.path.join(L.WORK, 'kg_sv_haplotypes.fa')
kg_rows = []
for h, s in (L.read_fasta(kg_path).items() if os.path.exists(kg_path) else []):
    sample, hp = h.rsplit('_hap', 1)
    nm = [name_unit(u) for u in L.split_units(s)]
    kg_rows.append(dict(haplotype=h, sample=sample, hap=hp, length=len(s), units=len(nm), structure='-'.join(nm),
                        lineage=LINEAGE.get('-'.join(nm), 'other')))
L.write_tsv(os.path.join(L.TAB, 'kg_sv_alleles.tsv'), sorted(kg_rows, key=lambda r: r['haplotype']),
            ['haplotype', 'sample', 'hap', 'length', 'units', 'structure', 'lineage'])

# the catalogue of units
human = L.read_fasta(os.path.join(L.DATA, 'primates', 'PER3_proteins.fa'))['Homo_sapiens|NP_001276791.1']
cat = []
for u, n in unit_count.most_common():
    p = L.translate(u)
    cat.append(dict(unit=names[u], dna=u, length=len(u), protein=p, ser_thr=sum(p.count(x) for x in 'ST'),
                    occurrences=n))
cat.sort(key=lambda r: (r['unit'].split('*')[0], '*' in r['unit'], -r['occurrences']))
L.write_tsv(os.path.join(L.TAB, 'unit_catalogue.tsv'), cat, ['unit', 'dna', 'length', 'protein', 'ser_thr', 'occurrences'])

# allele structures
st = collections.Counter((r['units'], r['structure'], r['lineage'], r['ser_thr']) for r in rows)
L.write_tsv(os.path.join(L.TAB, 'allele_structures.tsv'),
            [dict(units=k[0], structure=k[1], lineage=k[2], ser_thr=k[3], haplotypes=v) for k, v in st.most_common()],
            ['units', 'structure', 'lineage', 'ser_thr', 'haplotypes'])

# where the units are in the protein
ref_prot = ''.join(L.translate(u) for u in ref_units)
i = human.find(ref_prot)
print('02: GRCh38 repeat protein (90 aa) at NP_001276791.1 residues %d-%d of %d' % (i + 1, i + len(ref_prot), len(human)))
for u in ref_units:
    p = L.translate(u)
    print('    %s %s  %s  S/T=%d' % (names[u], u, p, sum(p.count(x) for x in 'ST')))
c = collections.Counter(r['units'] for r in rows)
print('02: haplotypes by units:', dict(sorted(c.items())), '| length method disagrees for', len(disagree),
      [(r['haplotype'], r['length']) for r in disagree])
print('02: lineages:', dict(collections.Counter(r['lineage'] for r in rows).most_common()))
print('02: 1000 Genomes haplotypes with a deletion call in the repeat, by structure:',
      dict(collections.Counter(r['structure'] for r in kg_rows).most_common()))
for k, v in st.most_common():
    print('    %4d  %s units  %-32s %-6s S/T=%d' % (v, k[0], k[1], k[2], k[3]))
