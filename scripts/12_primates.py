"""12. The repeat in other primates: units in the exon of 35 primate genome sequences (28 species; both haplotypes
of the T2T ape genomes) and in the PER3 proteins of 53 species.

The repeat region is found in the translation: it starts at the first L after a stretch that every simian has
([SN]P[AT][AT]T..), and ends at the end of the last unit (GS[AS]AS); its length in amino acids divided by 18 (one
unit) is the number of units. DNA is translated in the three forward frames (the sequences are in the direction of
the gene) and the frame in which the region is found is used. Units are also named after the closest human unit.

Outputs (outputs/tables): primate_units_dna.tsv, primate_units_protein.tsv, primate_units_by_species.tsv
"""
import collections
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import per3lib as L

START = re.compile(r'[SN]P[AT][AT]T[A-Z]{2}L')
END = re.compile(r'GS[AS]AS')
HUMAN_UNITS = {'R1': 'LSTGSPPRENPSHPTASA', 'R2': 'LSTGSPPMKNPSHPTASA', 'R3': 'LSTGSPPMKNPSHPTAST',
               'R4': 'LSMGLPPSRTPSHPTATV', 'R5': 'LSTGSPPSESPSRTGSAA'}


def repeat_region(p):
    m = START.search(p)
    if not m:
        return None
    rest = p[m.end() - 1:]                     # from the L onwards
    e = END.search(rest)
    if not e:
        return None
    return rest[:e.start() + 4]                # up to and including GSAA/GSSA (the last unit's end)


def closest(u):
    def d(a, b):
        return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
    best = min(HUMAN_UNITS, key=lambda k: d(HUMAN_UNITS[k], u))
    return '%s(%d)' % (best, d(HUMAN_UNITS[best], u))


species = {r['species']: r for r in L.read_tsv(os.path.join(L.DATA, 'primates', 'species.tsv'))}
asm_path = os.path.join(L.DATA, 'primates', 'assemblies.tsv')
assemblies = {r['sequence']: r for r in L.read_tsv(asm_path)} if os.path.exists(asm_path) else {}


def long_read(tech, name):
    if re.search(r'PacBio|Nanopore|HiFi|Revio|Sequel', tech or '', re.I):
        return 'yes'
    if re.search(r'^(T2T|NHGRI_m)', name or ''):       # telomere-to-telomere assemblies (long reads)
        return 'yes'
    return 'no' if tech else 'unknown'
desc = L.fasta_descriptions(os.path.join(L.DATA, 'primates', 'PER3_repeat_exon_DNA.fa'))
dna = L.read_fasta(os.path.join(L.DATA, 'primates', 'PER3_repeat_exon_DNA.fa'))
rows = []
for name, s in dna.items():
    sp = name.split('|')[0].replace('_', ' ')
    found = None
    for f in range(3):
        reg = repeat_region(L.translate(s[f:]))
        if reg:
            found = (f, reg)
            break
    f, reg = found if found else (None, '')
    units = [reg[i:i + 18] for i in range(0, len(reg), 18)]
    rows.append(dict(sequence=name, species=sp, common_name=species.get(sp, {}).get('common_name', ''),
                     group=species.get(sp, {}).get('group', ''), source=name.split('|')[2], accession=name.split('|')[1],
                     description=desc.get(name, ''), frame=f, aa=len(reg), units=round(len(reg) / 18, 2),
                     assembly=assemblies.get(name, {}).get('assembly_name', ''),
                     sequencing=assemblies.get(name, {}).get('sequencing_technology', ''),
                     contig_n50=assemblies.get(name, {}).get('contig_n50', ''),
                     long_read_assembly=long_read(assemblies.get(name, {}).get('sequencing_technology', ''),
                                                  assemblies.get(name, {}).get('assembly_name', '')),
                     units_closest_human=' '.join(closest(u) for u in units), region=reg))
L.write_tsv(os.path.join(L.TAB, 'primate_units_dna.tsv'), rows,
            ['sequence', 'species', 'common_name', 'group', 'source', 'accession', 'assembly', 'sequencing',
             'contig_n50', 'long_read_assembly', 'frame', 'aa', 'units', 'units_closest_human', 'region', 'description'])

prot = L.read_fasta(os.path.join(L.DATA, 'primates', 'PER3_proteins.fa'))
prows = []
for name, p in prot.items():
    sp = name.split('|')[0].replace('_', ' ')
    reg = repeat_region(p)
    # outside the simians there is no repeat: report the single unit-like stretch, if the anchors are found
    prows.append(dict(protein=name, species=sp, common_name=species.get(sp, {}).get('common_name', ''),
                      group=species.get(sp, {}).get('group', ''), length=len(p),
                      aa=len(reg) if reg else '', units=round(len(reg) / 18, 2) if reg else '', region=reg or ''))
L.write_tsv(os.path.join(L.TAB, 'primate_units_protein.tsv'), prows,
            ['protein', 'species', 'common_name', 'group', 'length', 'aa', 'units', 'region'])

# by species: units in each genome sequence, and in the protein
by = collections.OrderedDict()
for r in rows:
    by.setdefault(r['species'], dict(species=r['species'], common_name=r['common_name'], group=r['group'],
                                     dna_units=[], protein_units=''))['dna_units'].append('%g' % r['units'])
for r in prows:
    d = by.setdefault(r['species'], dict(species=r['species'], common_name=r['common_name'], group=r['group'],
                                         dna_units=[], protein_units=''))
    d['protein_units'] = '' if r['units'] == '' else '%g' % r['units']
out = []
for d in by.values():
    d = dict(d); d['dna_units'] = '/'.join(d['dna_units']); out.append(d)
L.write_tsv(os.path.join(L.TAB, 'primate_units_by_species.tsv'), out, ['species', 'common_name', 'group', 'dna_units', 'protein_units'])

print('12: units in the genome sequences:')
for r in rows:
    print('    %-30s %-18s %-5s long-read: %-7s %s' % (r['species'], r['source'], r['units'], r['long_read_assembly'],
                                                     r['units_closest_human']))
g = collections.defaultdict(set)
for r in prows:
    g[r['group']].add(r['units'] if r['units'] != '' else 'none')
print('12: units in the proteins, by group:')
for k, v in g.items():
    print('    %-45s %s' % (k, sorted(map(str, v))))

# The region of each protein aligned to the human repeat (residues 981-1070 of NP_001276791.1), from a multiple
# alignment of the 53 proteins with kalign: outside the simians there is no repeat, only a short unit-like segment.
import subprocess
afa = os.path.join(L.WORK, 'PER3_proteins.afa')
with open(os.devnull) as dn:
    subprocess.run(['kalign', '-i', os.path.join(L.DATA, 'primates', 'PER3_proteins.fa'), '-o', afa], stdin=dn,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
aln = L.read_fasta(afa)
hs = aln['Homo_sapiens|NP_001276791.1']
cols, k = [], 0
for c, ch in enumerate(hs):
    if ch != '-':
        k += 1
        if 981 <= k <= 1070:
            cols.append(c)
c0, c1 = cols[0], cols[-1]
arows = []
for n, s_ in aln.items():
    seg = s_[c0:c1 + 1].replace('-', '')
    sp = n.split('|')[0].replace('_', ' ')
    arows.append(dict(protein=n, species=sp, group=species.get(sp, {}).get('group', ''), residues_aligned_to_human_repeat=len(seg),
                      segment=seg, upstream=s_[max(0, c0 - 40):c0].replace('-', '')[-10:],
                      downstream=s_[c1 + 1:c1 + 40].replace('-', '')[:10]))
L.write_tsv(os.path.join(L.TAB, 'primate_aligned_region.tsv'), arows,
            ['protein', 'species', 'group', 'residues_aligned_to_human_repeat', 'upstream', 'segment', 'downstream'])
print('12: residues aligned to the 90-residue human repeat (kalign), by group:')
gg = collections.defaultdict(list)
for r in arows:
    gg[r['group']].append(r['residues_aligned_to_human_repeat'])
for k_, v in gg.items():
    print('    %-45s %s' % (k_, sorted(v)))
