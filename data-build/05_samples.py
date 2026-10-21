#!/usr/bin/env python3
"""hprc/hprc_samples.tsv and 1000g/1kgp_samples.tsv: who the people are.

Sources: the HPRC release 2 sample table (GitHub, human-pangenomics/hprc_intermediate_assembly), and for the 1000
Genomes Project the pedigree and population file of the 3,202 high-coverage samples, the list of the 2,504 unrelated
people of phase 3 and the phase 3 table of populations (EBI FTP). The populations are named as in that table (with
its "Bejing" corrected); the HPRC's own names are used only for its people outside the 1000 Genomes populations.

Run:  python3 05_samples.py HPRC_VCF OUTDIR      (HPRC_VCF: hprc_per3_vntr.vcf.gz, for its list of samples)
"""
import csv, io, os, subprocess, sys, urllib.request

VCF, OUT = sys.argv[1], sys.argv[2]
FTP = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/'
URLS = {
    'hprc': 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/sample/hprc_release2_sample_metadata.csv',
    'ped': FTP + 'data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt',
    'panel': FTP + 'release/20130502/integrated_call_samples_v3.20130502.ALL.panel',
    'pops': FTP + 'phase3/20131219.populations.tsv',
}
get = lambda k: urllib.request.urlopen(URLS[k], timeout=120).read().decode()

pops = {}
for r in csv.DictReader(io.StringIO(get('pops')), delimiter='\t'):
    if r.get('Population Code'):
        pops[r['Population Code']] = (r['Population Description'].replace('Bejing', 'Beijing'), r['Super Population'])
ped = [l.split() for l in get('ped').splitlines()[1:] if l.strip()]
unrelated = {l.split('\t')[0] for l in get('panel').splitlines()[1:] if l.strip()}
SEX = {'1': 'male', '2': 'female'}

# the 3,202 people of the 1000 Genomes Project, high coverage
kg = []
for fam, s, father, mother, sex, pop, sup in ped:
    kg.append([s, pop, pops[pop][0], sup, SEX.get(sex, 'unknown'), 'yes' if s in unrelated else 'no', father, mother])
with open(os.path.join(OUT, '1kgp_samples.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['sample', 'population', 'population_name', 'superpopulation', 'sex', 'unrelated_2504', 'father', 'mother'])
    w.writerows(kg)
in_kg = {r[0] for r in kg}

# the people of the HPRC VCF
samples = subprocess.run(['bcftools', 'query', '-l', VCF], capture_output=True, text=True, check=True).stdout.split()
meta = {r['sample_id']: r for r in csv.DictReader(io.StringIO(get('hprc')))}
EXTRA = {'ASL': 'AFR', 'MKK': 'AFR'}                     # HPRC populations outside the 1000 Genomes Project
GIAB = {'HG002': ('AJ', 'Ashkenazi Jewish (Genome in a Bottle)', 'EUR'), 'HG005': ('CHN', 'Han Chinese (Genome in a Bottle)', 'EAS')}
rows = []
for s in samples:
    m = meta[s]
    if s in GIAB:
        code, name, sup = GIAB[s]
    else:
        code = m['population_abbreviation']
        name, sup = pops[code] if code in pops else (m['population_descriptor'], EXTRA[code])
    rows.append([s, code, name, sup, m['sex'], 'yes' if s in in_kg else 'no'])
with open(os.path.join(OUT, 'hprc_samples.tsv'), 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['sample', 'population', 'population_name', 'superpopulation', 'sex', 'in_1000G_high_coverage'])
    w.writerows(rows)
print(len(kg), '1000 Genomes samples;', len(rows), 'HPRC samples')
