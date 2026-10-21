"""primates/: find the PER3 VNTR exon in the second haplotype of each T2T ape genome (Yoo et al. 2025).

Run:  python3 10_ape_haplotypes.py OUTDIR HUMAN_EXON19.fa      (needs BLAST+: blastn)
HUMAN_EXON19.fa: the human exon 19, chr1:7829834-7830161 of GRCh38 (samtools faidx of 01_reference.sh's file).
The BLAST hits, extended by about 150 bases on each side and turned to the direction of the gene, are the ape
sequences of PER3_repeat_exon_DNA.fa named "second_haplotype".

For each species: the RefSeq primary chromosome that carries PER3 and the gene's position on it are known
(from the NCBI orthologue download). The other haplotype's matching chromosome is found by name; a window around
the proportional position is fetched with E-utilities and searched with BLAST for the human exon."""
import json, ssl, subprocess, sys, time, urllib.request, os, re

ctx = ssl.create_default_context()
OUT, QUERY = sys.argv[1], os.path.abspath(sys.argv[2])
os.makedirs(OUT, exist_ok=True)
os.chdir(OUT)


def get(url, timeout=300):
    for i in range(4):
        try:
            with urllib.request.urlopen(url, context=ctx, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            print('  retry', i, e)
            time.sleep(3 + 3 * i)
    raise RuntimeError(url)


def seqreport(acc):
    d = json.loads(get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/sequence_reports?page_size=1000'))
    return d.get('reports', [])


# species: (name, primary RefSeq assembly, PER3 chromosome accession, gene start-end on it, other haplotypes)
APES = [
    ('Pan_troglodytes', 'GCF_028858775.2', 'NC_072398.2', 224152972, 224213878, [('alt', 'GCA_028858805.2')]),
    ('Pan_paniscus', 'GCF_029289425.2', 'NC_073249.2', 220381214, 220442021, [('pat', 'GCA_028858825.2'), ('mat', 'GCA_028858845.2')]),
    ('Gorilla_gorilla', 'GCF_029281585.2', 'NC_073224.2', 232772881, 232833717, [('mat', 'GCA_028885495.2'), ('pat', 'GCA_028885475.2')]),
    ('Pongo_pygmaeus', 'GCF_028885625.2', 'NC_072373.2', 223376637, 223437706, [('alt', 'GCA_028885525.2')]),
    ('Pongo_abelii', 'GCF_028885655.2', 'NC_071985.2', 226387495, 226448352, [('alt', 'GCA_028885685.2')]),
    ('Symphalangus_syndactylus', 'GCF_028878055.3', 'NC_072444.2', 29030089, 29092767, [('alt', 'GCA_028878085.3')]),
]

results = []
for sp, pri, chrom, gs, ge, others in APES:
    rep = seqreport(pri)
    row = next(r for r in rep if r.get('refseq_accession') == chrom)
    chr_name, plen = row['chr_name'], int(row['length'])
    print(sp, 'primary', chrom, 'chr', chr_name, plen)
    for tag, acc in others:
        orep = seqreport(acc)
        cand = [r for r in orep if r.get('chr_name') == chr_name and r.get('role') == 'assembled-molecule']
        if not cand:
            print('  no chromosome', chr_name, 'in', acc); continue
        oacc, olen = cand[0]['genbank_accession'], int(cand[0]['length'])
        mid = (gs + ge) // 2
        for scale in (olen / plen, ):
            centre = int(mid * scale)
            lo, hi = max(1, centre - 3_000_000), min(olen, centre + 3_000_000)
            fa = f'{sp}_{tag}_{oacc}_{lo}_{hi}.fa'
            if not os.path.exists(fa):
                url = (f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={oacc}'
                       f'&seq_start={lo}&seq_stop={hi}&rettype=fasta&retmode=text')
                open(fa, 'wb').write(get(url, timeout=900))
            out = subprocess.run(['blastn', '-task', 'blastn', '-query', QUERY, '-subject', fa,
                                  '-outfmt', '6 sseqid pident length qstart qend sstart send evalue bitscore', '-evalue', '1e-20'],
                                 capture_output=True, text=True).stdout.strip()
            print('  ', tag, oacc, olen, 'window', lo, hi, '\n   ', out.replace('\n', '\n    ') or 'NO HIT')
            results.append(dict(species=sp, hap=tag, assembly=acc, chrom=oacc, window_start=lo, window_end=hi, fasta=fa, blast=out))
json.dump(results, open('ape_second_haplotypes.json', 'w'), indent=1)
