"""Shared helpers for the PER3 VNTR analysis (Python 3 standard library only)."""
import gzip
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(ROOT, 'outputs')
TAB = os.path.join(OUT, 'tables')
WORK = os.path.join(OUT, 'work')

# GRCh38 coordinates (1-based, inclusive)
EXON19 = (7829834, 7830161)          # PER3 exon 19 of ENST00000377532 (MANE Select)
REPEAT = (7829888, 7830157)          # the five 54-bp units of GRCh38
REGION = (7829701, 7830300)          # 600 bp rebuilt for every haplotype
UNIT_LEN = 54
MOTIF = 'CTGTCC'                     # first six bases of every unit (codons Leu-Ser)

CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]


def translate(seq):
    seq = seq.upper()
    return ''.join(CODON.get(seq[i:i + 3], 'X') for i in range(0, len(seq) - 2, 3))


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    seqs, order, name = {}, [], None
    with op(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                order.append(name)
                seqs[name] = []
            else:
                seqs[name].append(line.strip().upper())
    return {k: ''.join(seqs[k]) for k in order}


def fasta_descriptions(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].rstrip('\n').split(None, 1)
                out[parts[0]] = parts[1] if len(parts) > 1 else ''
    return out


def unit_starts(seq, motif=MOTIF, period=UNIT_LEN, tol=3):
    """Start positions (0-based) of the longest chain of motif occurrences spaced period +- tol apart."""
    pos = [m.start() for m in re.finditer('(?=' + motif + ')', seq)]
    best = []
    for p in pos:
        chain = [p]
        while True:
            nxt = [q for q in pos if abs(q - chain[-1] - period) <= tol]
            if not nxt:
                break
            chain.append(nxt[0])
        if len(chain) > len(best):
            best = chain
    return best


def split_units(seq):
    """The units of a VNTR sequence: from each unit start to the next; the last is UNIT_LEN long."""
    st = unit_starts(seq)
    out = []
    for j, s in enumerate(st):
        e = st[j + 1] if j + 1 < len(st) else s + UNIT_LEN
        out.append(seq[s:e])
    return out


def run(cmd, **kw):
    """Run a command (list), return stdout as text; raise on failure."""
    r = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kw)
    return r.stdout


def bcftools_query(vcf, fmt, region=None, samples=None, extra=None):
    cmd = ['bcftools', 'query', '-f', fmt]
    if region:
        cmd += ['-r', region]
    if samples:
        cmd += ['-s', samples]
    if extra:
        cmd += extra
    cmd.append(vcf)
    return run(cmd).splitlines()


def vcf_samples(vcf):
    return run(['bcftools', 'query', '-l', vcf]).split()


def read_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        return [dict(zip(header, line.rstrip('\n').split('\t'))) for line in f if line.strip()]


def write_tsv(path, rows, cols):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')


def dbsnp_ids(region=None):
    """{(pos, ref, alt): rsID} from the bundled dbSNP build 157 file (GRCh38)."""
    vcf = os.path.join(DATA, 'reference', 'dbsnp157_per3_500kb.vcf.gz')
    ids = {}
    for line in bcftools_query(vcf, '%POS\t%REF\t%ALT\t%ID\n', region=region):
        pos, ref, alts, rs = line.split('\t')
        for alt in alts.split(','):
            ids[(int(pos), ref, alt)] = rs
    return ids


def ensure_dirs():
    for d in (OUT, TAB, WORK, os.path.join(OUT, 'figures'), os.path.join(OUT, 'logs')):
        os.makedirs(d, exist_ok=True)
