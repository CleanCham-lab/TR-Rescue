#!/usr/bin/env python3
"""Synthetic data for tr-rescue `trace` tests.

Builds a small reference with repeated Alu-like elements (one inside a flank), three
samples (biallelic GAA expansion with partial survival in fail reads; normal; biallelic
AAGGG expansion that mostly vanishes), an aligned + haplotagged HiFi BAM per sample, and a
Revio-like run folder with unaligned fail_reads / hifi_reads BAMs carrying ff/rq/np tags.

Usage: make_synthetic.py OUTDIR
"""
import os
import random
import sys

import pysam

R = random.Random(7)
RX = random.Random(11)  # chrX uses its own stream so the autosomal data stay identical
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_seq(n):
    return "".join(R.choice("ACGT") for _ in range(n))


# ---------------------------------------------------------------- reference
CLEN = 400_000
alu = rand_seq(300)
contigs = {"chrA": list(rand_seq(CLEN)), "chrB": list(rand_seq(CLEN)), "chrC": list(rand_seq(CLEN))}
# scatter 60 Alu copies (some reverse-complemented)
for i in range(60):
    c = R.choice(list(contigs))
    p = R.randrange(1000, CLEN - 1000)
    e = alu if i % 2 else rc(alu)
    contigs[c][p:p + 300] = list(e)
# one Alu copy inside the left flank of the GAA locus (must be filtered out of the panel)
LOC_A = 200_000
contigs["chrA"][LOC_A - 1500:LOC_A - 1200] = list(alu)
# GAA locus: (GAA)9 at chrA:200000-200027 ; AAAAG locus: (AAAAG)11 at chrB:200000-200055
gaa_ref = "GAA" * 9
contigs["chrA"][LOC_A:LOC_A + len(gaa_ref)] = list(gaa_ref)
LOC_B = 200_000
rfc_ref = "AAAAG" * 11
contigs["chrB"][LOC_B:LOC_B + len(rfc_ref)] = list(rfc_ref)
# X-linked CGG locus (normal 20 repeats) for sex / copy-number tests
contigs["chrX"] = [RX.choice("ACGT") for _ in range(CLEN)]
LOC_X = 200_000
cgg_ref = "CGG" * 20
contigs["chrX"][LOC_X:LOC_X + len(cgg_ref)] = list(cgg_ref)
contigs = {k: "".join(v) for k, v in contigs.items()}
# soft-mask a stretch to check case handling
contigs["chrC"] = contigs["chrC"][:5000] + contigs["chrC"][5000:9000].lower() + contigs["chrC"][9000:]

# an ALT contig carrying the TRFC neighbourhood (flank k-mers must be reported as ALT-shared)
contigs["chrB_KI000001v1_alt"] = contigs["chrB"][LOC_B - 8000:LOC_B + 8000]

ref_path = os.path.join(OUT, "ref.fa")
with open(ref_path, "w") as f:
    for name, s in contigs.items():
        f.write(f">{name}\n")
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
pysam.faidx(ref_path)
with open(os.path.join(OUT, "catalog.bed"), "w") as f:
    f.write(f"chrA\t{LOC_A}\t{LOC_A + len(gaa_ref)}\tID=TGAA;MOTIFS=GAA;NORMAL_MAX=33;PATH_MIN=66\n")
    f.write(f"chrX\t{LOC_X}\t{LOC_X + len(cgg_ref)}\tID=TX;MOTIFS=CGG;NORMAL_MAX=44;PATH_MIN=200\n")
    f.write(f"chrB\t{LOC_B}\t{LOC_B + len(rfc_ref)}\tID=TRFC;MOTIFS=AAAAG,AAAGG,AAGGG;NORMAL_MAX=20;PATH_MIN=400;PATH_MOTIFS=AAGGG\n")

LOCI = {
    "chrA": (LOC_A, LOC_A + len(gaa_ref)),
    "chrB": (LOC_B, LOC_B + len(rfc_ref)),
    "chrX": (LOC_X, LOC_X + len(cgg_ref)),
}

# ---------------------------------------------------------------- samples
# alleles per (contig) -> (hap1 repeat seq, hap2 repeat seq), plus survival for expanded molecules
SAMPLES = {
    "P1": {"sex": "M", "bc": "bc2001", "alleles": {"chrA": ("GAA" * 600, "GAA" * 800)}, "exp_fate": (0.15, 0.45, 0.10)},
    "N1": {"sex": "F", "bc": "bc2002", "alleles": {}, "exp_fate": None},
    "P2": {"sex": "F", "bc": "bc2003", "alleles": {"chrB": ("AAGGG" * 1000, "AAGGG" * 900)}, "exp_fate": (0.03, 0.05, 0.02)},
    # one-allele expansion, mostly lost; HiPhase block only +-40 kb around the locus
    "P3": {"sex": "M", "bc": "bc2004", "alleles": {"chrA": ("GAA" * 700, "GAA" * 9)}, "exp_fate": (0.05, 0.10, 0.02),
           "hp_window": ("chrA", 160_000, 240_000)},
}
MOVIE = "m84000_250101_000000_s1"
# P2 sits alone in its own SMRT Cell (single-plex), so its unassigned BAMs are safe to use
MOVIE2 = "m84000_250102_000000_s2"
SINGLE = {"P2": MOVIE2}


def mutate(seq, rate, R=R):
    out = []
    for ch in seq:
        r = R.random()
        if r < rate * 0.4:
            continue  # deletion
        if r < rate * 0.8:
            out.append(ch)
            out.append(R.choice("ACGT"))  # insertion
            continue
        if r < rate:
            out.append(R.choice("ACGT"))  # substitution
            continue
        out.append(ch)
    return "".join(out)


def hap_contig(sample, contig, hap):
    s = contigs[contig].upper()
    al = SAMPLES[sample]["alleles"].get(contig)
    if not al:
        return s, None
    st, en = LOCI[contig]
    a = al[hap - 1]
    return s[:st] + a + s[en:], (st, st + len(a), en - st)


def write_header(path_like_sample, movie=MOVIE):
    return {
        "HD": {"VN": "1.6", "SO": "unknown", "pb": "5.0.0"},
        "RG": [{"ID": "abcd1234/0--0", "PL": "PACBIO", "PU": movie, "SM": path_like_sample["sm"], "BC": path_like_sample["bcseq"], "PM": "REVIO"}],
    }


run_root = os.path.join(OUT, "run", "r84000_20250101_000000", "1_A01")
run_root2 = os.path.join(OUT, "run", "r84000_20250102_000000", "1_B01")
for rr in (run_root, run_root2):
    for sub in ("fail_reads", "hifi_reads"):
        os.makedirs(os.path.join(rr, sub), exist_ok=True)
aln_dir = os.path.join(OUT, "aln")
os.makedirs(aln_dir, exist_ok=True)

zmw = [1000]
BCSEQ = {"bc2001": "ACACACAGACTGTGAG", "bc2002": "ACACATCTCGTGAGAG", "bc2003": "ACACGCATGACACACT", "bc2004": "ACAGTCGAGCGCTGCG"}
unassigned_records = []

for sname, cfg in SAMPLES.items():
    rg = {"sm": f"Synthetic_{sname}", "bcseq": BCSEQ[cfg["bc"]]}
    movie = SINGLE.get(sname, MOVIE)
    hdr_un = write_header(rg, movie)
    hdr_al = dict(write_header(rg, movie))
    hdr_al["HD"] = {"VN": "1.6", "SO": "coordinate"}
    hdr_al["SQ"] = [{"SN": n, "LN": len(s)} for n, s in contigs.items()]
    hifi_raw, fail_raw, aligned = [], [], []
    for contig in [c for c in contigs if "_alt" not in c]:
        haps = {h: hap_contig(sname, contig, h) for h in (1, 2)}
        one_x = contig == "chrX" and cfg["sex"] == "M"
        rng = RX if contig == "chrX" else R
        target_bases = (15 if one_x else 30) * len(contigs[contig])
        made = 0
        while made < target_bases:
            hap = 1 if one_x else rng.choice((1, 2))
            hseq, alle = haps[hap]
            L = max(3000, int(rng.lognormvariate(9.65, 0.3)))  # ~16 kb
            s = rng.randrange(0, len(hseq) - L)
            e = s + L
            made += L
            mol = hseq[s:e]
            # fate
            expanded = alle is not None and (alle[1] - alle[0]) > alle[2] + 30 and s < alle[1] and e > alle[0]
            if expanded and cfg["exp_fate"]:
                ph, p8, p1 = cfg["exp_fate"]
            else:
                ph = 0.95 - 0.012 * max(0, L / 1000 - 10)
                p8 = (1 - ph) * 0.8
                p1 = (1 - ph) * 0.15
            u = rng.random()
            fate = "hifi" if u < ph else "f08" if u < ph + p8 else "f01" if u < ph + p8 + p1 else "gone"
            if fate == "gone":
                continue
            zmw[0] += 1
            minus = rng.random() < 0.5
            native = rc(mol) if minus else mol
            if fate == "hifi":
                seq_n = native
                name = f"{movie}/{zmw[0]}/ccs"
                hifi_raw.append((name, seq_n, None, 0.999, 12))
                # aligned record in reference orientation
                fwd = rc(seq_n) if minus else seq_n
                # coordinates on the reference
                if alle is None:
                    rs, cig = s, [(0, len(fwd))]
                else:
                    a0, a1, reflen = alle
                    delta = (a1 - a0) - reflen
                    if e <= a0:
                        rs, cig = s, [(0, len(fwd))]
                    elif s >= a1:
                        rs, cig = s - delta, [(0, len(fwd))]
                    elif s < a0 and e > a1:
                        left = a0 - s
                        right = e - a1
                        mid = len(fwd) - left - right
                        rs = s
                        if delta >= 0:
                            cig = [(0, left), (0, reflen), (1, mid - reflen), (0, right)] if mid > reflen else [(0, left), (0, max(mid, 0)), (2, reflen - max(mid, 0)), (0, right)]
                        else:
                            cig = [(0, left), (0, mid), (2, reflen - mid), (0, right)]
                    elif s < a0:  # ends inside allele -> soft clip
                        left = a0 - s
                        rs, cig = s, [(0, left), (4, len(fwd) - left)]
                    elif e > a1:  # starts inside allele
                        right = e - a1
                        rs, cig = a0 + reflen, [(4, len(fwd) - right), (0, right)]
                    else:
                        continue  # wholly inside repeat: unmapped, dropped from aligned BAM
                cig = [(op, n) for op, n in cig if n > 0]
                hp = hap if rng.random() < 0.9 and not one_x else None
                hw = cfg.get("hp_window")
                if hw and not (contig == hw[0] and hw[1] <= rs <= hw[2]):
                    hp = None
                aligned.append((name, fwd, contig, rs, cig, minus, hp))
            else:
                rate = 0.10 if fate == "f08" else 0.03
                ff = 8 if fate == "f08" else 1
                seq_n = mutate(native, rate, rng)
                name = f"{movie}/{zmw[0]}/{'0_' + str(len(seq_n)) if ff == 8 else 'ccs'}"
                fail_raw.append((name, seq_n, ff, -1.0 if ff == 8 else 0.97, 1 if ff == 8 else 5))

    def write_unaligned(path, recs, header):
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            for name, seq, ff, rq, np_ in recs:
                a = pysam.AlignedSegment(out.header)
                a.query_name = name
                a.query_sequence = seq
                a.flag = 4
                a.query_qualities = pysam.qualitystring_to_array("~" * len(seq))
                tags = [("rq", float(rq)), ("np", int(np_)), ("RG", "abcd1234/0--0")]
                if ff is not None:
                    tags.insert(0, ("ff", int(ff)))
                a.set_tags(tags)
                out.write(a)

    if sname in SINGLE:
        # single-plex cell: every 3rd fail read lost its barcode; plus a few expanded HiFi reads without barcode
        un_fail = fail_raw[::3]
        fail_raw = [r for i, r in enumerate(fail_raw) if i % 3]
        R2 = random.Random(13)
        hseq, alle = hap_contig(sname, "chrB", 1)
        un_hifi = []
        for j in range(3):
            s0 = alle[0] - R2.randrange(4000, 6000)
            e0 = alle[1] + R2.randrange(4000, 6000)
            seq = hseq[s0:e0]
            zmw[0] += 1
            un_hifi.append((f"{movie}/{zmw[0]}/ccs", rc(seq) if j % 2 else seq, None, 0.995, 8))
        rr = run_root2
        write_unaligned(os.path.join(rr, "fail_reads", f"{movie}.fail_reads.unassigned.bam"), un_fail, hdr_un)
        write_unaligned(os.path.join(rr, "hifi_reads", f"{movie}.hifi_reads.unassigned.bam"), un_hifi, hdr_un)
    else:
        rr = run_root
        unassigned_records.extend(fail_raw[:50])
    write_unaligned(os.path.join(rr, "fail_reads", f"{movie}.fail_reads.{cfg['bc']}.bam"), fail_raw, hdr_un)
    write_unaligned(os.path.join(rr, "hifi_reads", f"{movie}.hifi_reads.{cfg['bc']}.bam"), hifi_raw, hdr_un)

    tmp = os.path.join(aln_dir, f"{sname}.unsorted.bam")
    with pysam.AlignmentFile(tmp, "wb", header=hdr_al) as out:
        for name, seq, contig, rs, cig, minus, hp in aligned:
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            a.query_sequence = seq
            a.flag = 16 if minus else 0
            a.reference_id = out.get_tid(contig)
            a.reference_start = rs
            a.mapping_quality = 60
            a.cigartuples = cig
            a.query_qualities = pysam.qualitystring_to_array("~" * len(seq))
            tags = [("RG", "abcd1234/0--0")]
            if hp:
                tags.append(("HP", hp))
                tags.append(("PS", 1))
            a.set_tags(tags)
            out.write(a)
    final = os.path.join(aln_dir, f"{sname}.haplotagged.bam")
    pysam.sort("-o", final, tmp)
    pysam.index(final)
    os.remove(tmp)
    print(sname, "hifi", len(hifi_raw), "fail", len(fail_raw), "aligned", len(aligned))

# unassigned BAM (mixed reads) - must not be auto-included
write_unaligned(os.path.join(run_root, "fail_reads", f"{MOVIE}.fail_reads.unassigned.bam"), unassigned_records,
                write_header({"sm": "", "bcseq": ""}))

with open(os.path.join(OUT, "manifest.tsv"), "w") as f:
    f.write("sample_id\tgroup\tlocus\thifi_bam\tfail_bams\thifi_raw_bams\tnote\n")
    f.write(f"P1\tpositive\tTGAA\t{aln_dir}/P1.haplotagged.bam\tauto\tauto\tsynthetic biallelic GAA 600/800\n")
    f.write(f"N1\tnegative\t-\t{aln_dir}/N1.haplotagged.bam\t{run_root}/fail_reads/{MOVIE}.fail_reads.bc2002.bam\t-\tsynthetic normal\n")
    f.write(f"P2\tpositive\tTRFC\t{aln_dir}/P2.haplotagged.bam\tauto\t-\tsynthetic biallelic AAGGG 1000/900 mostly lost\n")
    f.write(f"P3\tpositive\tTGAA\t{aln_dir}/P3.haplotagged.bam\tauto\t-\tsynthetic one-allele GAA 700 (HP1), narrow phase block\n")
print("ok", OUT)
