#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From Snippy core VCF + cleaned alignment (with reference row kept) + Gubbins-masked alignment,
produce SNP-only MSAs:
  1) intersect_snps.fasta : Snippy core SNPs that survived cleaning (may include Ns after masking)
  2) strict_core_after_masking.fasta : only columns where all samples are A/C/G/T after masking
"""

import sys, argparse, gzip
from collections import OrderedDict

AMBIG = set("Nn?-RYSWKMBDHVryswkmbdhv")

def read_fasta(fn):
    opener = gzip.open if fn.endswith(".gz") else open
    name, seq = None, []
    with opener(fn, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def build_contig_offsets(ref_fa):
    offsets = OrderedDict()
    total = 0
    for cid, seq in read_fasta(ref_fa):
        offsets[cid] = (total, len(seq))
        total += len(seq)
    if not offsets:
        sys.exit("reference.fasta empty or unreadable")
    return offsets, total

def pick_ref_row(clean_aln, expected_len=None):
    # choose the sequence with the longest ungapped length as reference
    best = (None, "", -1)
    for name, seq in read_fasta(clean_aln):
        ung = len(seq.replace("-", ""))
        if ung > best[2]:
            best = (name, seq, ung)
    if best[0] is None:
        sys.exit("no sequences found in cleaned alignment")
    if expected_len and best[2] != expected_len:
        sys.stderr.write(f"[warn] ungapped length in cleaned ref row ({best[2]}) != reference length ({expected_len})\n")
    return best[0], best[1]  # (name, seq)

def build_maps_from_clean_ref(clean_ref_seq):
    """
    For each alignment column j, if ref has a base (not '-'), increment linear genome position.
    Returns:
      genomepos_to_col : dict {linear_pos0 -> aln_col_index}
      col_to_genomepos : list length = aln_len, with pos0 or None
    """
    genomepos_to_col = {}
    col_to_genomepos = [None] * len(clean_ref_seq)
    pos0 = 0
    for j, b in enumerate(clean_ref_seq):
        if b == "-":
            col_to_genomepos[j] = None
        else:
            col_to_genomepos[j] = pos0
            genomepos_to_col[pos0] = j
            pos0 += 1
    return genomepos_to_col, col_to_genomepos

def parse_vcf_core_positions(vcf_fn, contig_offsets):
    pos_set = set()
    opener = gzip.open if vcf_fn.endswith(".gz") else open
    with opener(vcf_fn, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.split("\t")
            if len(f) < 2: continue
            ctg = f[0]
            try:
                pos1 = int(f[1])
            except ValueError:
                continue
            if ctg not in contig_offsets:
                sys.stderr.write(f"[warn] VCF contig {ctg} not in reference; skipping\n")
                continue
            off, _ = contig_offsets[ctg]
            pos0 = off + (pos1 - 1)
            pos_set.add(pos0)
    return sorted(pos_set)

def extract_columns(msa_fn, keep_cols_sorted, out_fn_all, out_fn_strict):
    # read masked alignment (no reference row) and slice columns
    names, seqs = [], []
    for n, s in read_fasta(msa_fn):
        names.append(n)
        seqs.append(s)
    if not names:
        sys.exit("masked alignment has no sequences")

    # slice to kept columns
    sliced = []
    for s in seqs:
        # guard against length mismatches
        if len(s) < (keep_cols_sorted[-1] + 1):
            sys.exit("masked alignment shorter than requested column index; check you used the same cleaned alignment")
        sliced.append("".join(s[i] for i in keep_cols_sorted))

    # write intersect (may include Ns)
    with open(out_fn_all, "w") as out:
        for n, s in zip(names, sliced):
            out.write(f">{n}\n{s}\n")

    # strict core after masking: drop any column with any ambiguous/missing
    cols = list(zip(*sliced))
    ok_cols = [k for k, col in enumerate(cols) if not any(c in AMBIG or c == "-" for c in col)]
    strict = ["".join(s[i] for i in ok_cols) for s in sliced]

    with open(out_fn_strict, "w") as out:
        for n, s in zip(names, strict):
            out.write(f">{n}\n{s}\n")

def main():
    ap = argparse.ArgumentParser(description="Extract recombination-masked core SNPs from Gubbins using Snippy core VCF")
    ap.add_argument("--reference", required=True, help="reference.fasta used by Snippy (contigs in original order)")
    ap.add_argument("--clean_with_ref", required=True, help="cleaned full alignment WITH reference row kept (same cleaning as used for Gubbins)")
    ap.add_argument("--snippy_core_vcf", required=True, help="Snippy core VCF (core.vcf)")
    ap.add_argument("--gubbins_masked_aln", required=True, help="Gubbins masked full-genome alignment (cleaned, reference removed)")
    ap.add_argument("--out_prefix", default="core_after_gubbins", help="output prefix")
    args = ap.parse_args()

    contig_offsets, ref_len = build_contig_offsets(args.reference)
    ref_name, clean_ref_seq = pick_ref_row(args.clean_with_ref, expected_len=ref_len)
    genomepos_to_col, _ = build_maps_from_clean_ref(clean_ref_seq)

    # linearize VCF positions
    vcf_positions = parse_vcf_core_positions(args.snippy_core_vcf, contig_offsets)

    # keep only SNPs whose columns survived cleaning
    keep_cols = [genomepos_to_col[p] for p in vcf_positions if p in genomepos_to_col]
    if not keep_cols:
        sys.exit("no SNP columns survived cleaning; did you pass the correct files?")

    keep_cols_sorted = sorted(set(keep_cols))

    # extract from Gubbins-masked alignment
    extract_columns(
        args.gubbins_masked_aln,
        keep_cols_sorted,
        out_fn_all = f"{args.out_prefix}.intersect_snps.fasta",
        out_fn_strict = f"{args.out_prefix}.strict_core_after_masking.fasta"
    )

if __name__ == "__main__":
    main()

