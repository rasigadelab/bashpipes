#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Title: extract_core_snps_after_gubbins.py
# Description: Produces SNP-only Multi-Sequence Alignments (MSAs) based on Snippy core VCF, cleaned alignment and Gubbins masked alignment.
# Author: Aurélie Fischer
# Research team: PHE3ID, CIRI, Hospices Civils de Lyon
# License: GNU Affero General Public License v3.0 (AGPL-3.0)
# SPDX-License-Identifier: AGPL-3.0-only
#
# Copyright (C) 2026 Aurélie Fischer


# From Snippy core VCF + cleaned alignment (with reference row kept) + Gubbins-masked alignment,
# produce SNP-only MSAs:
#   1) intersect_snps.fasta : Snippy core SNPs that survived cleaning (may include Ns after masking)
#   2) strict_core_after_masking.fasta : only columns where all samples are A/C/G/T after masking

import sys, argparse, gzip
from collections import OrderedDict

AMBIG = set("Nn?-RYSWKMBDHVryswkmbdhv")

def read_fasta(fn):
    """
    Read FASTA file.

    Args:
        fn (str): FASTA filename.

    Returns:
        None
    """
    # Opening mechanism
    opener = gzip.open if fn.endswith(".gz") else open
    name, seq = None, []
    with opener(fn, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n\r") # Remove any linebreak
            if not line: continue 
            # First, get sequence ID
            if line.startswith(">"):
                if name is not None: 
                    # Save previous name + sequence
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            # Then, get DNA sequence
            else:
                seq.append(line)
    if name is not None:
        # Finish reading by saving last name + sequence
        yield name, "".join(seq)

def build_contig_offsets(ref_fa):
    """
    Count cumulative contigs length in contig appearance in FASTA order.

    Args:
        ref_fa (str): FASTA filename.

    Returns:
        Ordered dictionary of {contig: (cumulated contig length, contig length)}
        Total contigs length
    """
    offsets = OrderedDict()
    total = 0
    # Get contig name (key) and cumulated contigs length + length of current contig sequence (tuple) in an ordered dict
    for cid, seq in read_fasta(ref_fa):
        offsets[cid] = (total, len(seq))
        total += len(seq)
    if not offsets:
        sys.exit("reference.fasta empty or unreadable")
    return offsets, total

def pick_ref_row(clean_aln, expected_len=None):
    """
    Choose the sequence with the longest ungapped length as reference.

    Args:
        clean_aln (str): FASTA clean alignment filename.
        expected_length (int): Expected best ungapped length in bp

    Returns:
        Best ungapped contig name
        Best ungapped contig sequence
    """
    best = (None, "", -1)
    for name, seq in read_fasta(clean_aln):
        # Get contig length without gaps ("-")
        ung = len(seq.replace("-", ""))
        # If contig ungapped length is longer than already stored one
        if ung > best[2]:
            # Then replace best with current contig info
            best = (name, seq, ung)
    if best[0] is None:
        sys.exit("no sequences found in cleaned alignment")
    if expected_len and best[2] != expected_len:
        sys.stderr.write(f"[warn] ungapped length in cleaned ref row ({best[2]}) != reference length ({expected_len})\n")
    return best[0], best[1]  # (name, seq)

def build_maps_from_clean_ref(clean_ref_seq):
    """
    For each alignment column idx, if reference has a base (not '-'), increment linear genome position.

    Args:
        clean_ref_seq (str): clean DNA sequence of reference genome extracted from alignment file.

    Returns:
      genomepos_to_col : dictionary indicating the position of informative bases (no "-")
      col_to_genomepos : list length = aln_len, with pos0 or None
    """
    genomepos_to_col = {} # {linear_pos0 -> ref_col_index}
    col_to_genomepos = [None] * len(clean_ref_seq) # (list length = aln_len, with pos0 or None)
    pos0 = 0 # position of next informative base
    for idx, base in enumerate(clean_ref_seq):
        if base == "-":
            # At this position in reference there is no base (-).
            # pos0 does not change
            col_to_genomepos[idx] = None
        else:
            # At this position in reference there is a base (ATCG).
            # Info: at this position, there is an informative base 
            col_to_genomepos[idx] = pos0
            # Info: informative base is at idx position in reference
            genomepos_to_col[pos0] = idx
            pos0 += 1
    return genomepos_to_col, col_to_genomepos

def parse_vcf_core_positions(vcf_fn, contig_offsets):
    """
    Transform positions of core variants from position on contig to position on reference genome.

    Args:
        vcf_fn (str): VCF file containing positions of all core variants called
        contig_offsets (dict): cumulative length of contigs and contig length of reference genome per contig

    Returns:
        sorted(pos_set) : sorted set of core variant positions (reference-based)
    """
    pos_set = set()
    opener = gzip.open if vcf_fn.endswith(".gz") else open
    # Read VCF line / line
    with opener(vcf_fn, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue # avoid comments
            f = line.split("\t")
            if len(f) < 2: continue # avoid uninformative lines
            ctg = f[0] # contig name "NODE_1_length_400000_cov_45.202521_pilon"
            try:
                pos1 = int(f[1]) # variant position on reference
            except ValueError:
                continue
            if ctg not in contig_offsets:
                sys.stderr.write(f"[warn] VCF contig {ctg} not in reference; skipping\n")
                continue
            off, _ = contig_offsets[ctg] # Get start and length of contig in reference fasta
            pos0 = off + (pos1 - 1) # Convert position on contig to position on whole assembly
            pos_set.add(pos0) # Add to set of variant positions
    return sorted(pos_set)

def extract_columns(msa_fn, keep_cols_sorted, out_fn_all, out_fn_strict):
    """
    Filter core sites in MSA file based on core positions and write as output the new alignment file with and without masking ambiguous sites.

    Args:
        msa_fn (str): a Multi Sequence Alignment fasta file
        keep_cols_sorted (list): List of positions to keep in MSA (ex: core SNPs positions) 
        out_fn_all (str): Output file for filtered alignment
        out_fn_strict (str): Output file for filtered alignment without ambiguous sites.

    Returns:
        None.
    """
    # Read masked alignment (no reference row)
    names, seqs = [], []
    for n, s in read_fasta(msa_fn):
        names.append(n)
        seqs.append(s)
    if not names:
        sys.exit("Masked alignment has no sequences.")

    # Slice every sequence of input alignment file to kept columns only
    sliced = []
    for s in seqs:
        # Guard against length mismatches
        if len(s) < (keep_cols_sorted[-1] + 1):
            sys.exit("Masked alignment shorter than requested column index; check you used the same cleaned alignment.")
        sliced.append("".join(s[i] for i in keep_cols_sorted))

    # Re-write new alignment file with only kept positions (may include Ns)
    with open(out_fn_all, "w") as out:
        for n, s in zip(names, sliced):
            out.write(f">{n}\n{s}\n")

    # Strict core after masking: drop any column with any ambiguous/missing
    cols = list(zip(*sliced))
    ok_cols = [k for k, col in enumerate(cols) if not any(c in AMBIG or c == "-" for c in col)]
    strict = ["".join(s[i] for i in ok_cols) for s in sliced]

    # Write strict core alignment after masking (without any "-" or "N")
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

    # Step1- Get contigs start all along reference genome
    contig_offsets, ref_len = build_contig_offsets(args.reference)
    # Step2- Identify reference entry in core alignment file
    ref_name, clean_ref_seq = pick_ref_row(args.clean_with_ref, expected_len=ref_len)
    # Step3- Get informative base sites positions of reference in alignment file
    genomepos_to_col, _ = build_maps_from_clean_ref(clean_ref_seq)
    # Step4- Linearize VCF positions (from contig-based to reference-based positions)
    vcf_positions = parse_vcf_core_positions(args.snippy_core_vcf, contig_offsets)
    # Step5- Keep only SNPs whose columns survived cleaning
    keep_cols = [genomepos_to_col[p] for p in vcf_positions if p in genomepos_to_col]
    if not keep_cols:
        sys.exit("No SNP columns survived cleaning; did you pass the correct files?")
    keep_cols_sorted = sorted(set(keep_cols))
    # Step6- Extract from Gubbins-masked alignment
    extract_columns(
        args.gubbins_masked_aln,
        keep_cols_sorted,
        out_fn_all = f"{args.out_prefix}.intersect_snps.fasta",
        out_fn_strict = f"{args.out_prefix}.strict_core_after_masking.fasta"
    )

if __name__ == "__main__":
    main()

