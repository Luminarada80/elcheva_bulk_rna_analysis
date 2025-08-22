#!/usr/bin/env python3
"""
Reads mapping summary: counts reads in exon, intron, and long-junction categories.

Input:
  - BAM (indexed alongside)
  - BED12 gene annotation (exon blocks derived from blockSizes/blockStarts)
Output:
  - Text file with summary counts

Original author: Long Gao
Modernized for Python 3 / pysam
"""

import sys
import argparse
import time
import gc
from collections import defaultdict

import pysam

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam", help="Input BAM (indexed)")
    ap.add_argument("bed", help="Annotation in BED12 format")
    ap.add_argument("output", help="Output text file")
    ap.add_argument("--division", type=int, default=100,
                    help="Number of bins for speed index (default: 100)")
    ap.add_argument("--skip-chr", nargs="*", default=["random", "chrUn", "chrM"],
                    help="Skip chromosomes containing any of these substrings")
    ap.add_argument("--junction-min-skip", type=int, default=200,
                    help="Minimum total 'N' skip length to count as junction (default: 200)")
    return ap.parse_args()

def bed12_to_blocks(line):
    """
    Parse a BED12 line into exon and intron blocks.
    Returns (chrom, exon_blocks, intron_blocks) where blocks are [(start, end), ...]
    """
    # BED12 columns (0-based): chrom, start, end, name, score, strand, thickStart, thickEnd,
    # itemRgb, blockCount, blockSizes, blockStarts
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise ValueError("BED must be BED12 with at least 12 columns.")

    chrom = fields[0]
    tx_start = int(fields[1])
    # tx_end   = int(fields[2])  # not needed below

    block_sizes = [int(x) for x in fields[10].strip(",").split(",") if x]
    block_starts = [int(x) for x in fields[11].strip(",").split(",") if x]

    if len(block_sizes) != len(block_starts):
        raise ValueError("blockSizes and blockStarts length mismatch in BED line:\n" + line)

    exon_blocks = []
    for bs, rel_start in zip(block_sizes, block_starts):
        exon_blocks.append((tx_start + rel_start, tx_start + rel_start + bs))

    intron_blocks = []
    for i in range(len(exon_blocks) - 1):
        intron_blocks.append((exon_blocks[i][1], exon_blocks[i+1][0]))

    return chrom, exon_blocks, intron_blocks

def should_skip_chrom(chrom, skip_terms):
    for term in skip_terms:
        if term in chrom:
            return True
    return False

def build_speed_index(blocks_by_chr, division):
    """
    Build a speed/anchor index to narrow binary-like searches.
    For each chromosome, returns list of tuples:
        (first_block_start_of_bin, start_index, end_index)
    Handles small lists gracefully.
    """
    speed = {}
    for chrom, blocks in blocks_by_chr.items():
        # blocks are already sorted by (start, end)
        n = len(blocks)
        if n == 0:
            speed[chrom] = [(0, 0, 0)]
            continue

        # number of bins cannot exceed n; ensure at least size 1 bins
        bins = min(division, n)
        # compute roughly equal-sized bins
        base = n // bins
        rem = n % bins

        index_list = []
        start = 0
        for i in range(bins):
            size = base + (1 if i < rem else 0)
            end = start + size
            # guard (shouldn't happen if size>0)
            if start >= n:
                start = n - 1
            first_start = blocks[start][0]
            index_list.append((first_start, start, end))
            start = end

        # just in case, ensure last bin ends at n
        if index_list and index_list[-1][2] != n:
            first_start = blocks[index_list[-1][1]][0]
            index_list[-1] = (first_start, index_list[-1][1], n)

        speed[chrom] = index_list

    return speed

def pick_sublist(chrom, pos, speed_index):
    """
    Given a chrom, a genomic coordinate 'pos', and that chrom's speed index,
    return (start_idx, end_idx) for the block sublist to scan.
    """
    anchors = speed_index.get(chrom)
    if not anchors:
        return (0, 0)

    # linear scan over anchors (anchors length <= division, typically 100)
    for i in range(len(anchors) - 1):
        left_anchor = anchors[i][0]
        right_anchor = anchors[i + 1][0]
        if left_anchor <= pos <= right_anchor:
            return anchors[i][1], anchors[i][2]
    # if not found, use last bin
    return anchors[-1][1], anchors[-1][2]

def read_overlaps_blocks(read_left, read_right, blocks):
    """
    Check if [read_left, read_right) overlaps any interval in 'blocks' (sorted list of (start,end)).
    Early-exit once block.start > read_right.
    """
    for (left, right) in blocks:
        # three overlap cases: left edge inside, fully inside, right edge inside
        if (read_left <= left <= read_right) or \
           (left <= read_left and read_right <= right) or \
           (read_left <= right <= read_right):
            return True
        if left > read_right:
            break
    return False

def sum_cigar_N(op_tuples):
    """
    Sum lengths of REF_SKIP ('N') ops in a pysam cigar tuple list.
    CIGAR op codes: 3 == N (ref skip)
    """
    total = 0
    if not op_tuples:
        return 0
    for op, length in op_tuples:
        if op == 3:
            total += length
    return total

def main():
    args = parse_args()
    t1 = time.time()

    # Build exon/intron interval indexes from BED12
    exon_index = defaultdict(list)   # chrom -> [(start,end), ...]
    intron_index = defaultdict(list)

    with open(args.bed, "r") as bedf:
        for line in bedf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, exons, introns = bed12_to_blocks(line)
            if should_skip_chrom(chrom, args.skip_chr):
                continue
            # accumulate unique blocks (optional; most BEDs already have unique blocks per transcript)
            for blk in exons:
                exon_index[chrom].append(blk)
            for blk in introns:
                intron_index[chrom].append(blk)

    # deduplicate & sort
    for chrom in list(exon_index.keys()):
        exon_index[chrom] = sorted(set(exon_index[chrom]))
    for chrom in list(intron_index.keys()):
        intron_index[chrom] = sorted(set(intron_index[chrom]))

    # build speed indices
    exon_speed = build_speed_index(exon_index, args.division)
    intron_speed = build_speed_index(intron_index, args.division)

    # free large temporary structures from original script (no-op here but keep parity)
    gc.collect()

    # Open BAM
    bam = pysam.AlignmentFile(args.bam, "rb")

    exon_num = 0
    intron_num = 0
    junction_num = 0
    all_num = 0

    # Iterate reads
    for read in bam.fetch(until_eof=False):
        # Skip unmapped/secondary/supplementary
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        ref_id = read.reference_id
        if ref_id < 0:
            continue  # safety
        chrom = bam.get_reference_name(ref_id)
        if should_skip_chrom(chrom, args.skip_chr):
            continue

        all_num += 1

        # Junction: sum 'N' operations in CIGAR
        n_skip = sum_cigar_N(read.cigartuples)
        if n_skip >= args.junction_min_skip:
            junction_num += 1

        # Coordinates for overlap test
        read_left = read.reference_start  # 0-based left-most ref coordinate
        # Use query-aligned length; if None, fall back to aligned_length
        qlen = read.query_alignment_length if read.query_alignment_length is not None else read.query_length
        if qlen is None:
            # last resort: compute from CIGAR (M,=,X,I,D contribute differently; but for a window,
            # we mainly need a span—use reference span)
            read_right = read.reference_end
        else:
            # convert to a half-open end; for overlap logic we use reference end instead (safer)
            read_right = read.reference_end

        # If we lack a reference end for some reason, fallback to left + qlen
        if read_right is None:
            read_right = read_left + (qlen or 0)

        # Exon / intron classification
        at_exon = False
        at_intron = False

        if chrom in exon_index and exon_index[chrom]:
            s, e = pick_sublist(chrom, read_left, exon_speed)
            if s < e:
                at_exon = read_overlaps_blocks(read_left, read_right, exon_index[chrom][s:e])

        if not at_exon and chrom in intron_index and intron_index[chrom]:
            s, e = pick_sublist(chrom, read_left, intron_speed)
            if s < e:
                at_intron = read_overlaps_blocks(read_left, read_right, intron_index[chrom][s:e])

        if at_exon:
            exon_num += 1
        elif at_intron:
            intron_num += 1
        # else: intergenic (counted later)

    bam.close()
    t2 = time.time()

    intergenic = all_num - exon_num - intron_num

    with open(args.output, "w") as out:
        out.write(f"It took {t2 - t1:.2f} sec\n")
        out.write(f"Exon Reads : {exon_num}\n")
        out.write(f"Intron Reads : {intron_num}\n")
        out.write(f"Junction Reads : {junction_num}\n")
        out.write(f"Intergenic Reads : {intergenic}\n")
        out.write(f"Number of uniquely mapped reads : {all_num}\n")


if __name__ == "__main__":
    sys.exit(main())
