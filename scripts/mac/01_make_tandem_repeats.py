#!/usr/bin/env python3
"""Insert tandem repeat blocks into a genome FASTA.

Modes:
  - tandem: duplicate random genome segments in tandem
  - motif:  insert motif repeat blocks into random locations
  - both:   apply motif then tandem (same output FASTA)

Genome source:
  - default: use --in_fa reference
  - --random_genome: ignore --in_fa and generate random genome

Notes:
  - expects a single-contig FASTA. If multiple contigs are present,
    pass a pre-merged FASTA or set --contig to select one record.
  - --min_spacing applies per mode to reduce clustering of insertions.
  - --spacing_distribution can generate positions with Poisson/exponential gaps.
  - --gc_preserve attempts to keep GC% similar when inserting motif blocks
    by removing a GC-matched segment of equal length.
  - --motif_gc_target lets you generate a random motif with target GC.
"""

import argparse
import random
import math
from pathlib import Path


def read_fasta(path, contig=None):
    name = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                hdr = line[1:].split()[0]
                if name is None:
                    name = hdr
                else:
                    if contig is None:
                        raise RuntimeError(
                            "FASTA has multiple records; provide single-contig FASTA or use --contig."
                        )
                    if hdr == contig:
                        name = hdr
                        seq_parts = []
                    else:
                        continue
            else:
                if contig is None or name == contig:
                    seq_parts.append(line.upper())

    if name is None:
        raise RuntimeError("No FASTA header found.")
    if contig is not None and name != contig:
        raise RuntimeError(f"Contig '{contig}' not found in FASTA.")
    return name, "".join(seq_parts)


def write_fasta(path, name, seq, width=80):
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i : i + width] + "\n")


def gc_fraction(s):
    if not s:
        return 0.0
    s = s.upper()
    return (s.count("G") + s.count("C")) / len(s)


def pick_insert_positions_uniform(length, n, min_spacing):
    if n <= 0:
        return []
    if min_spacing <= 0:
        return [random.randint(0, length) for _ in range(n)]

    positions = []
    attempts = 0
    max_attempts = max(1000, n * 500)
    while len(positions) < n and attempts < max_attempts:
        attempts += 1
        pos = random.randint(0, length)
        if all(abs(pos - p) >= min_spacing for p in positions):
            positions.append(pos)
    if len(positions) < n:
        raise RuntimeError(
            f"Could not place {n} positions with min_spacing={min_spacing}. "
            "Try lowering min_spacing or n."
        )
    return positions


def pick_insert_positions_poisson(length, n, mean_gap, min_spacing):
    # simple 1D process: start at random position and add exponential gaps
    if n <= 0:
        return []
    if mean_gap <= 0:
        return pick_insert_positions_uniform(length, n, min_spacing)

    positions = []
    pos = random.randint(0, length)
    positions.append(pos)
    attempts = 0
    max_attempts = max(1000, n * 500)
    while len(positions) < n and attempts < max_attempts:
        attempts += 1
        gap = random.expovariate(1.0 / mean_gap)
        pos = int(pos + gap)
        if pos > length:
            # wrap around
            pos = pos % max(1, length)
        if all(abs(pos - p) >= min_spacing for p in positions):
            positions.append(pos)
    if len(positions) < n:
        raise RuntimeError(
            f"Could not place {n} positions with mean_gap={mean_gap}, min_spacing={min_spacing}. "
            "Try lowering min_spacing or n."
        )
    return positions


def remove_gc_matched_segment(seq, target_gc, seg_len, tol=0.05, max_tries=2000):
    if seg_len <= 0 or seg_len > len(seq):
        return seq
    for _ in range(max_tries):
        start = random.randint(0, len(seq) - seg_len)
        seg = seq[start : start + seg_len]
        if abs(gc_fraction(seg) - target_gc) <= tol:
            return seq[:start] + seq[start + seg_len :]
    start = random.randint(0, len(seq) - seg_len)
    return seq[:start] + seq[start + seg_len :]


def make_random_genome(length, gc_target):
    if length <= 0:
        raise ValueError("--random_length must be > 0")
    if not (0.0 <= gc_target <= 1.0):
        raise ValueError("--random_gc must be between 0 and 1")
    seq = []
    for _ in range(length):
        if random.random() < gc_target:
            seq.append(random.choice(["G", "C"]))
        else:
            seq.append(random.choice(["A", "T"]))
    return "".join(seq)


def make_gc_matched_motif(length, gc_target):
    if length <= 0:
        raise ValueError("--motif_len must be > 0 for GC-matched motif generation")
    if not (0.0 <= gc_target <= 1.0):
        raise ValueError("--motif_gc_target must be between 0 and 1")
    motif = []
    for _ in range(length):
        if random.random() < gc_target:
            motif.append(random.choice(["G", "C"]))
        else:
            motif.append(random.choice(["A", "T"]))
    return "".join(motif)


def main():
    ap = argparse.ArgumentParser(
        description="Insert tandem repeats and/or motif repeats into a single-contig FASTA."
    )
    ap.add_argument("--in_fa", required=False, default=None)
    ap.add_argument("--out_fa", required=True)
    ap.add_argument(
        "--mode",
        choices=["tandem", "motif", "both"],
        default="tandem",
        help="repeat mode to apply (default: tandem)",
    )
    ap.add_argument("--n_events", type=int, default=10, help="tandem: number of insertion events")
    ap.add_argument("--seg_len", type=int, default=1000, help="tandem: duplicated segment length (bp)")
    ap.add_argument("--copies", type=int, default=5, help="tandem: total copies in tandem (>=2)")
    ap.add_argument("--motif", default="ATTA", help="motif: repeat string, e.g. ATTA")
    ap.add_argument("--motif_repeat", type=int, default=1000, help="motif: repeats per block")
    ap.add_argument("--motif_events", type=int, default=5, help="motif: number of blocks inserted")
    ap.add_argument(
        "--motif_mode",
        choices=["insert", "replace"],
        default="insert",
        help="motif: insert (length grows) or replace (length preserved)",
    )
    ap.add_argument("--min_spacing", type=int, default=0, help="minimum spacing between insertion positions")
    ap.add_argument(
        "--spacing_distribution",
        choices=["uniform", "poisson"],
        default="uniform",
        help="how to sample insertion positions (default: uniform)",
    )
    ap.add_argument(
        "--spacing_mean",
        type=float,
        default=10000.0,
        help="mean gap for poisson spacing (bp)",
    )
    ap.add_argument(
        "--gc_preserve",
        action="store_true",
        help="attempt to preserve GC when inserting motifs by removing a GC-matched segment",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contig", default=None, help="optional contig name if FASTA has multiple records")

    # random genome options
    ap.add_argument("--random_genome", action="store_true", help="generate random genome instead of using --in_fa")
    ap.add_argument("--random_length", type=int, default=100000, help="random genome length (bp)")
    ap.add_argument("--random_gc", type=float, default=0.5, help="random genome GC fraction (0-1)")

    # GC-matched motif generation
    ap.add_argument("--motif_gc_target", type=float, default=None, help="generate a random motif with this GC")
    ap.add_argument("--motif_len", type=int, default=12, help="length of generated motif if using --motif_gc_target")

    args = ap.parse_args()

    random.seed(args.seed)

    if args.random_genome:
        name = "RandomGenome"
        seq = make_random_genome(args.random_length, args.random_gc)
    else:
        if not args.in_fa:
            raise RuntimeError("--in_fa is required unless --random_genome is set")
        name, seq = read_fasta(args.in_fa, contig=args.contig)

    if args.mode in ("tandem", "both"):
        if args.copies < 2:
            raise ValueError("--copies must be >= 2")
        if args.seg_len <= 0 or args.seg_len > len(seq):
            raise ValueError("--seg_len must be between 1 and genome length")

    if args.mode in ("motif", "both"):
        if args.motif_gc_target is not None:
            args.motif = make_gc_matched_motif(args.motif_len, args.motif_gc_target)
        if not args.motif or any(b not in "ACGTacgt" for b in args.motif):
            raise ValueError("--motif must be a DNA string (A/C/G/T)")
        if args.motif_repeat < 1 or args.motif_events < 1:
            raise ValueError("--motif_repeat and --motif_events must be >= 1")

    def get_positions(length, n):
        if args.spacing_distribution == "poisson":
            return pick_insert_positions_poisson(length, n, args.spacing_mean, args.min_spacing)
        return pick_insert_positions_uniform(length, n, args.min_spacing)

    s = seq
    if args.mode in ("motif", "both"):
        motif_block = args.motif.upper() * args.motif_repeat
        positions = get_positions(len(s), args.motif_events)
        for insert_pos in sorted(positions, reverse=True):
            if args.motif_mode == "replace":
                s = s[:insert_pos] + motif_block + s[insert_pos + len(motif_block) :]
            else:
                s = s[:insert_pos] + motif_block + s[insert_pos:]
                if args.gc_preserve:
                    s = remove_gc_matched_segment(s, gc_fraction(motif_block), len(motif_block))

    if args.mode in ("tandem", "both"):
        positions = get_positions(len(s), args.n_events)
        for insert_pos in sorted(positions, reverse=True):
            start = random.randint(0, len(s) - args.seg_len)
            seg = s[start : start + args.seg_len]
            block = seg * args.copies
            s = s[:insert_pos] + block + s[insert_pos:]

    out_name = (
        f"{name}|mode{args.mode}"
        f"|tandem_events{args.n_events}_len{args.seg_len}_copies{args.copies}"
        f"|motif_{args.motif}_rep{args.motif_repeat}_events{args.motif_events}"
        f"|motifmode{args.motif_mode}_gc_preserve{int(args.gc_preserve)}"
        f"|spacing_{args.spacing_distribution}_mean{args.spacing_mean}_min{args.min_spacing}"
    )
    if args.random_genome:
        out_name += f"|random_len{args.random_length}_gc{args.random_gc}"
    if args.motif_gc_target is not None:
        out_name += f"|motif_gc_target{args.motif_gc_target}_len{args.motif_len}"

    Path(args.out_fa).parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.out_fa, out_name, s)


if __name__ == "__main__":
    main()
