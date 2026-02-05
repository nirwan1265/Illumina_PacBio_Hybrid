#!/usr/bin/env python3
"""Insert tandem repeat blocks into a genome FASTA.

Note: expects a single-contig FASTA. If multiple contigs are present,
pass a pre-merged FASTA or set --contig to select one record.
"""

import argparse
import random
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
                    # stop if already got one contig
                    if contig is None:
                        raise RuntimeError(
                            "FASTA has multiple records; provide single-contig FASTA or use --contig."
                        )
                    if hdr == contig:
                        name = hdr
                        seq_parts = []
                    else:
                        # ignore other contigs
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


def main():
    ap = argparse.ArgumentParser(
        description="Insert tandem repeats (segment duplicated k times) into a single-contig FASTA."
    )
    ap.add_argument("--in_fa", required=True)
    ap.add_argument("--out_fa", required=True)
    ap.add_argument("--n_events", type=int, default=10, help="number of repeat insertion events")
    ap.add_argument("--seg_len", type=int, default=1000, help="length of each duplicated segment (bp)")
    ap.add_argument("--copies", type=int, default=5, help="number of total copies in tandem (>=2)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contig", default=None, help="optional contig name if FASTA has multiple records")
    args = ap.parse_args()

    random.seed(args.seed)
    name, seq = read_fasta(args.in_fa, contig=args.contig)

    if args.copies < 2:
        raise ValueError("--copies must be >= 2")
    if args.seg_len <= 0 or args.seg_len > len(seq):
        raise ValueError("--seg_len must be between 1 and genome length")

    s = seq
    for _ in range(args.n_events):
        start = random.randint(0, len(s) - args.seg_len)
        seg = s[start : start + args.seg_len]
        insert_pos = random.randint(0, len(s))
        block = seg * args.copies
        s = s[:insert_pos] + block + s[insert_pos:]

    out_name = f"{name}|tandemrep_events{args.n_events}_len{args.seg_len}_copies{args.copies}"
    Path(args.out_fa).parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.out_fa, out_name, s)


if __name__ == "__main__":
    main()
