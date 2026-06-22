from __future__ import annotations

import argparse
import csv

from .core import compute_ds_force


def _read_sequence(path: str, contig: str) -> str:
    current_name: str | None = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as fasta:
        for raw_line in fasta:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name == contig and chunks:
                    break
                current_name = line[1:].split()[0]
                chunks = []
            elif current_name == contig:
                chunks.append(line)
    if not chunks:
        raise ValueError(f"Contig {contig} not found, please try again with a different contig name.")
    return "".join(chunks)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute DS forces over sliding genome windows.")
    parser.add_argument("infile", help="Input FASTA file with the genome assembly.")
    parser.add_argument("contig", help="Contig to be considered for DS force computation.")
    parser.add_argument("slice_start", type=int, help="Start position, 1-based, on the contig sequence.")
    parser.add_argument("slice_end", type=int, help="End position, 1-based, on the contig sequence.")
    parser.add_argument("--outfile", default="", help="Name of the output CSV file.")
    parser.add_argument("--window_length", type=int, default=3000, help="Length of the sliding window.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.slice_start < 1:
        raise ValueError("'slice_start' parameter must be positive.")

    outfile = args.outfile or f"{args.contig}_{args.slice_start}-{args.slice_end}.csv"
    contig_seq = _read_sequence(args.infile, args.contig)
    seq = contig_seq[args.slice_start - 1 : min(args.slice_end + args.window_length - 1, len(contig_seq))]
    if len(seq) < args.window_length:
        open(outfile, "w", encoding="utf-8").close()
        return 0

    result = compute_ds_force(
        seq,
        return_lcs_positions=True,
        sliding_window_length=args.window_length,
    )
    ds_forces, range1s, range2s = result
    window_starts = [args.slice_start + i for i in range(len(seq) - args.window_length + 1)]

    with open(outfile, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["contig", "window start:end", "DS_force", "seqA start:end", "seqB start:end"],
        )
        writer.writeheader()
        for idx, ds_force in enumerate(ds_forces):
            window_start = window_starts[idx]
            window_end = window_start + args.window_length - 1
            seq_a = range1s[idx][0] + args.slice_start - 1, range1s[idx][1] + args.slice_start - 1
            seq_b = range2s[idx][0] + args.slice_start - 1, range2s[idx][1] + args.slice_start - 1
            writer.writerow(
                {
                    "contig": args.contig,
                    "window start:end": f"{window_start}:{window_end}",
                    "DS_force": round(ds_force, 6),
                    "seqA start:end": f"{seq_a[0]}:{seq_a[1]}",
                    "seqB start:end": f"{seq_b[0]}:{seq_b[1]}",
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
