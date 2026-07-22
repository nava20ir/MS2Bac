#!/usr/bin/env python3

from pathlib import Path
import math


def read_fasta(fasta_file):
    records = []
    header = None
    seq = []

    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq)))

    return records


def split_fasta(fasta_file, num_splits, output_prefix="chunk"):
    records = read_fasta(fasta_file)

    if num_splits <= 0:
        raise ValueError("num_splits must be > 0")

    total = len(records)
    chunk_size = math.ceil(total / num_splits)

    for i in range(num_splits):
        start = i * chunk_size
        end = min(start + chunk_size, total)

        if start >= total:
            break

        out_file = f"{output_prefix}_{i+1}.fasta"

        with open(out_file, "w") as out:
            for header, seq in records[start:end]:
                out.write(f"{header}\n{seq}\n")

        print(f"Wrote {end-start} sequences to {out_file}")


if __name__ == "__main__":
    import sys
    split_fasta(
        fasta_file=sys.argv[1],
        num_splits=int(sys.argv[2]),
        output_prefix="peptides"
    )
