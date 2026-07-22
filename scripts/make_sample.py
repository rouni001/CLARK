#!/usr/bin/env python3
"""Sample reads from the genomes configured by scripts/set_targets.sh."""

import argparse
import gzip
import random
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Sample reads from CLARK target genomes.")
    parser.add_argument("--targets", required=True, help="targets.txt written by set_targets.sh")
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=None)
    return parser.parse_args()


def read_targets(path):
    targets = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            targets.append((fields[0], fields[1]))
    return targets


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def load_sequence(path):
    """Concatenate every record in a FASTA file into one sequence string."""
    chunks = []
    try:
        with open_maybe_gzip(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    continue
                chunks.append(line.strip())
    except OSError as exc:
        print("Warning: failed to read %s: %s" % (path, exc), file=sys.stderr)
        return ""
    return "".join(chunks)


def main():
    args = parse_args()
    if args.count <= 0:
        raise SystemExit("--count must be a positive integer")
    if args.length <= 0:
        raise SystemExit("--length must be a positive integer")

    rng = random.Random(args.seed)

    targets = read_targets(args.targets)
    if not targets:
        raise SystemExit("no targets found in %s" % args.targets)

    usable = []
    for path, taxid in targets:
        seq = load_sequence(path)
        if len(seq) >= args.length:
            usable.append((path, taxid, seq))
        else:
            print("Warning: skipping %s (shorter than requested read length)" % path, file=sys.stderr)

    if not usable:
        raise SystemExit("no target genome is long enough for a %d-base read" % args.length)

    with open(args.output, "w", encoding="utf-8") as out:
        for i in range(1, args.count + 1):
            path, taxid, seq = rng.choice(usable)
            start = rng.randint(0, len(seq) - args.length)
            read = seq[start:start + args.length]
            name = path.rsplit("/", 1)[-1]
            out.write(">sample_%d|source=%s|taxid=%s|pos=%d-%d\n" % (i, name, taxid, start, start + args.length))
            out.write(read + "\n")

    print("wrote %d read(s) to %s" % (args.count, args.output), file=sys.stderr)


if __name__ == "__main__":
    main()
