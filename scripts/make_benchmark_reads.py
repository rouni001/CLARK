#!/usr/bin/env python3
"""Simulate classification-benchmark reads with published, citable specs.

Ad-hoc uniform-random sampling (scripts/make_sample.py) is fine for smoke
tests, but a performance/accuracy paper needs an input whose *shape* traces
back to a dataset other papers already used. CLARK's own BMC Genomics 2015
paper benchmarked against three read sets from the Kraken project ("HiSeq",
"MiSeq", "simBA-5") and one from the FAMeS project ("simHC"); the ISCA 2021
"Sieve" paper reuses the same lineage. The exact original FASTQ/FASTA files
are hosted on decade-old project pages (ccb.jhu.edu, fames.jgi-psf.org) that
may no longer serve the literal bytes, so this script instead simulates
reads matching each dataset's *published characteristics* -- read length,
source-genome count, and relative error rate -- from whatever reference
genomes scripts/set_targets.sh already configured. This is a documented
approximation, not a re-download of the original files; say so in any
paper that cites it.

Profile specs (see README/paper notes for citations):
  hiseq   92 bp reads,  equal proportion from 10 genomes, baseline error rate
  miseq  156 bp reads,  equal proportion from 10 genomes, baseline error rate
  simba5 100 bp reads,  broad draw across all available genomes, 5x error rate
  simhc  800 bp reads,  power-law abundance across up to 113 genomes, elevated error rate

Every emitted read's header carries its true taxid, and a companion TSV
ground-truth file maps read_id -> true taxid for scripts/eval_accuracy.py.
"""

import argparse
import gzip
import random
import sys

BASES = "ACGT"

PROFILES = {
    "hiseq": dict(read_length=92, n_genomes=10, count=10000, error_rate=0.001, abundance="equal"),
    "miseq": dict(read_length=156, n_genomes=10, count=10000, error_rate=0.001, abundance="equal"),
    "simba5": dict(read_length=100, n_genomes=None, count=10000, error_rate=0.005, abundance="broad"),
    "simhc": dict(read_length=800, n_genomes=113, count=10000, error_rate=0.01, abundance="powerlaw"),
    "custom": dict(read_length=None, n_genomes=None, count=None, error_rate=None, abundance="equal"),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", required=True, help="targets.txt written by set_targets.sh")
    parser.add_argument("--profile", required=True, choices=sorted(PROFILES.keys()))
    parser.add_argument("--output", required=True, help="output FASTA path")
    parser.add_argument("--ground-truth", required=True, help="output TSV path (read_id, true_taxid, source_file, profile)")
    parser.add_argument("--count", type=int, default=None, help="override the profile's default read count")
    parser.add_argument("--length", type=int, default=None, help="override the profile's default read length")
    parser.add_argument("--error-rate", type=float, default=None, help="override the profile's default per-base substitution rate")
    parser.add_argument("--n-genomes", type=int, default=None, help="override the profile's default number of source genomes")
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


def mutate(read, error_rate, rng):
    if error_rate <= 0:
        return read
    bases = list(read)
    for i, base in enumerate(bases):
        if base not in "ACGT":
            continue
        if rng.random() < error_rate:
            choices = [b for b in BASES if b != base]
            bases[i] = rng.choice(choices)
    return "".join(bases)


def select_genomes(usable, n_genomes, rng):
    if n_genomes is None or n_genomes >= len(usable):
        return list(usable)
    return rng.sample(usable, n_genomes)


def build_weights(genomes, abundance, rng):
    n = len(genomes)
    if abundance == "equal":
        return [1.0] * n
    if abundance == "broad":
        # Draw broadly across every available genome, roughly proportional
        # to genome length (a longer reference contributes more distinct
        # source positions), approximating simBA-5's wide 1,967-taxa draw.
        return [max(len(seq), 1) for _, _, seq in genomes]
    if abundance == "powerlaw":
        # Steep power-law-like abundance, approximating FAMeS simHC.
        order = list(range(n))
        rng.shuffle(order)
        weights = [0.0] * n
        for rank, idx in enumerate(order, start=1):
            weights[idx] = 1.0 / (rank ** 1.5)
        return weights
    raise ValueError("unknown abundance model: %s" % abundance)


def weighted_choice(genomes, weights, total, rng):
    x = rng.random() * total
    acc = 0.0
    for genome, w in zip(genomes, weights):
        acc += w
        if x <= acc:
            return genome
    return genomes[-1]


def main():
    args = parse_args()
    spec = dict(PROFILES[args.profile])

    if args.profile == "custom":
        if args.count is None or args.length is None or args.error_rate is None:
            raise SystemExit("--profile custom requires --count, --length, and --error-rate")

    count = args.count if args.count is not None else spec["count"]
    length = args.length if args.length is not None else spec["read_length"]
    error_rate = args.error_rate if args.error_rate is not None else spec["error_rate"]
    n_genomes = args.n_genomes if args.n_genomes is not None else spec["n_genomes"]

    if count is None or count <= 0:
        raise SystemExit("--count must be a positive integer")
    if length is None or length <= 0:
        raise SystemExit("--length must be a positive integer")
    if error_rate is None or error_rate < 0:
        raise SystemExit("--error-rate must be zero or a positive fraction")

    rng = random.Random(args.seed)

    targets = read_targets(args.targets)
    if not targets:
        raise SystemExit("no targets found in %s" % args.targets)

    usable = []
    for path, taxid in targets:
        seq = load_sequence(path)
        if len(seq) >= length:
            usable.append((path, taxid, seq))
        else:
            print("Warning: skipping %s (shorter than requested read length)" % path, file=sys.stderr)

    if not usable:
        raise SystemExit("no target genome is long enough for a %d-base read" % length)

    genomes = select_genomes(usable, n_genomes, rng)
    weights = build_weights(genomes, spec["abundance"], rng)
    total_weight = sum(weights)

    with open(args.output, "w", encoding="utf-8") as out, open(args.ground_truth, "w", encoding="utf-8") as gt:
        gt.write("read_id\ttrue_taxid\tsource_file\tprofile\n")
        for i in range(1, count + 1):
            path, taxid, seq = weighted_choice(genomes, weights, total_weight, rng)
            start = rng.randint(0, len(seq) - length)
            read = mutate(seq[start:start + length], error_rate, rng)
            name = path.rsplit("/", 1)[-1]
            read_id = "%s_%d" % (args.profile, i)
            out.write(">%s|source=%s|taxid=%s|pos=%d-%d|profile=%s\n" % (read_id, name, taxid, start, start + length, args.profile))
            out.write(read + "\n")
            gt.write("%s\t%s\t%s\t%s\n" % (read_id, taxid, name, args.profile))

    print(
        "wrote %d %s-profile read(s) (length=%d, error_rate=%.4f, genomes=%d) to %s (ground truth: %s)"
        % (count, args.profile, length, error_rate, len(genomes), args.output, args.ground_truth),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
