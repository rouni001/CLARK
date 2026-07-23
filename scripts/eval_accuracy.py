#!/usr/bin/env python3
"""Compute classification accuracy (sensitivity/precision) for a CLARK run.

Compares a CLARK results CSV (header "Object_ID, Length, Assignment", one
row per classified read) against the ground-truth TSV written alongside a
scripts/make_benchmark_reads.py FASTA (read_id, true_taxid, source_file,
profile). Reports both a micro-averaged (pooled) summary and a per-taxid
breakdown, so a paper can cite sensitivity/precision rather than only
throughput.

Reads CLARK leaves unclassified (Assignment is "NA" or "-1") count as a
false negative for their true taxid, never as a false positive.
"""

import argparse
import csv
import sys

UNCLASSIFIED = {"NA", "-1", ""}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, help="CLARK results CSV")
    parser.add_argument("--ground-truth", required=True, help="ground-truth TSV from make_benchmark_reads.py")
    parser.add_argument("--per-taxid", action="store_true", help="also print a per-taxid breakdown")
    return parser.parse_args()


def load_ground_truth(path):
    truth = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            truth[row["read_id"]] = row["true_taxid"]
    return truth


def load_predictions(path):
    predictions = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            read_id = row[0].strip()
            assignment = row[2].strip()
            predictions[read_id] = assignment
    return predictions


def main():
    args = parse_args()
    truth = load_ground_truth(args.ground_truth)
    if not truth:
        raise SystemExit("no ground-truth entries found in %s" % args.ground_truth)

    predictions = load_predictions(args.results)

    tp = fp = fn = 0
    per_taxid = {}

    def bucket(taxid):
        return per_taxid.setdefault(taxid, {"tp": 0, "fp": 0, "fn": 0})

    for read_id, true_taxid in truth.items():
        predicted = predictions.get(read_id)
        if predicted is None:
            print("Warning: no prediction found for read %s (treated as unclassified)" % read_id, file=sys.stderr)
            predicted = "NA"

        if predicted in UNCLASSIFIED:
            fn += 1
            bucket(true_taxid)["fn"] += 1
            continue

        if predicted == true_taxid:
            tp += 1
            bucket(true_taxid)["tp"] += 1
        else:
            fn += 1
            fp += 1
            bucket(true_taxid)["fn"] += 1
            bucket(predicted)["fp"] += 1

    sensitivity = tp / (tp + fn) if (tp + fn) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0

    print("reads_total=%d tp=%d fp=%d fn=%d sensitivity=%.6f precision=%.6f" % (
        len(truth), tp, fp, fn, sensitivity, precision,
    ))

    if args.per_taxid:
        for taxid in sorted(per_taxid.keys()):
            counts = per_taxid[taxid]
            t_tp, t_fp, t_fn = counts["tp"], counts["fp"], counts["fn"]
            t_sens = t_tp / (t_tp + t_fn) if (t_tp + t_fn) else 0.0
            t_prec = t_tp / (t_tp + t_fp) if (t_tp + t_fp) else 0.0
            print("taxid=%s tp=%d fp=%d fn=%d sensitivity=%.6f precision=%.6f" % (
                taxid, t_tp, t_fp, t_fn, t_sens, t_prec,
            ))


if __name__ == "__main__":
    main()
