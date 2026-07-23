#!/usr/bin/env python3
"""Repair absolute-path metadata after moving/copying a CLARK database directory.

scripts/set_targets.sh caches several hidden metadata files under the
database directory (e.g. ".protozoa", ".protozoa.fileToAccssnTaxID",
".protozoa.fileToTaxIDs") whose *first field on every line is an absolute
path* to a genome FASTA file, recorded at download/scan time. If you copy
or rsync the whole database directory to a different machine (or a
different path on the same machine), those cached files still contain the
old absolute paths -- set_targets.sh sees the cache is non-empty, skips
re-downloading/re-scanning, and CLARK ends up looking for genomes at a
path that no longer exists (e.g. "/home/h/.../Protozoa/foo.fna" when the
directory is now under "/home/u/...").

This script repairs those cached files in place: for each first-field path
that no longer exists on disk, it looks for the same relative suffix (type
subdirectory + filename, or just the filename) under the *current* database
directory and rewrites the line with the corrected path. It never touches
lines whose path already resolves correctly, and it leaves -- with a
warning -- any line whose file genuinely cannot be found anywhere under the
new database directory (a real missing/deleted genome, not a path issue).

After running this, re-run scripts/set_targets.sh as usual (no network
access needed) to regenerate a clean targets.txt from the repaired cache.
"""

import argparse
import os
import sys

DATABASES = ("bacteria", "viruses", "plasmid", "plastid", "protozoa", "fungi", "human", "custom")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("db_dir", help="database directory previously passed to scripts/set_targets.sh")
    parser.add_argument("--dry-run", action="store_true", help="report what would change without writing anything")
    return parser.parse_args()


def candidate_files(db_dir):
    files = []
    for db in DATABASES:
        for suffix in ("", ".fileToAccssnTaxID", ".fileToTaxIDs"):
            path = os.path.join(db_dir, "." + db + suffix)
            if os.path.isfile(path):
                files.append(path)
    return files


def resolve_path(old_path, db_dir):
    """Find the same file under db_dir by trying progressively shorter
    path suffixes of old_path, longest (most specific) first."""
    parts = [p for p in old_path.split("/") if p]
    for k in range(len(parts) - 1, 0, -1):
        candidate = os.path.join(db_dir, *parts[-k:])
        if os.path.isfile(candidate):
            return candidate
    return None


def rehome_file(path, db_dir, dry_run):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    fixed = 0
    unresolved = 0
    out_lines = []
    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped:
            out_lines.append(line)
            continue
        fields = stripped.split("\t")
        old_path = fields[0]
        if os.path.isfile(old_path):
            out_lines.append(line)
            continue
        new_path = resolve_path(old_path, db_dir)
        if new_path is None:
            print("Warning: %s: could not resolve '%s' under %s (leaving unchanged)" % (path, old_path, db_dir), file=sys.stderr)
            unresolved += 1
            out_lines.append(line)
            continue
        fields[0] = new_path
        out_lines.append("\t".join(fields) + "\n")
        fixed += 1

    if fixed and not dry_run:
        tmp_path = path + ".rehome.tmp"
        with open(tmp_path, "w", encoding="utf-8") as handle:
            handle.writelines(out_lines)
        os.replace(tmp_path, path)

    return fixed, unresolved


def main():
    args = parse_args()
    db_dir = os.path.abspath(args.db_dir)
    if not os.path.isdir(db_dir):
        raise SystemExit("database directory '%s' does not exist" % db_dir)

    files = candidate_files(db_dir)
    if not files:
        raise SystemExit("no cached metadata files (.bacteria, .viruses, ...) found under %s" % db_dir)

    total_fixed = 0
    total_unresolved = 0
    for path in files:
        fixed, unresolved = rehome_file(path, db_dir, args.dry_run)
        total_fixed += fixed
        total_unresolved += unresolved
        if fixed:
            verb = "would fix" if args.dry_run else "fixed"
            print("%s: %s %d path(s)" % (path, verb, fixed), file=sys.stderr)

    if total_fixed == 0 and total_unresolved == 0:
        print("No stale paths found under %s; nothing to do." % db_dir, file=sys.stderr)
    else:
        print(
            "Total: %d path(s) %s, %d unresolved."
            % (total_fixed, "would be fixed" if args.dry_run else "fixed", total_unresolved),
            file=sys.stderr,
        )
    if total_unresolved:
        print("Some genome files could not be found under the new database directory; re-download those, or investigate before rerunning set_targets.sh.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
