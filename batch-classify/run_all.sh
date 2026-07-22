#!/usr/bin/env bash
# Classifies a freshly sampled objects.fa against each of several already-built
# CLARK databases living under one shared database directory (as produced by
# running scripts/set_targets.sh once per database type against that same
# directory). See README.md in this folder for details.

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: batch-classify/run_all.sh <db-directory> [type ...]

<db-directory> is the directory previously passed to scripts/set_targets.sh
(the one containing each type's genomes, taxonomy/, and <type>_0/ database
folders). Types default to: viruses plasmid plastid fungi human.

For each type, this:
  1. Re-runs scripts/set_targets.sh <db-directory> <type> --species to
     regenerate targets.txt/.settings for that type (no network access or
     rebuilding needed -- sequences and taxonomy are already downloaded).
  2. Samples a fresh objects.fa from that type's own targets via
     scripts/make_sample.sh.
  3. Classifies it with exe/CLARK -k <kmer> -T <targets> -D <DBD/> -O
     <objects.fa> -R <results> -n <threads>.

Results land in batch-classify/results/<type>/.

Environment overrides:
  CLARK_HOME          repository root (default: parent of this script's dir)
  CLARK_VARIANT_EXE   CLARK executable name (default: CLARK)
  CLARK_KMER          k-mer size (default: 31)
  CLARK_THREADS       thread count (default: 8)
  CLARK_RANK          taxonomy rank flag for set_targets.sh (default: --species)
  CLARK_SAMPLE_COUNT  number of sampled reads per type (default: 20)
  CLARK_SAMPLE_LEN    sampled read length (default: 150)
  CLARK_BATCH_RESULTS_DIR  where per-type results land (default: this script's results/)
USAGE
}

script_dir() {
	local source="${BASH_SOURCE[0]}"
	while [ -h "$source" ]; do
		local dir
		dir="$(cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd)"
		source="$(readlink "$source")"
		[[ "$source" != /* ]] && source="$dir/$source"
	done
	cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd
}

die() {
	echo "Error: $*" >&2
	exit 1
}

if [ "$#" -lt 1 ]; then
	usage
	exit 1
fi

DBDR_INPUT="$1"
shift

TYPES=("$@")
if [ "${#TYPES[@]}" -eq 0 ]; then
	TYPES=(viruses plasmid plastid fungi human)
fi

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"
VARIANT_EXE="${CLARK_VARIANT_EXE:-CLARK}"
KMER="${CLARK_KMER:-31}"
THREADS="${CLARK_THREADS:-8}"
RANK_FLAG="${CLARK_RANK:---species}"
SAMPLE_COUNT="${CLARK_SAMPLE_COUNT:-20}"
SAMPLE_LEN="${CLARK_SAMPLE_LEN:-150}"

[ -d "$DBDR_INPUT" ] || die "database directory '$DBDR_INPUT' does not exist"
DBDR="$(cd -P "$DBDR_INPUT" >/dev/null 2>&1 && pwd)"

CLARK_EXE="$LDIR/exe/$VARIANT_EXE"
[ -x "$CLARK_EXE" ] || die "missing executable '$CLARK_EXE'. Run 'make all' from the CLARK repository root first."

SETTINGS_FILE="$LDIR/.settings"
OUT_DIR="${CLARK_BATCH_RESULTS_DIR:-$SCRIPT_DIR/results}"
mkdir -p "$OUT_DIR"

for type in "${TYPES[@]}"; do
	echo "== $type =="

	CLARK_HOME="$LDIR" "$LDIR/scripts/set_targets.sh" "$DBDR" "$type" "$RANK_FLAG"

	targets="$DBDR/targets.txt"
	dbd="$(awk '$1 == "-D" { print $2; exit }' "$SETTINGS_FILE")"
	[ -n "$dbd" ] || die "failed to determine the database directory for '$type' from $SETTINGS_FILE"

	type_out="$OUT_DIR/$type"
	mkdir -p "$type_out"
	objects="$type_out/objects.fa"
	results="$type_out/results"

	"$LDIR/scripts/make_sample.sh" -n "$SAMPLE_COUNT" -l "$SAMPLE_LEN" -o "$objects" -T "$targets"

	"$CLARK_EXE" -k "$KMER" -T "$targets" -D "$dbd" -O "$objects" -R "$results" -n "$THREADS"

	echo "-- $type done: $results.csv"
done
