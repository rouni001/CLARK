#!/usr/bin/env bash
# Classifies a freshly sampled objects.fa against each of several already-built
# CLARK databases living under one shared database directory (as produced by
# running scripts/set_targets.sh once per database type against that same
# directory). See README.md in this folder for details.

set -euo pipefail

# EDIT THIS if you use -g/--gpu: absolute path to the cuCLARK GPU executable.
# Everything else (targets, db dir, objects, results, -k/-n) is passed to it
# exactly as it would be to CLARK, since cuCLARK accepts the same CLI.
CUCLARK_EXE="${CUCLARK_EXE:-/path/to/cuCLARK}"

usage() {
	cat <<'USAGE'
Usage: batch-classify/run_all.sh -d <db-dir> [options]

Required:
  -d <path>   Database directory previously passed to scripts/set_targets.sh
              (the one containing each type's genomes, taxonomy/, and
              <type>_0/ database folders).

Options:
  -t <list>   Comma- or space-separated database types to run.
              (default: viruses,plasmid,plastid,fungi,human)
  -n <n>      Thread count passed to -n. (default: 8)
  -c <n>      Number of reads to sample per type. (default: 20)
  -l <n>      Sampled read length. (default: 150)
  -k <n>      k-mer size passed to -k. (default: 31)
  -p          Enable the opt-in CuD performance breakdown (build/load/
              match/write time, plus match/hits-update/classify timing)
              by setting CLARK_CUD_PROFILE=1 and CLARK_CUD_PROFILE_FINE=1
              for the run. CPU (CLARK) only; ignored with -g. (default: off)
  -g          Run on GPU with cuCLARK instead of CPU CLARK. Edit
              CUCLARK_EXE at the top of this script to point at your
              cuCLARK binary. Overrides -x. (default: off)
  -x <name>   CPU executable variant: CLARK, CLARK-l, CLARK-S.
              Ignored if -g is set. (default: CLARK)
  -r <flag>   Taxonomy rank flag for set_targets.sh. (default: --species)
  -o <path>   Where per-type results land.
              (default: this script's results/ directory)
  -h          Show this message.

For each type, this:
  1. Re-runs scripts/set_targets.sh <db-dir> <type> <rank> to regenerate
     targets.txt/.settings for that type (no network access or rebuilding
     needed -- sequences and taxonomy are already downloaded).
  2. Samples a fresh objects.fa from that type's own targets via
     scripts/make_sample.sh.
  3. Classifies it with the chosen executable: -k <kmer> -T <targets>
     -D <DBD/> -O <objects.fa> -R <results> -n <threads>.

Results land in <results-dir>/<type>/.

CLARK_HOME (env var) overrides the repository root; defaults to the
parent of this script's directory.
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

require_positive_int() {
	case "$2" in
		''|*[!0-9]*) die "$1 requires a positive integer, got '$2'" ;;
	esac
	[ "$2" -gt 0 ] || die "$1 requires a positive integer, got '$2'"
}

DBDR_INPUT=""
TYPES_RAW=""
THREADS=8
SAMPLE_COUNT=20
SAMPLE_LEN=150
KMER=31
PROFILE=0
GPU=0
VARIANT_EXE="CLARK"
RANK_FLAG="--species"
RESULTS_DIR=""

if [ "$#" -eq 0 ]; then
	usage
	exit 1
fi

while [ "$#" -gt 0 ]; do
	case "$1" in
		-d)
			[ "$#" -ge 2 ] || die "-d requires a path"
			DBDR_INPUT="$2"
			shift 2
			;;
		-t)
			[ "$#" -ge 2 ] || die "-t requires a comma- or space-separated list"
			TYPES_RAW="$2"
			shift 2
			;;
		-n)
			[ "$#" -ge 2 ] || die "-n requires a positive integer"
			require_positive_int "-n" "$2"
			THREADS="$2"
			shift 2
			;;
		-c)
			[ "$#" -ge 2 ] || die "-c requires a positive integer"
			require_positive_int "-c" "$2"
			SAMPLE_COUNT="$2"
			shift 2
			;;
		-l)
			[ "$#" -ge 2 ] || die "-l requires a positive integer"
			require_positive_int "-l" "$2"
			SAMPLE_LEN="$2"
			shift 2
			;;
		-k)
			[ "$#" -ge 2 ] || die "-k requires a positive integer"
			require_positive_int "-k" "$2"
			KMER="$2"
			shift 2
			;;
		-p)
			PROFILE=1
			shift
			;;
		-g)
			GPU=1
			shift
			;;
		-x)
			[ "$#" -ge 2 ] || die "-x requires an executable name"
			VARIANT_EXE="$2"
			shift 2
			;;
		-r)
			[ "$#" -ge 2 ] || die "-r requires a taxonomy rank flag"
			RANK_FLAG="$2"
			shift 2
			;;
		-o)
			[ "$#" -ge 2 ] || die "-o requires a path"
			RESULTS_DIR="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			die "unrecognized option: $1 (see -h)"
			;;
	esac
done

[ -n "$DBDR_INPUT" ] || die "-d <db-dir> is required (see -h)"

TYPES=()
if [ -n "$TYPES_RAW" ]; then
	IFS=', ' read -r -a TYPES <<< "$TYPES_RAW"
fi
if [ "${#TYPES[@]}" -eq 0 ]; then
	TYPES=(viruses plasmid plastid fungi human)
fi

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"

[ -d "$DBDR_INPUT" ] || die "database directory '$DBDR_INPUT' does not exist"
DBDR="$(cd -P "$DBDR_INPUT" >/dev/null 2>&1 && pwd)"

if [ "$GPU" -eq 1 ]; then
	CLASSIFY_EXE="$CUCLARK_EXE"
	[ -x "$CLASSIFY_EXE" ] || die "missing GPU executable '$CLASSIFY_EXE'. Edit CUCLARK_EXE at the top of this script to point at your cuCLARK binary."
else
	CLASSIFY_EXE="$LDIR/exe/$VARIANT_EXE"
	[ -x "$CLASSIFY_EXE" ] || die "missing executable '$CLASSIFY_EXE'. Run 'make all' from the CLARK repository root first."
fi

SETTINGS_FILE="$LDIR/.settings"
OUT_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
mkdir -p "$OUT_DIR"

PROFILE_ENV=()
if [ "$PROFILE" -eq 1 ] && [ "$GPU" -eq 0 ]; then
	PROFILE_ENV=(CLARK_CUD_PROFILE=1 CLARK_CUD_PROFILE_FINE=1)
fi

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

	env "${PROFILE_ENV[@]}" "$CLASSIFY_EXE" -k "$KMER" -T "$targets" -D "$dbd" -O "$objects" -R "$results" -n "$THREADS"

	echo "-- $type done: $results.csv"
done
