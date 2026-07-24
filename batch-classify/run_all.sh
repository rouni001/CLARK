#!/usr/bin/env bash
# Classifies a freshly sampled objects.fa against each of several already-built
# CLARK databases living under one shared database directory (as produced by
# running scripts/set_targets.sh once per database type against that same
# directory). See README.md in this folder for details.

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: batch-classify/run_all.sh --db-dir <path> [options]

Required:
  --db-dir <path>        Database directory previously passed to
                          scripts/set_targets.sh (the one containing each
                          type's genomes, taxonomy/, and <type>_0/ database
                          folders).

Options:
  --types <t1,t2,...>     Comma- or space-separated database types to run.
                          (default: viruses,plasmid,plastid,fungi,human)
  --threads <n>           Thread count passed to CLARK's -n. (default: 8)
  --sample-count <n>      Number of reads to sample per type. (default: 20)
  --sample-len <n>        Sampled read length. (default: 150)
  --kmer <n>              k-mer size passed to CLARK's -k. (default: 31)
  --profile               Enable CLARK's opt-in CuD performance breakdown
                          (build/load/match/write time, plus the fine
                          match/hits-update/classify timing) by setting
                          CLARK_CUD_PROFILE=1 and CLARK_CUD_PROFILE_FINE=1
                          for the CLARK invocation. (default: off)
  --variant-exe <name>    CLARK executable to run: CLARK, CLARK-l, CLARK-S.
                          (default: CLARK)
  --rank <flag>           Taxonomy rank flag for set_targets.sh.
                          (default: --species)
  --results-dir <path>    Where per-type results land.
                          (default: this script's results/ directory)
  -h, --help              Show this message.

For each type, this:
  1. Re-runs scripts/set_targets.sh <db-dir> <type> <rank> to regenerate
     targets.txt/.settings for that type (no network access or rebuilding
     needed -- sequences and taxonomy are already downloaded).
  2. Samples a fresh objects.fa from that type's own targets via
     scripts/make_sample.sh.
  3. Classifies it with exe/<variant> -k <kmer> -T <targets> -D <DBD/> -O
     <objects.fa> -R <results> -n <threads>.

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
VARIANT_EXE="CLARK"
RANK_FLAG="--species"
RESULTS_DIR=""

if [ "$#" -eq 0 ]; then
	usage
	exit 1
fi

while [ "$#" -gt 0 ]; do
	case "$1" in
		--db-dir)
			[ "$#" -ge 2 ] || die "--db-dir requires a path"
			DBDR_INPUT="$2"
			shift 2
			;;
		--types)
			[ "$#" -ge 2 ] || die "--types requires a comma- or space-separated list"
			TYPES_RAW="$2"
			shift 2
			;;
		--threads)
			[ "$#" -ge 2 ] || die "--threads requires a positive integer"
			require_positive_int "--threads" "$2"
			THREADS="$2"
			shift 2
			;;
		--sample-count)
			[ "$#" -ge 2 ] || die "--sample-count requires a positive integer"
			require_positive_int "--sample-count" "$2"
			SAMPLE_COUNT="$2"
			shift 2
			;;
		--sample-len)
			[ "$#" -ge 2 ] || die "--sample-len requires a positive integer"
			require_positive_int "--sample-len" "$2"
			SAMPLE_LEN="$2"
			shift 2
			;;
		--kmer)
			[ "$#" -ge 2 ] || die "--kmer requires a positive integer"
			require_positive_int "--kmer" "$2"
			KMER="$2"
			shift 2
			;;
		--profile)
			PROFILE=1
			shift
			;;
		--variant-exe)
			[ "$#" -ge 2 ] || die "--variant-exe requires an executable name"
			VARIANT_EXE="$2"
			shift 2
			;;
		--rank)
			[ "$#" -ge 2 ] || die "--rank requires a taxonomy rank flag"
			RANK_FLAG="$2"
			shift 2
			;;
		--results-dir)
			[ "$#" -ge 2 ] || die "--results-dir requires a path"
			RESULTS_DIR="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			die "unrecognized option: $1 (see --help)"
			;;
	esac
done

[ -n "$DBDR_INPUT" ] || die "--db-dir is required (see --help)"

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

CLARK_EXE="$LDIR/exe/$VARIANT_EXE"
[ -x "$CLARK_EXE" ] || die "missing executable '$CLARK_EXE'. Run 'make all' from the CLARK repository root first."

SETTINGS_FILE="$LDIR/.settings"
OUT_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
mkdir -p "$OUT_DIR"

PROFILE_ENV=()
if [ "$PROFILE" -eq 1 ]; then
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

	env "${PROFILE_ENV[@]}" "$CLARK_EXE" -k "$KMER" -T "$targets" -D "$dbd" -O "$objects" -R "$results" -n "$THREADS"

	echo "-- $type done: $results.csv"
done
