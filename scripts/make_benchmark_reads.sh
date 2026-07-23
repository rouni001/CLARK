#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/make_benchmark_reads.sh -p <profile> [options]

Builds a FASTA of simulated reads matching the published characteristics of
an established metagenome-classifier benchmark (see scripts/make_benchmark_reads.py
for the full rationale and citations), sampled from the genomes CLARK is
currently configured to classify against (as set by scripts/set_targets.sh).
Also writes a ground-truth TSV (read_id -> true taxid) for
scripts/eval_accuracy.py.

Profiles: hiseq, miseq, simba5, simhc, custom

Options:
  -p <profile>        Benchmark profile (required).
  -o <file>           Output FASTA path (default: ./benchmark.fa).
  -g <file>           Ground-truth TSV path (default: <output>.truth.tsv).
  -T <targets.txt>    Override the targets file (default: read from .settings).
  -n <count>          Override the profile's default read count.
  -l <length>         Override the profile's default read length.
  -e <error-rate>     Override the profile's default per-base substitution rate.
  -g-n <n-genomes>    Override the profile's default number of source genomes.
  --seed <n>          Random seed for reproducible sampling.

Required for -p custom: -n, -l, -e.
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

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"
SETTINGS_FILE="${CLARK_SETTINGS_FILE:-$LDIR/.settings}"
PYTHON_CMD="${CLARK_PYTHON:-python3}"

PROFILE=""
OUTPUT="benchmark.fa"
GROUND_TRUTH=""
TARGETS=""
COUNT=""
LENGTH=""
ERROR_RATE=""
NGENOMES=""
SEED=""

while [ "$#" -gt 0 ]; do
	case "$1" in
		-p)
			[ "$#" -ge 2 ] || die "-p requires a profile name"
			PROFILE="$2"
			shift 2
			;;
		-o)
			[ "$#" -ge 2 ] || die "-o requires a file path"
			OUTPUT="$2"
			shift 2
			;;
		-g)
			[ "$#" -ge 2 ] || die "-g requires a file path"
			GROUND_TRUTH="$2"
			shift 2
			;;
		-T)
			[ "$#" -ge 2 ] || die "-T requires a targets file"
			TARGETS="$2"
			shift 2
			;;
		-n)
			[ "$#" -ge 2 ] || die "-n requires a positive integer"
			COUNT="$2"
			shift 2
			;;
		-l)
			[ "$#" -ge 2 ] || die "-l requires a positive integer"
			LENGTH="$2"
			shift 2
			;;
		-e)
			[ "$#" -ge 2 ] || die "-e requires a number"
			ERROR_RATE="$2"
			shift 2
			;;
		-g-n)
			[ "$#" -ge 2 ] || die "-g-n requires a positive integer"
			NGENOMES="$2"
			shift 2
			;;
		--seed)
			[ "$#" -ge 2 ] || die "--seed requires an integer"
			SEED="$2"
			shift 2
			;;
		--help|-h)
			usage
			exit 0
			;;
		*)
			die "unrecognized option: $1"
			;;
	esac
done

[ -n "$PROFILE" ] || die "-p <profile> is required"

if [ -z "$GROUND_TRUTH" ]; then
	GROUND_TRUTH="$OUTPUT.truth.tsv"
fi

if [ -z "$TARGETS" ]; then
	[ -s "$SETTINGS_FILE" ] || die "targets are not configured; run scripts/set_targets.sh first, or pass -T <targets.txt>"
	TARGETS="$(awk '$1 == "-T" { $1 = ""; sub(/^ /, ""); print; exit }' "$SETTINGS_FILE")"
	[ -n "$TARGETS" ] || die "failed to find a targets file entry in '$SETTINGS_FILE'"
fi

[ -s "$TARGETS" ] || die "targets file '$TARGETS' is missing or empty"
command -v "$PYTHON_CMD" >/dev/null 2>&1 || die "Python 3 is required to build the benchmark reads"

EXTRA_ARGS=()
[ -n "$COUNT" ] && EXTRA_ARGS+=(--count "$COUNT")
[ -n "$LENGTH" ] && EXTRA_ARGS+=(--length "$LENGTH")
[ -n "$ERROR_RATE" ] && EXTRA_ARGS+=(--error-rate "$ERROR_RATE")
[ -n "$NGENOMES" ] && EXTRA_ARGS+=(--n-genomes "$NGENOMES")
[ -n "$SEED" ] && EXTRA_ARGS+=(--seed "$SEED")

"$PYTHON_CMD" "$SCRIPT_DIR/make_benchmark_reads.py" \
	--targets "$TARGETS" \
	--profile "$PROFILE" \
	--output "$OUTPUT" \
	--ground-truth "$GROUND_TRUTH" \
	"${EXTRA_ARGS[@]}"
