#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/make_sample.sh -n <count> [options]

Builds a small FASTA file by sampling reads from the genomes CLARK is
currently configured to classify against (as set by
scripts/set_targets.sh), for smoke-testing scripts/classify_metagenome.sh.

Options:
  -n <count>          Number of reads to sample (required).
  -l <length>         Read length in bases (default: 150).
  -o <file>           Output FASTA path (default: ./sample.fa).
  -T <targets.txt>    Override the targets file (default: read from .settings).
  --seed <n>          Random seed for reproducible sampling.
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

COUNT=""
LENGTH=150
OUTPUT="sample.fa"
TARGETS=""
SEED=""

while [ "$#" -gt 0 ]; do
	case "$1" in
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
		-o)
			[ "$#" -ge 2 ] || die "-o requires a file path"
			OUTPUT="$2"
			shift 2
			;;
		-T)
			[ "$#" -ge 2 ] || die "-T requires a targets file"
			TARGETS="$2"
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

case "$COUNT" in
	''|*[!0-9]*) die "-n requires a positive integer" ;;
esac
[ "$COUNT" -gt 0 ] || die "-n requires a positive integer"

case "$LENGTH" in
	''|*[!0-9]*) die "-l requires a positive integer" ;;
esac
[ "$LENGTH" -gt 0 ] || die "-l requires a positive integer"

if [ -z "$TARGETS" ]; then
	[ -s "$SETTINGS_FILE" ] || die "targets are not configured; run scripts/set_targets.sh first, or pass -T <targets.txt>"
	TARGETS="$(awk '$1 == "-T" { $1 = ""; sub(/^ /, ""); print; exit }' "$SETTINGS_FILE")"
	[ -n "$TARGETS" ] || die "failed to find a targets file entry in '$SETTINGS_FILE'"
fi

[ -s "$TARGETS" ] || die "targets file '$TARGETS' is missing or empty"
command -v "$PYTHON_CMD" >/dev/null 2>&1 || die "Python 3 is required to build the sample"

SEED_ARGS=()
if [ -n "$SEED" ]; then
	SEED_ARGS=(--seed "$SEED")
fi

"$PYTHON_CMD" "$SCRIPT_DIR/make_sample.py" \
	--targets "$TARGETS" \
	--count "$COUNT" \
	--length "$LENGTH" \
	--output "$OUTPUT" \
	"${SEED_ARGS[@]}"
