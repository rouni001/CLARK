#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/classify_metagenome.sh -O <fileObjects> -R <fileResults> [options]
       scripts/classify_metagenome.sh -P <file1> <file2> -R <fileResults> [options]

Common options:
  -k <kmerSize>        k-mer length for CLARK
  -t <minFreqTarget>   minimum k-mer frequency in targets
  -o <minFreqObject>   minimum k-mer frequency in objects
  -m <mode>            0 full, 1 default, 2 express, or 3 spectrum
  -n <threads>         number of threads
  -g <iteration>       gap for CLARK-l
  -s <factor>          sampling factor
  --long               enable long-read memory mode
  --ldm                load database by memory-mapped file
  --kso                preliminary k-spectrum analysis for mode 3
  --extended           extended full-mode output
  --light              run CLARK-l
  --spaced             run CLARK-S
  --gzipped            decompress input reads before classification

Run scripts/set_targets.sh <DB_DIR> <database choice...> before classification.
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

decompress_input() {
	local src="$1"
	local label="$2"
	local dst="$TMPDIR_CLARK/$label"

	[ -f "$src" ] || die "failed to open input file '$src'"
	gzip -dc "$src" > "$dst" || die "failed to decompress '$src'"
	[ -s "$dst" ] || die "decompressed file is empty for '$src'"
	printf '%s\n' "$dst"
}

trim_leading_space() {
	local value="$1"
	value="${value#"${value%%[![:space:]]*}"}"
	printf '%s\n' "$value"
}

if [ "$#" -lt 2 ]; then
	usage
	exit 1
fi

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"
SETTINGS_FILE="${CLARK_SETTINGS_FILE:-$LDIR/.settings}"
EXE_DIR="${CLARK_EXE_DIR:-$LDIR/exe}"

[ -s "$SETTINGS_FILE" ] || die "targets are not configured. Run scripts/set_targets.sh first."

PARAMS=()
while IFS= read -r line || [ -n "$line" ]; do
	[ -n "$line" ] || continue
	key="${line%%[[:space:]]*}"
	value="${line#"$key"}"
	value="$(trim_leading_space "$value")"
	[ -n "$key" ] || continue
	[ -n "$value" ] || die "malformed settings line in '$SETTINGS_FILE': $line"
	PARAMS+=("$key" "$value")
done < "$SETTINGS_FILE"

VARIANT="DEFAULT"
GZIPPED=0

for arg in "$@"; do
	case "$arg" in
		--gzipped)
			GZIPPED=1
			;;
		--light)
			[ "$VARIANT" != "SPACED" ] || die "--light and --spaced cannot be used together"
			VARIANT="LIGHT"
			;;
		--spaced)
			[ "$VARIANT" != "LIGHT" ] || die "--light and --spaced cannot be used together"
			VARIANT="SPACED"
			;;
	esac
done

TMPDIR_CLARK=""
cleanup() {
	if [ -n "$TMPDIR_CLARK" ] && [ -d "$TMPDIR_CLARK" ]; then
		rm -rf "$TMPDIR_CLARK"
	fi
}
trap cleanup EXIT INT TERM

if [ "$GZIPPED" -eq 1 ]; then
	TMPDIR_CLARK="$(mktemp -d "${TMPDIR:-/tmp}/CLARKGZP.XXXXXX")"
fi

ARGS=("$@")
i=0
while [ "$i" -lt "${#ARGS[@]}" ]; do
	var="${ARGS[$i]}"
	case "$var" in
		-T)
			die "classify_metagenome.sh does not accept -T. Targets are set by set_targets.sh."
			;;
		-D)
			die "classify_metagenome.sh does not accept -D. The database directory is set by set_targets.sh."
			;;
		-O)
			PARAMS+=("-O")
			i=$((i + 1))
			[ "$i" -lt "${#ARGS[@]}" ] || die "-O requires a file"
			input="${ARGS[$i]}"
			[ -f "$input" ] || die "failed to open input file '$input'"
			if [ "$GZIPPED" -eq 1 ]; then
				PARAMS+=("$(decompress_input "$input" "objects.fa")")
			else
				PARAMS+=("$input")
			fi
			;;
		-P)
			PARAMS+=("-P")
			i=$((i + 1))
			[ "$i" -lt "${#ARGS[@]}" ] || die "-P requires two files"
			input1="${ARGS[$i]}"
			i=$((i + 1))
			[ "$i" -lt "${#ARGS[@]}" ] || die "-P requires two files"
			input2="${ARGS[$i]}"
			[ -f "$input1" ] || die "failed to open paired input file '$input1'"
			[ -f "$input2" ] || die "failed to open paired input file '$input2'"
			if [ "$GZIPPED" -eq 1 ]; then
				PARAMS+=("$(decompress_input "$input1" "mate1.fq")")
				PARAMS+=("$(decompress_input "$input2" "mate2.fq")")
			else
				PARAMS+=("$input1" "$input2")
			fi
			;;
		--gzipped|--light|--spaced)
			;;
		*)
			PARAMS+=("$var")
			;;
	esac
	i=$((i + 1))
done

case "$VARIANT" in
	LIGHT)
		CLARK_BINARY="CLARK-l"
		;;
	SPACED)
		CLARK_BINARY="CLARK-S"
		;;
	*)
		CLARK_BINARY="CLARK"
		;;
esac

CLARK_EXE="$EXE_DIR/$CLARK_BINARY"
[ -x "$CLARK_EXE" ] || die "missing executable '$CLARK_EXE'. Run scripts/install.sh first."

"$CLARK_EXE" "${PARAMS[@]}"
