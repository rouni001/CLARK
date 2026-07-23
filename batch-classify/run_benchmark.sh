#!/usr/bin/env bash
# Like run_all.sh, but samples reads from scripts/make_benchmark_reads.sh
# (a benchmark profile with published, citable specs -- see that script's
# docstring) instead of scripts/make_sample.sh's ad-hoc uniform sampling,
# and reports classification sensitivity/precision via scripts/eval_accuracy.py
# alongside the usual CLARK results. See README.md in this folder for details.

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: batch-classify/run_benchmark.sh <db-directory> <profile> [type ...]

<db-directory> is the directory previously passed to scripts/set_targets.sh
(the one containing each type's genomes, taxonomy/, and <type>_0/ database
folders). <profile> is one of: hiseq, miseq, simba5, simhc, custom (see
scripts/make_benchmark_reads.py). Types default to: viruses plasmid plastid
fungi human.

For each type, this:
  1. Re-runs scripts/set_targets.sh <db-directory> <type> --species to
     regenerate targets.txt/.settings for that type (no network access or
     rebuilding needed -- sequences and taxonomy are already downloaded).
  2. Samples a benchmark objects.fa + ground-truth TSV from that type's own
     targets via scripts/make_benchmark_reads.sh.
  3. Classifies it with exe/CLARK -k <kmer> -T <targets> -D <DBD/> -O
     <objects.fa> -R <results> -n <threads>.
  4. Scores the result against the ground truth with scripts/eval_accuracy.py,
     writing accuracy.txt and echoing it to stdout.

Results land in batch-classify/results/<type>/.

Environment overrides:
  CLARK_HOME              repository root (default: parent of this script's dir)
  CLARK_VARIANT_EXE       CLARK executable name (default: CLARK)
  CLARK_KMER              k-mer size (default: 31)
  CLARK_THREADS           thread count (default: 8)
  CLARK_RANK              taxonomy rank flag for set_targets.sh (default: --species)
  CLARK_BENCHMARK_COUNT   override the profile's default read count
  CLARK_BENCHMARK_LEN     override the profile's default read length
  CLARK_BENCHMARK_ERROR_RATE  override the profile's default error rate
  CLARK_BENCHMARK_NGENOMES    override the profile's default source-genome count
  CLARK_BENCHMARK_SEED    random seed passed to make_benchmark_reads.sh
  CLARK_BATCH_RESULTS_DIR where per-type results land (default: this script's results/)
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

if [ "$#" -lt 2 ]; then
	usage
	exit 1
fi

DBDR_INPUT="$1"
PROFILE="$2"
shift 2

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

[ -d "$DBDR_INPUT" ] || die "database directory '$DBDR_INPUT' does not exist"
DBDR="$(cd -P "$DBDR_INPUT" >/dev/null 2>&1 && pwd)"

CLARK_EXE="$LDIR/exe/$VARIANT_EXE"
[ -x "$CLARK_EXE" ] || die "missing executable '$CLARK_EXE'. Run 'make all' from the CLARK repository root first."

SETTINGS_FILE="$LDIR/.settings"
OUT_DIR="${CLARK_BATCH_RESULTS_DIR:-$SCRIPT_DIR/results}"
mkdir -p "$OUT_DIR"

BENCHMARK_ARGS=()
[ -n "${CLARK_BENCHMARK_COUNT:-}" ] && BENCHMARK_ARGS+=(-n "$CLARK_BENCHMARK_COUNT")
[ -n "${CLARK_BENCHMARK_LEN:-}" ] && BENCHMARK_ARGS+=(-l "$CLARK_BENCHMARK_LEN")
[ -n "${CLARK_BENCHMARK_ERROR_RATE:-}" ] && BENCHMARK_ARGS+=(-e "$CLARK_BENCHMARK_ERROR_RATE")
[ -n "${CLARK_BENCHMARK_NGENOMES:-}" ] && BENCHMARK_ARGS+=(-g-n "$CLARK_BENCHMARK_NGENOMES")
[ -n "${CLARK_BENCHMARK_SEED:-}" ] && BENCHMARK_ARGS+=(--seed "$CLARK_BENCHMARK_SEED")

for type in "${TYPES[@]}"; do
	echo "== $type ($PROFILE) =="

	CLARK_HOME="$LDIR" "$LDIR/scripts/set_targets.sh" "$DBDR" "$type" "$RANK_FLAG"

	targets="$DBDR/targets.txt"
	dbd="$(awk '$1 == "-D" { print $2; exit }' "$SETTINGS_FILE")"
	[ -n "$dbd" ] || die "failed to determine the database directory for '$type' from $SETTINGS_FILE"

	type_out="$OUT_DIR/$type"
	mkdir -p "$type_out"
	objects="$type_out/objects.fa"
	truth="$type_out/objects.fa.truth.tsv"
	results="$type_out/results"
	accuracy="$type_out/accuracy.txt"

	"$LDIR/scripts/make_benchmark_reads.sh" -p "$PROFILE" -T "$targets" -o "$objects" -g "$truth" "${BENCHMARK_ARGS[@]}"

	"$CLARK_EXE" -k "$KMER" -T "$targets" -D "$dbd" -O "$objects" -R "$results" -n "$THREADS"

	python3 "$LDIR/scripts/eval_accuracy.py" --results "$results.csv" --ground-truth "$truth" --per-taxid | tee "$accuracy"

	echo "-- $type done: $results.csv, $accuracy"
done
