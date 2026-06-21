#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/set_targets.sh <Directory_path> <database choice...> [taxonomy rank]

Database choices:
  bacteria viruses plasmid plastid protozoa fungi human custom

Taxonomy rank:
  --species (default), --genus, --family, --order, --class, --phylum

Download options for RefSeq databases:
  --download-threads <N>              Download up to N sequence files at a time.
  --resume-downloads                  Keep existing sequence files and resume partial downloads.
  --refseq-category <all|representative|reference>
                                      Select all, representative, or reference RefSeq assemblies.
  --assembly-level <level|all>        Select a RefSeq assembly_level (default: Complete Genome).
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

rank_value() {
	case "$1" in
		--species) echo 0 ;;
		--genus) echo 1 ;;
		--family) echo 2 ;;
		--order) echo 3 ;;
		--class) echo 4 ;;
		--phylum) echo 5 ;;
		*) return 1 ;;
	esac
}

if [ "$#" -lt 2 ]; then
	usage
	exit 1
fi

DBDR_INPUT="$1"
shift
RANK=0
DATABASES=()
DOWNLOAD_THREADS="${CLARK_REFSEQ_THREADS:-1}"
RESUME_DOWNLOADS="${CLARK_REFSEQ_RESUME:-0}"
REFSEQ_CATEGORY="${CLARK_REFSEQ_CATEGORY:-all}"
ASSEMBLY_LEVEL="${CLARK_REFSEQ_ASSEMBLY_LEVEL:-Complete Genome}"

while [ "$#" -gt 0 ]; do
	case "$1" in
		--species|--genus|--family|--order|--class|--phylum)
			RANK="$(rank_value "$1")"
			shift
			;;
		--download-threads)
			[ "$#" -ge 2 ] || die "--download-threads requires a positive integer"
			DOWNLOAD_THREADS="$2"
			shift 2
			;;
		--resume-downloads)
			RESUME_DOWNLOADS=1
			shift
			;;
		--refseq-category)
			[ "$#" -ge 2 ] || die "--refseq-category requires all, representative, or reference"
			REFSEQ_CATEGORY="$2"
			shift 2
			;;
		--assembly-level)
			[ "$#" -ge 2 ] || die "--assembly-level requires a value"
			ASSEMBLY_LEVEL="$2"
			shift 2
			;;
		--*)
			die "unrecognized option '$1'"
			;;
		*)
			DATABASES+=("$1")
			shift
			;;
	esac
done

[ "${#DATABASES[@]}" -gt 0 ] || die "choose at least one database"

case "$DOWNLOAD_THREADS" in
	''|*[!0-9]*) die "--download-threads must be a positive integer" ;;
esac
[ "$DOWNLOAD_THREADS" -gt 0 ] || die "--download-threads must be a positive integer"

case "$REFSEQ_CATEGORY" in
	all|representative|reference) ;;
	*) die "--refseq-category must be all, representative, or reference" ;;
esac

export CLARK_REFSEQ_THREADS="$DOWNLOAD_THREADS"
export CLARK_REFSEQ_RESUME="$RESUME_DOWNLOADS"
export CLARK_REFSEQ_CATEGORY="$REFSEQ_CATEGORY"
export CLARK_REFSEQ_ASSEMBLY_LEVEL="$ASSEMBLY_LEVEL"

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"
mkdir -p "$DBDR_INPUT" || die "failed to create database directory '$DBDR_INPUT'"
DBDR="$(cd -P "$DBDR_INPUT" >/dev/null 2>&1 && pwd)" || die "failed to resolve database directory '$DBDR_INPUT'"

echo "$DBDR" > "$LDIR/.DBDirectory"
: > "$DBDR/targets.txt"
rm -f "$DBDR/.tmp" "$LDIR/.settings" "$LDIR/files_excluded.txt" "$DBDR/files_excluded.txt"

subDB=""
for db in "${DATABASES[@]}"; do
	echo -n "Collecting metadata of $db... "
	"$LDIR/scripts/make_metadata.sh" "$db" "$DBDR"
	[ -s "$DBDR/.$db" ] || die "metadata list for '$db' was not created"
	[ -f "$DBDR/.taxondata" ] || die "taxonomy data is missing in '$DBDR'"
	echo "done."

	if [ -s "$DBDR/.$db.fileToTaxIDs" ]; then
		"$LDIR/exe/getTargetsDef" "$DBDR/.$db.fileToTaxIDs" "$RANK" >> "$DBDR/targets.txt"
		subDB="${subDB}${db}_"
		if [ -f "$LDIR/files_excluded.txt" ]; then
			cat "$LDIR/files_excluded.txt" >> "$DBDR/.tmp"
			rm -f "$LDIR/files_excluded.txt"
		fi
	fi
done

[ -s "$DBDR/targets.txt" ] || die "no targets were generated"

subDB="${subDB}${RANK}"
{
	printf -- '-T %s\n' "$DBDR/targets.txt"
	printf -- '-D %s\n' "$DBDR/$subDB/"
} > "$LDIR/.settings"

if [ ! -d "$DBDR/$subDB" ]; then
	echo "Creating directory to store discriminative k-mers: $DBDR/$subDB"
	mkdir -p "$DBDR/$subDB"
fi

if [ -s "$DBDR/.tmp" ]; then
	mv "$DBDR/.tmp" "$DBDR/files_excluded.txt"
fi

echo "$DBDR/$subDB" > "$LDIR/.dbAddress"
echo "Targets configured in $DBDR/targets.txt"
