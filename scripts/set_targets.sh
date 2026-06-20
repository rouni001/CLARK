#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: ./set_targets.sh <Directory_path> <database choice...> [taxonomy rank]

Database choices:
  bacteria viruses plasmid plastid protozoa fungi human custom

Taxonomy rank:
  --species (default), --genus, --family, --order, --class, --phylum
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

DBDR="$1"
shift
RANK=0
DATABASES=()

for arg in "$@"; do
	case "$arg" in
		--species|--genus|--family|--order|--class|--phylum)
			RANK="$(rank_value "$arg")"
			;;
		--*)
			die "unrecognized taxonomy rank '$arg'"
			;;
		*)
			DATABASES+=("$arg")
			;;
	esac
done

[ "${#DATABASES[@]}" -gt 0 ] || die "choose at least one database"

SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"
mkdir -p "$DBDR"

echo "$DBDR" > "$LDIR/.DBDirectory"
: > "$DBDR/targets.txt"
rm -f "$DBDR/.tmp" "$LDIR/.settings" "$LDIR/files_excluded.txt" "$DBDR/files_excluded.txt"

subDB=""
for db in "${DATABASES[@]}"; do
	echo -n "Collecting metadata of $db... "
	"$LDIR/make_metadata.sh" "$db" "$DBDR"
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
