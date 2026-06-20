#!/usr/bin/env bash

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage: scripts/make_metadata.sh <database name> <database directory>

Supported database names:
  bacteria viruses plasmid plastid protozoa fungi human custom
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

supported_database() {
	case "$1" in
		bacteria|viruses|plasmid|plastid|protozoa|fungi|human|custom) return 0 ;;
		*) return 1 ;;
	esac
}

if [ "$#" -lt 2 ]; then
	usage
	exit 1
fi

DB="$1"
DBDR="$2"
TAXDR="taxonomy"
SCRIPT_DIR="$(script_dir)"
LDIR="${CLARK_HOME:-$(cd "$SCRIPT_DIR/.." >/dev/null 2>&1 && pwd)}"

supported_database "$DB" || die "unsupported database '$DB'. Supported: bacteria, viruses, plasmid, plastid, protozoa, fungi, human, custom."

mkdir -p "$DBDR/Custom"

if [ ! -d "$DBDR/$TAXDR" ]; then
	echo "Taxonomy data missing. The program will download data to $DBDR/$TAXDR."
	mkdir -p "$DBDR/$TAXDR"
	"$LDIR/scripts/download_taxondata.sh" "$DBDR/$TAXDR"
fi

if [ ! -f "$DBDR/.taxondata" ]; then
	echo "Failed to find taxonomy files. The program will try to download them..."
	"$LDIR/scripts/download_taxondata.sh" "$DBDR/$TAXDR"
	[ -f "$DBDR/.taxondata" ] || die "failed to find taxonomy files"
fi

if [ ! -s "$DBDR/.$DB" ]; then
	if [ "$DB" != "custom" ]; then
		echo "Sequences for $DB not found. The program will download them."
		"$LDIR/scripts/download_RefSeqDB.sh" "$DBDR" "$DB"
	else
		find "$DBDR/Custom" -type f -name '*.f*' > "$DBDR/.$DB"
		if [ ! -s "$DBDR/.$DB" ]; then
			die "the custom database directory '$DBDR/Custom' is empty. Add FASTA files, then rerun with 'custom'."
		fi
	fi
fi

if [ ! -x "$LDIR/exe/getfilesToTaxNodes" ] || [ ! -x "$LDIR/exe/getAccssnTaxID" ]; then
	die "required helper executables are missing. Run scripts/install.sh first."
fi

[ -s "$DBDR/.$DB" ] || die "failed to find $DB sequences"

if [ "$DB" = "human" ]; then
	if [ ! -s "$DBDR/.$DB" ]; then
		find "$DBDR/Human" -type f > "$DBDR/.$DB"
	fi
	if [ ! -s "$DBDR/.$DB.fileToTaxIDs" ]; then
		while IFS= read -r file || [ -n "$file" ]; do
			[ -n "$file" ] || continue
			printf '%s X 9606 9605 9604 9443 40674 7711\n' "$file"
		done < "$DBDR/.$DB" > "$DBDR/.$DB.fileToTaxIDs"
	fi
	exit 0
fi

if [ ! -s "$DBDR/.$DB.fileToAccssnTaxID" ]; then
	echo "Re-building $DB.fileToAccssnTaxID"
	"$LDIR/exe/getAccssnTaxID" "$DBDR/.$DB" "$DBDR/$TAXDR/nucl_accss" "$DBDR/$TAXDR/merged.dmp" > "$DBDR/.$DB.fileToAccssnTaxID"
fi

if [ ! -s "$DBDR/.$DB.fileToTaxIDs" ]; then
	echo "$DB: Retrieving taxonomy nodes for each sequence based on taxon ID..."
	"$LDIR/exe/getfilesToTaxNodes" "$DBDR/$TAXDR/nodes.dmp" "$DBDR/.$DB.fileToAccssnTaxID" > "$DBDR/.$DB.fileToTaxIDs"
fi
