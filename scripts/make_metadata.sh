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

build_refseq_taxid_map_from_provenance() {
	local provenance="$DBDR/.$DB.provenance.tsv"
	[ -s "$provenance" ] || return 1

	case "$DB" in
		bacteria|viruses|protozoa|fungi) ;;
		*) return 1 ;;
	esac

	awk -F '\t' -v provenance="$provenance" '
		BEGIN {
			while ((getline line < provenance) > 0) {
				ncols = split(line, cols, FS)
				if (line ~ /^database\t/) {
					for (i = 1; i <= ncols; i++) {
						header[cols[i]] = i
					}
					continue
				}
				if (!("accession" in header) || !("taxid" in header) || !("url" in header)) {
					exit 2
				}
				taxid = cols[header["taxid"]]
				if (taxid !~ /^[0-9]+$/ || taxid == "0") {
					continue
				}
				file_name = cols[header["url"]]
				sub(/^.*\//, "", file_name)
				sub(/\.gz$/, "", file_name)
				taxid_by_file[file_name] = cols[header["accession"]] "\t" taxid
			}
			close(provenance)
		}
		{
			file_name = $0
			sub(/^.*\//, "", file_name)
			if (file_name in taxid_by_file) {
				print $0 "\t" taxid_by_file[file_name]
				matched++
			} else {
				missing++
			}
		}
		END {
			if (matched == 0 || missing > 0) {
				exit 1
			}
		}
	' "$DBDR/.$DB"
}

ensure_taxonomy_data() {
	local require_accession_maps="$1"
	local downloader_args=()
	if [ "$require_accession_maps" != "1" ]; then
		downloader_args+=(--skip-accession-maps)
	fi
	if [ "${CLARK_REFSEQ_INSECURE_TLS:-0}" = "1" ]; then
		downloader_args+=(--insecure-tls)
	fi

	if [ ! -d "$DBDR/$TAXDR" ]; then
		echo "Taxonomy data missing. The program will download data to $DBDR/$TAXDR."
		mkdir -p "$DBDR/$TAXDR"
		"$LDIR/scripts/download_taxondata.sh" "${downloader_args[@]}" "$DBDR/$TAXDR"
	fi

	if [ ! -f "$DBDR/.taxondata" ] || { [ "$require_accession_maps" = "1" ] && [ ! -s "$DBDR/$TAXDR/nucl_accss" ]; }; then
		echo "Failed to find required taxonomy files. The program will try to download them..."
		"$LDIR/scripts/download_taxondata.sh" "${downloader_args[@]}" "$DBDR/$TAXDR"
	fi

	[ -f "$DBDR/.taxondata" ] || die "failed to find taxonomy files"
	if [ "$require_accession_maps" = "1" ]; then
		[ -s "$DBDR/$TAXDR/nucl_accss" ] || die "failed to find accession-to-taxid maps"
	fi
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
	die "required helper executables are missing. Run 'make all' from the CLARK repository root first."
fi

[ -s "$DBDR/.$DB" ] || die "failed to find $DB sequences"

if [ "$DB" != "human" ] && [ ! -s "$DBDR/.$DB.fileToAccssnTaxID" ]; then
	tmp_accss="$DBDR/.$DB.fileToAccssnTaxID.tmp.$$"
	if build_refseq_taxid_map_from_provenance > "$tmp_accss"; then
		mv "$tmp_accss" "$DBDR/.$DB.fileToAccssnTaxID"
		echo "Built $DB.fileToAccssnTaxID from RefSeq provenance."
	else
		rm -f "$tmp_accss"
	fi
fi

require_accession_maps=0
if [ "$DB" != "human" ] && [ ! -s "$DBDR/.$DB.fileToAccssnTaxID" ]; then
	require_accession_maps=1
fi
ensure_taxonomy_data "$require_accession_maps"

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
