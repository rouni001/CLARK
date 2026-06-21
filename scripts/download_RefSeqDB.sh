#!/bin/sh

# 
#   CLARK, CLAssifier based on Reduced K-mers.
# 
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#   Copyright @ The Regents of the University of California. All rights reserved.
#
#   download_RefseqDB.sh: To download complete reference genomes of NCBI/RefSeq
#   (i.e., Bacteria/Archaea, Viruses, Plasmid, Plastid, Protozoa, Fungi and Human).

set -eu

usage() {
	echo "Usage: $0 [--dry-run] <Directory for the sequences> <Database: bacteria, viruses, plasmid, plastid, protozoa, fungi or human> "
}

die() {
	echo "Error: $*" >&2
	exit 1
}

if [ "${1:-}" = "--dry-run" ]; then
	DRY_RUN=1
	shift
else
	DRY_RUN=${CLARK_REFSEQ_DRY_RUN:-0}
fi

if [ "$#" -ne 2 ]; then
	usage
	exit 1
fi

DIR=${CLARK_HOME:-$(CDPATH= cd "$(dirname "$0")/.." && pwd -P)}
DBDR="$1"
DB="$2"
STATIC_URLS_FILE=${CLARK_REFSEQ_STATIC_URLS:-"$DIR/scripts/refseq_static_urls.tsv"}
RUN_STARTED_UTC=$(date -u '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date '+%Y-%m-%dT%H:%M:%SZ')

case "$DB" in
	bacteria|viruses|plasmid|plastid|protozoa|fungi|human) ;;
	*) die "failed to recognize parameter: $DB. Please choose between: bacteria, viruses, plasmid, plastid, protozoa, fungi or human." ;;
esac

mkdir -p "$DBDR"
DBDR=$(CDPATH= cd "$DBDR" && pwd -P)
MARKER="$DBDR/.$DB"
MANIFEST="$DBDR/.$DB.download_manifest.tsv"
PROVENANCE="$DBDR/.$DB.provenance.tsv"

db_directory_name() {
	case "$1" in
		bacteria) echo "Bacteria" ;;
		viruses) echo "Viruses" ;;
		plasmid) echo "Plasmid" ;;
		plastid) echo "Plastid" ;;
		protozoa) echo "Protozoa" ;;
		fungi) echo "Fungi" ;;
		human) echo "Human" ;;
	esac
}

display_name() {
	case "$1" in
		bacteria) echo "Bacteria/Archaea" ;;
		viruses) echo "Viruses" ;;
		plasmid) echo "Plasmid" ;;
		plastid) echo "Plastid" ;;
		protozoa) echo "Protozoa" ;;
		fungi) echo "Fungi" ;;
		human) echo "Human" ;;
	esac
}

assembly_sources() {
	case "$DB" in
		bacteria)
			printf '%s\n' bacteria archaea
			;;
		viruses)
			printf '%s\n' viral
			;;
		protozoa|fungi)
			printf '%s\n' "$DB"
			;;
	esac
}

sequence_pattern() {
	case "$DB" in
		plasmid|plastid) echo "*.fa" ;;
		*) echo "*.fna" ;;
	esac
}

needs_sequence_split() {
	case "$DB" in
		plasmid|plastid) return 0 ;;
		*) return 1 ;;
	esac
}

init_reports() {
	printf 'timestamp_utc\tdatabase\taction\tsource\turl\tlocal_path\tstatus\n' > "$MANIFEST"
	printf 'database\tsource\taccession\tseq_rel_date\tassembly_level\tversion_status\turl\n' > "$PROVENANCE"
}

record_manifest() {
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$RUN_STARTED_UTC" "$DB" "$1" "$2" "$3" "$4" "$5" >> "$MANIFEST"
}

record_provenance() {
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$DB" "$1" "$2" "$3" "$4" "$5" "$6" >> "$PROVENANCE"
}

accession_from_url() {
	file_name=${1##*/}
	accession=${file_name%_genomic.fna.gz}
	accession=${accession%.genomic.fna.gz}
	accession=${accession%.fna.gz}
	echo "$accession"
}

fetch_to_file() {
	url="$1"
	output="$2"
	source="$3"

	if [ "$DRY_RUN" = "1" ]; then
		record_manifest "plan" "$source" "$url" "$output" "dry-run"
		return 0
	fi

	if command -v wget >/dev/null 2>&1; then
		wget -O "$output" "$url"
	elif command -v curl >/dev/null 2>&1; then
		curl -fL -o "$output" "$url"
	else
		die "neither wget nor curl is available"
	fi

	[ -s "$output" ] || die "failed to download $url"
	record_manifest "download" "$source" "$url" "$(pwd)/$output" "downloaded"
}

fetch_url() {
	url="$1"
	source="$2"
	output=${url##*/}

	if [ "$DRY_RUN" = "1" ]; then
		record_manifest "plan" "$source" "$url" "$output" "dry-run"
		return 0
	fi

	if command -v wget >/dev/null 2>&1; then
		wget "$url"
	elif command -v curl >/dev/null 2>&1; then
		curl -fLO "$url"
	else
		die "neither wget nor curl is available"
	fi

	[ -s "$output" ] || die "failed to download $url"
	record_manifest "download" "$source" "$url" "$(pwd)/$output" "downloaded"
}

download_assembly_source() {
	source="$1"
	summary_url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/$source/assembly_summary.txt"
	summary_file="assembly_summary.$source.txt"
	urls_file=".$DB.$source.urls"

	fetch_to_file "$summary_url" "$summary_file" "$source"
	if [ "$DRY_RUN" = "1" ]; then
		record_provenance "$source" "assembly_summary" "NA" "Complete Genome" "latest" "$summary_url"
		return 0
	fi

	awk -F '\t' -v db="$DB" -v source="$source" -v provenance="$PROVENANCE" '
		BEGIN { OFS = "\t" }
		$12 == "Complete Genome" && $11 == "latest" && $20 != "" {
			n = split($20, path_parts, "/")
			url = $20 "/" path_parts[n] "_genomic.fna.gz"
			print db, source, $1, $15, $12, $11, url >> provenance
			print url
		}
	' "$summary_file" > "$urls_file"

	while IFS= read -r genome_url || [ -n "$genome_url" ]; do
		[ -n "$genome_url" ] || continue
		fetch_url "$genome_url" "$source"
	done < "$urls_file"

	rm -f "$summary_file" "$urls_file"
}

download_static_urls() {
	if [ ! -f "$STATIC_URLS_FILE" ]; then
		case "$DB" in
			plasmid|plastid|fungi|human)
				die "missing static URL manifest: $STATIC_URLS_FILE"
				;;
			*)
				return 0
				;;
		esac
	fi

	awk -v db="$DB" '
		BEGIN { FS = "\t" }
		$0 !~ /^#/ && NF >= 3 && $1 == db { print $2 " " $3 }
	' "$STATIC_URLS_FILE" | while read -r source genome_url; do
		[ -n "$genome_url" ] || continue
		record_provenance "$source" "$(accession_from_url "$genome_url")" "NA" "static" "latest" "$genome_url"
		fetch_url "$genome_url" "$source"
	done
}

decompress_gz_files() {
	find "$(pwd)" -type f -name '*.gz' -print | while IFS= read -r gz_file || [ -n "$gz_file" ]; do
		[ -n "$gz_file" ] || continue
		gunzip "$gz_file"
	done
}

split_fna_records() {
	find "$(pwd)" -type f -name '*.fna' -print | while IFS= read -r fasta_file || [ -n "$fasta_file" ]; do
		[ -n "$fasta_file" ] || continue
		"$DIR/exe/exeSeq" "$fasta_file" ./
	done
}

write_sequence_marker() {
	find "$(pwd)" -name "$(sequence_pattern)" > "$MARKER"
}

if [ "$DRY_RUN" != "1" ] && [ -s "$MARKER" ]; then
	echo "$(display_name "$DB") sequences already in $DBDR."
	exit 0
fi

DATA_DIR="$DBDR/$(db_directory_name "$DB")"
if [ "$DRY_RUN" = "1" ]; then
	echo "Dry run: planning $(display_name "$DB") RefSeq downloads."
else
	rm -Rf "$DATA_DIR" "$DBDR"/".$DB".*
	mkdir -m 775 "$DATA_DIR"
	cd "$DATA_DIR" || exit 1
	echo "Downloading now $(display_name "$DB") RefSeq genomes."
fi

init_reports

for source in $(assembly_sources); do
	download_assembly_source "$source"
done
download_static_urls

if [ "$DRY_RUN" = "1" ]; then
	echo "Dry run complete."
	echo "Download manifest: $MANIFEST"
	echo "Provenance: $PROVENANCE"
	exit 0
fi

echo "Downloading done. Uncompressing files..."
decompress_gz_files

if needs_sequence_split; then
	echo "Processing sequences..."
	split_fna_records
fi

write_sequence_marker
[ -s "$MARKER" ] || die "failed to download $(display_name "$DB") sequences"
echo "$(display_name "$DB") sequences downloaded!"
