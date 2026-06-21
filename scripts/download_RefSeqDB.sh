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
	cat <<USAGE
Usage: $0 [options] <Directory for the sequences> <Database: bacteria, viruses, plasmid, plastid, protozoa, fungi or human>

Options:
  --dry-run                         Plan downloads without fetching sequence files.
  --threads <N>                     Download up to N sequence files at a time (default: 1).
  --resume                          Keep existing sequence files and resume partial downloads.
  --refseq-category <all|representative|reference>
                                    Filter RefSeq assembly summaries by refseq_category (default: all).
  --assembly-level <level|all>      Filter RefSeq assembly summaries by assembly_level (default: Complete Genome).
USAGE
}

die() {
	echo "Error: $*" >&2
	exit 1
}

DRY_RUN=${CLARK_REFSEQ_DRY_RUN:-0}
THREADS=${CLARK_REFSEQ_THREADS:-1}
RESUME=${CLARK_REFSEQ_RESUME:-0}
REFSEQ_CATEGORY=${CLARK_REFSEQ_CATEGORY:-all}
ASSEMBLY_LEVEL=${CLARK_REFSEQ_ASSEMBLY_LEVEL:-Complete Genome}

while [ "$#" -gt 0 ]; do
	case "$1" in
		--dry-run)
			DRY_RUN=1
			shift
			;;
		--threads)
			[ "$#" -ge 2 ] || die "--threads requires a positive integer"
			THREADS="$2"
			shift 2
			;;
		--resume)
			RESUME=1
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
			die "unrecognized option: $1"
			;;
		*)
			break
			;;
	esac
done

if [ "$#" -ne 2 ]; then
	usage
	exit 1
fi

case "$THREADS" in
	''|*[!0-9]*) die "--threads must be a positive integer" ;;
esac
[ "$THREADS" -gt 0 ] || die "--threads must be a positive integer"

case "$REFSEQ_CATEGORY" in
	all|representative|reference) ;;
	*) die "--refseq-category must be all, representative, or reference" ;;
esac

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
DOWNLOAD_LIST="$DBDR/.$DB.download_urls.tsv"

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
	printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n' > "$PROVENANCE"
	: > "$DOWNLOAD_LIST"
}

record_manifest_line() {
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$RUN_STARTED_UTC" "$DB" "$2" "$3" "$4" "$5" "$6" >> "$1"
}

record_manifest() {
	record_manifest_line "$MANIFEST" "$1" "$2" "$3" "$4" "$5"
}

record_provenance() {
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$DB" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$PROVENANCE"
}

append_download_url() {
	printf '%s\t%s\n' "$1" "$2" >> "$DOWNLOAD_LIST"
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

	if command -v wget >/dev/null 2>&1; then
		wget -O "$output" "$url"
	elif command -v curl >/dev/null 2>&1; then
		curl -fL -o "$output" "$url"
	else
		die "neither wget nor curl is available"
	fi

	[ -s "$output" ] || die "failed to download $url"
	record_manifest "metadata" "$source" "$url" "$output" "downloaded"
}

fetch_url_to_status() {
	url="$1"
	source="$2"
	status_file="$3"
	output=${url##*/}

	if [ -s "$output" ]; then
		record_manifest_line "$status_file" "download" "$source" "$url" "$(pwd)/$output" "skipped-existing"
		return 0
	fi
	if [ "${output%.gz}" != "$output" ] && [ -s "${output%.gz}" ]; then
		record_manifest_line "$status_file" "download" "$source" "$url" "$(pwd)/${output%.gz}" "skipped-existing"
		return 0
	fi

	partial="$output.part"
	if [ "$RESUME" != "1" ]; then
		rm -f "$partial"
	fi

	if command -v wget >/dev/null 2>&1; then
		if [ "$RESUME" = "1" ]; then
			wget -c -O "$partial" "$url"
		else
			wget -O "$partial" "$url"
		fi
	elif command -v curl >/dev/null 2>&1; then
		if [ "$RESUME" = "1" ]; then
			curl -fL -C - -o "$partial" "$url"
		else
			curl -fL -o "$partial" "$url"
		fi
	else
		die "neither wget nor curl is available"
	fi

	[ -s "$partial" ] || return 1
	mv "$partial" "$output"
	record_manifest_line "$status_file" "download" "$source" "$url" "$(pwd)/$output" "downloaded"
}

download_assembly_source() {
	source="$1"
	summary_url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/$source/assembly_summary.txt"
	summary_file="$DBDR/.$DB.assembly_summary.$source.txt"

	fetch_to_file "$summary_url" "$summary_file" "$source"

	awk -F '\t' \
		-v db="$DB" \
		-v source="$source" \
		-v provenance="$PROVENANCE" \
		-v category="$REFSEQ_CATEGORY" \
		-v assembly_level="$ASSEMBLY_LEVEL" \
		-v download_list="$DOWNLOAD_LIST" '
			BEGIN { OFS = "\t" }
			$0 !~ /^#/ && $11 == "latest" && $20 != "" {
				if (assembly_level != "all" && $12 != assembly_level) {
					next
				}
				if (category == "representative" && $5 != "representative genome") {
					next
				}
				if (category == "reference" && $5 != "reference genome") {
					next
				}
				n = split($20, path_parts, "/")
				url = $20 "/" path_parts[n] "_genomic.fna.gz"
				print db, source, $1, $6, $7, $15, $12, $11, url >> provenance
				print source, url >> download_list
			}
		' "$summary_file"

	rm -f "$summary_file"
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
		record_provenance "$source" "$(accession_from_url "$genome_url")" "NA" "NA" "NA" "static" "latest" "$genome_url"
		append_download_url "$source" "$genome_url"
	done
}

decompress_gz_files() {
	find "$(pwd)" -type f -name '*.gz' -print | while IFS= read -r gz_file || [ -n "$gz_file" ]; do
		[ -n "$gz_file" ] || continue
		gunzip -f "$gz_file"
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

download_url_list() {
	total=$(wc -l < "$DOWNLOAD_LIST" | tr -d ' ')
	echo "Selected $total $DB RefSeq genome(s) using assembly_level=$ASSEMBLY_LEVEL and refseq_category=$REFSEQ_CATEGORY."

	if [ "$DRY_RUN" = "1" ]; then
		while IFS='	' read -r source genome_url || [ -n "$genome_url" ]; do
			[ -n "$genome_url" ] || continue
			record_manifest "plan" "$source" "$genome_url" "$DATA_DIR/${genome_url##*/}" "dry-run"
		done < "$DOWNLOAD_LIST"
		return 0
	fi

	status_dir="$DBDR/.$DB.download-status.$$"
	rm -rf "$status_dir"
	mkdir -p "$status_dir"
	job_count=0
	job_index=0
	while IFS='	' read -r source genome_url || [ -n "$genome_url" ]; do
		[ -n "$genome_url" ] || continue
		job_index=$((job_index + 1))
		fetch_url_to_status "$genome_url" "$source" "$status_dir/$job_index.tsv" &
		job_count=$((job_count + 1))
		if [ "$job_count" -ge "$THREADS" ]; then
			wait || die "failed to download one or more RefSeq genomes"
			job_count=0
		fi
	done < "$DOWNLOAD_LIST"
	wait || die "failed to download one or more RefSeq genomes"

	if [ -n "$(find "$status_dir" -type f -print -quit)" ]; then
		cat "$status_dir"/*.tsv >> "$MANIFEST"
	fi
	rm -rf "$status_dir"
}

if [ "$DRY_RUN" != "1" ] && [ -s "$MARKER" ]; then
	echo "$(display_name "$DB") sequences already in $DBDR."
	exit 0
fi

DATA_DIR="$DBDR/$(db_directory_name "$DB")"
if [ "$DRY_RUN" = "1" ]; then
	echo "Dry run: planning $(display_name "$DB") RefSeq downloads."
else
	if [ "$RESUME" = "1" ]; then
		rm -f "$DBDR"/".$DB".download_manifest.tsv "$DBDR"/".$DB".provenance.tsv "$DBDR"/".$DB".download_urls.tsv
		mkdir -p "$DATA_DIR"
	else
		rm -Rf "$DATA_DIR" "$DBDR"/".$DB".*
		mkdir -m 775 "$DATA_DIR"
	fi
	cd "$DATA_DIR" || exit 1
	echo "Downloading now $(display_name "$DB") RefSeq genomes."
	echo "Download threads: $THREADS; resume: $RESUME; RefSeq category: $REFSEQ_CATEGORY; assembly level: $ASSEMBLY_LEVEL."
fi

init_reports

for source in $(assembly_sources); do
	download_assembly_source "$source"
done
download_static_urls
download_url_list

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
