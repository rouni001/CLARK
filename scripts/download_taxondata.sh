#!/usr/bin/env bash

set -euo pipefail

usage() {
	echo "Usage: $0 <Directory: directory to store taxonomy data>"
}

die() {
	echo "Error: $*" >&2
	exit 1
}

download_file() {
	local url="$1"
	local dest="$2"

	if command -v curl >/dev/null 2>&1; then
		curl -L --fail --retry 3 --output "$dest" "$url"
	elif command -v wget >/dev/null 2>&1; then
		wget -O "$dest" "$url"
	else
		die "curl or wget is required to download taxonomy data"
	fi
}

if [ "$#" -ne 1 ]; then
	usage
	exit 1
fi

TARGET_DIR="$1"
case "$TARGET_DIR" in
	""|"/")
		die "refusing to write taxonomy data to '$TARGET_DIR'"
		;;
esac

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

echo "Downloading NCBI taxonomy data..."
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz" "nucl_gb.accession2taxid.gz"
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz" "nucl_wgs.accession2taxid.gz"
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "taxdump.tar.gz"

if [ -s "nucl_gb.accession2taxid.gz" ] && [ -s "nucl_wgs.accession2taxid.gz" ] && [ -s "taxdump.tar.gz" ]; then
	echo "Uncompressing taxonomy data..."
	gunzip -f "nucl_wgs.accession2taxid.gz"
	gunzip -f "nucl_gb.accession2taxid.gz"
	tar -zxf "taxdump.tar.gz"
else
	die "failed to download one or more taxonomy files"
fi

if [ -s "nucl_gb.accession2taxid" ] && [ -s "nucl_wgs.accession2taxid" ] && [ -s "nodes.dmp" ] && [ -s "merged.dmp" ] && [ -s "names.dmp" ]; then
	cat "nucl_gb.accession2taxid" "nucl_wgs.accession2taxid" > "nucl_accss"
	touch "../.taxondata"
	echo "Taxonomy data ready in $TARGET_DIR"
else
	die "failed to unpack complete taxonomy data"
fi
