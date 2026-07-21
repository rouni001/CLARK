#!/usr/bin/env bash

set -euo pipefail

usage() {
	echo "Usage: $0 [--skip-accession-maps] [--insecure-tls] <Directory: directory to store taxonomy data>"
}

die() {
	echo "Error: $*" >&2
	exit 1
}

download_file() {
	local url="$1"
	local dest="$2"
	local curl_tls_opt=()
	local wget_tls_opt=()

	if [ "$INSECURE_TLS" = "1" ]; then
		curl_tls_opt=(-k)
		wget_tls_opt=(--no-check-certificate)
	fi

	if command -v curl >/dev/null 2>&1; then
		curl -L --fail --retry 3 "${curl_tls_opt[@]}" --output "$dest" "$url"
	elif command -v wget >/dev/null 2>&1; then
		wget "${wget_tls_opt[@]}" -O "$dest" "$url"
	else
		die "curl or wget is required to download taxonomy data"
	fi
}

SKIP_ACCESSION_MAPS=0
INSECURE_TLS="${CLARK_REFSEQ_INSECURE_TLS:-0}"

while [ "$#" -gt 0 ]; do
	case "$1" in
		--skip-accession-maps)
			SKIP_ACCESSION_MAPS=1
			shift
			;;
		--insecure-tls)
			INSECURE_TLS=1
			shift
			;;
		--*)
			die "unrecognized option: $1"
			;;
		*)
			break
			;;
	esac
done

case "$INSECURE_TLS" in
	0|1) ;;
	*) die "CLARK_REFSEQ_INSECURE_TLS must be 0 or 1" ;;
esac

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

if [ "$INSECURE_TLS" = "1" ]; then
	echo "Warning: TLS certificate verification is disabled for taxonomy downloads (--insecure-tls)."
fi
echo "Downloading NCBI taxonomy data..."
if [ "$SKIP_ACCESSION_MAPS" != "1" ]; then
	download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz" "nucl_gb.accession2taxid.gz"
	download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz" "nucl_wgs.accession2taxid.gz"
fi
download_file "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "taxdump.tar.gz"

if { [ "$SKIP_ACCESSION_MAPS" = "1" ] || { [ -s "nucl_gb.accession2taxid.gz" ] && [ -s "nucl_wgs.accession2taxid.gz" ]; }; } && [ -s "taxdump.tar.gz" ]; then
	echo "Uncompressing taxonomy data..."
	if [ "$SKIP_ACCESSION_MAPS" != "1" ]; then
		gunzip -f "nucl_wgs.accession2taxid.gz"
		gunzip -f "nucl_gb.accession2taxid.gz"
	fi
	tar -zxf "taxdump.tar.gz"
else
	die "failed to download one or more taxonomy files"
fi

if [ -s "nodes.dmp" ] && [ -s "merged.dmp" ] && [ -s "names.dmp" ]; then
	if [ "$SKIP_ACCESSION_MAPS" != "1" ]; then
		[ -s "nucl_gb.accession2taxid" ] && [ -s "nucl_wgs.accession2taxid" ] || die "failed to unpack accession-to-taxid maps"
		cat "nucl_gb.accession2taxid" "nucl_wgs.accession2taxid" > "nucl_accss"
	fi
	touch "../.taxondata"
	echo "Taxonomy data ready in $TARGET_DIR"
else
	die "failed to unpack complete taxonomy data"
fi
