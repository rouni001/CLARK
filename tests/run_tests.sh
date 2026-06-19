#!/usr/bin/env bash

set -euo pipefail

REPO_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"

fail() {
	echo "FAIL: $*" >&2
	exit 1
}

pass() {
	echo "PASS: $*"
}

require_file() {
	[ -e "$1" ] || fail "missing $1"
}

test_shell_syntax() {
	for script in \
		"$REPO_DIR/install.sh" \
		"$REPO_DIR/classify_metagenome.sh" \
		"$REPO_DIR/set_targets.sh" \
		"$REPO_DIR/make_metadata.sh" \
		"$REPO_DIR/download_taxondata.sh"; do
		bash -n "$script"
	done
	pass "modern shell entrypoints parse"
}

test_version_binary() {
	require_file "$REPO_DIR/exe/CLARK"
	version="$("$REPO_DIR/exe/CLARK" --version)"
	case "$version" in
		*"Version: 1.3.0.0"*) pass "CLARK version executable runs" ;;
		*) fail "unexpected version output: $version" ;;
	esac
}

create_fake_exe() {
	local dir="$1"
	mkdir -p "$dir"
	for binary in CLARK CLARK-l CLARK-S; do
		cat > "$dir/$binary" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
: "${CLARK_TEST_CAPTURE:?CLARK_TEST_CAPTURE is required}"
{
	echo "binary=$(basename "$0")"
	while [ "$#" -gt 0 ]; do
		echo "arg=<$1>"
		if [ "$1" = "-O" ]; then
			shift
			echo "arg=<$1>"
			if [ -f "$1" ]; then
				echo "object-bytes=$(wc -c < "$1" | tr -d ' ')"
			else
				echo "object-missing"
			fi
		fi
		shift || true
	done
} > "$CLARK_TEST_CAPTURE"
STUB
		chmod +x "$dir/$binary"
	done
}

test_classify_wrapper_quotes_paths() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings file"
	dbdir="$tmp/db dir"
	targets="$tmp/targets file.txt"
	input="$tmp/input reads.fastq"
	result="$tmp/result file.csv"
	fake_exe="$tmp/fake exe"
	capture="$tmp/capture.txt"

	mkdir -p "$dbdir"
	printf 'target.fa 12345\n' > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir/" > "$settings"
	printf '@r1\nACGT\n+\n!!!!\n' > "$input"
	create_fake_exe "$fake_exe"

	CLARK_SETTINGS_FILE="$settings" \
	CLARK_EXE_DIR="$fake_exe" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/classify_metagenome.sh" -O "$input" -R "$result" -m 2 -n 4

	grep -Fq "arg=<$targets>" "$capture" || fail "target path with spaces was not preserved"
	grep -Fq "arg=<$dbdir/>" "$capture" || fail "database path with spaces was not preserved"
	grep -Fq "arg=<$input>" "$capture" || fail "input path with spaces was not preserved"
	grep -Fq "arg=<$result>" "$capture" || fail "result path with spaces was not preserved"
	grep -Fq "object-bytes=16" "$capture" || fail "input file was not passed to classifier"
	pass "classify wrapper preserves paths with spaces"
}

test_classify_wrapper_gzip() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-gzip-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	input="$tmp/input reads.fastq.gz"
	result="$tmp/result.csv"
	fake_exe="$tmp/fake-exe"
	capture="$tmp/capture.txt"

	mkdir -p "$dbdir"
	printf 'target.fa 12345\n' > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir/" > "$settings"
	printf '@r1\nACGT\n+\n!!!!\n' | gzip > "$input"
	create_fake_exe "$fake_exe"

	TMPDIR="$tmp" \
	CLARK_SETTINGS_FILE="$settings" \
	CLARK_EXE_DIR="$fake_exe" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/classify_metagenome.sh" -O "$input" -R "$result" --gzipped

	grep -Fq "object-bytes=16" "$capture" || fail "gzipped input was not decompressed before classification"
	if find "$tmp" -maxdepth 1 -type d -name 'CLARKGZP.*' | grep -q .; then
		fail "temporary gzipped-input directory was not cleaned up"
	fi
	pass "classify wrapper decompresses gzipped input"
}

test_clark_l_label_bug_regression() {
	block="$(awk '
		/for\(size_t t = 0 ; t < m_labels.size\(\); t\+\+\)/ { capture=1 }
		capture { print }
		capture && /_filesHT.push_back\(fname\);/ { exit }
	' "$REPO_DIR/src/CLARK_hh.hh")"

	printf '%s\n' "$block" | grep -Fq 'm_labels[t]' || fail "regular label loop does not use m_labels[t]"
	if printf '%s\n' "$block" | grep -Fq 'm_labels_c[t]'; then
		fail "regular label loop still references m_labels_c[t]"
	fi
	pass "CLARK-l target filename regression is fixed"
}

test_ncbi_urls() {
	if grep -R "ftp://ftp.ncbi.nih.gov\|ftp://ftp.ncbi.nlm.nih.gov" "$REPO_DIR/download_RefSeqDB.sh" "$REPO_DIR/download_taxondata.sh" >/dev/null; then
		fail "legacy or misspelled NCBI FTP URL remains"
	fi
	grep -Fq "https://ftp.ncbi.nlm.nih.gov" "$REPO_DIR/download_taxondata.sh" || fail "taxonomy downloader does not use HTTPS NCBI URL"
	pass "NCBI download URLs use HTTPS host"
}

test_portable_script_paths() {
	if grep -R "[r]eadlink -f" "$REPO_DIR" --include='*.sh' --exclude-dir='.git' --exclude-dir='build' --exclude-dir='exe' >/dev/null; then
		fail "non-portable shell path resolver remains"
	fi
	pass "shell scripts use portable path resolution"
}

test_shell_syntax
test_version_binary
test_classify_wrapper_quotes_paths
test_classify_wrapper_gzip
test_clark_l_label_bug_regression
test_ncbi_urls
test_portable_script_paths
