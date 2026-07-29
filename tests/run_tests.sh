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

REQUIRED_EXECUTABLES=(
	CLARK
	CLARK-l
	CLARK-S
	converter
	dscriptMaker
	exeSeq
	extractSeqs
	getAbundance
	getAccssnTaxID
	getConfidenceDensity
	getGammaDensity
	getTargetSpecificKmersStat
	getTargetsDef
	getfilesToTaxNodes
	makeSummaryTables
)

test_shell_syntax() {
	for script in \
		"$REPO_DIR/scripts/classify_metagenome.sh" \
		"$REPO_DIR/scripts/set_targets.sh" \
		"$REPO_DIR/scripts/make_metadata.sh" \
		"$REPO_DIR/scripts/download_taxondata.sh" \
		"$REPO_DIR/tests/coverage_report.sh"; do
		bash -n "$script"
	done
	for script in "$REPO_DIR"/scripts/*.sh; do
		case "$(head -n 1 "$script")" in
			*"bash"*) bash -n "$script" ;;
			*) sh -n "$script" ;;
		esac
	done
	pass "modern shell entrypoints parse"
}

test_no_root_shell_scripts() {
	if find "$REPO_DIR" -maxdepth 1 -type f -name '*.sh' | grep -q .; then
		find "$REPO_DIR" -maxdepth 1 -type f -name '*.sh' >&2
		fail "root-level shell scripts remain; use scripts/ as the single script home"
	fi
	pass "repository root has no duplicated shell scripts"
}

test_required_executables() {
	for binary in "${REQUIRED_EXECUTABLES[@]}"; do
		[ -x "$REPO_DIR/exe/$binary" ] || fail "missing executable exe/$binary"
	done
	pass "all required executables are present"
}

test_version_binaries() {
	for binary in CLARK CLARK-l CLARK-S; do
		version="$("$REPO_DIR/exe/$binary" --version)"
		case "$version" in
			*"Version: 1.4.5.0-a"*) ;;
			*) fail "unexpected version output from $binary: $version" ;;
		esac
	done
	pass "CLARK variant version executables run"
}

test_file_hash_unit_binary() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-file-hash-unit-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	[ -x "$REPO_DIR/exe/unit_file_hash_tests" ] || fail "missing executable exe/unit_file_hash_tests"
	"$REPO_DIR/exe/unit_file_hash_tests" "$tmp"
}

test_kso_requires_spectrum_mode() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-kso-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	targets="$tmp/targets.txt"
	db="$tmp/db"
	objects="$tmp/objects.fa"
	result="$tmp/result.csv"
	output="$tmp/output.txt"

	printf 'target.fa 12345\n' > "$targets"
	printf '' > "$db"
	printf '>read1\nACGT\n' > "$objects"

	for binary in CLARK CLARK-l CLARK-S; do
		for args in "--kso -m 0" "-m 0 --kso"; do
			if "$REPO_DIR/exe/$binary" $args -T "$targets" -D "$db" -O "$objects" -R "$result" > "$output" 2>&1; then
				fail "$binary accepted --kso outside spectrum mode with args: $args"
			fi
			grep -Fq "option '--kso' is only for the spectrum mode" "$output" ||
				fail "$binary did not report --kso spectrum-mode validation for args: $args"
		done
	done
	pass "--kso requires spectrum mode regardless of argument order"
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
		case "$1" in
			-O)
				shift
				echo "arg=<$1>"
				if [ -f "$1" ]; then
					echo "object-bytes=$(wc -c < "$1" | tr -d ' ')"
				else
					echo "object-missing"
				fi
				;;
			-P)
				shift
				echo "arg=<$1>"
				if [ -f "$1" ]; then
					echo "pair1-bytes=$(wc -c < "$1" | tr -d ' ')"
				else
					echo "pair1-missing"
				fi
				shift
				echo "arg=<$1>"
				if [ -f "$1" ]; then
					echo "pair2-bytes=$(wc -c < "$1" | tr -d ' ')"
				else
					echo "pair2-missing"
				fi
				;;
		esac
		shift
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
		"$REPO_DIR/scripts/classify_metagenome.sh" -O "$input" -R "$result" -m 2 -n 4

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
		"$REPO_DIR/scripts/classify_metagenome.sh" -O "$input" -R "$result" --gzipped

	grep -Fq "object-bytes=16" "$capture" || fail "gzipped input was not decompressed before classification"
	if find "$tmp" -maxdepth 1 -type d -name 'CLARKGZP.*' | grep -q .; then
		fail "temporary gzipped-input directory was not cleaned up"
	fi
	pass "classify wrapper decompresses gzipped input"
}

test_classify_wrapper_paired_light_variant() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-paired-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	input1="$tmp/mate 1.fastq"
	input2="$tmp/mate 2.fastq"
	result="$tmp/result.csv"
	fake_exe="$tmp/fake-exe"
	capture="$tmp/capture.txt"

	mkdir -p "$dbdir"
	printf 'target.fa 12345\n' > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir/" > "$settings"
	printf '@r1\nACGT\n+\n!!!!\n' > "$input1"
	printf '@r2\nTGCA\n+\n!!!!\n' > "$input2"
	create_fake_exe "$fake_exe"

	CLARK_SETTINGS_FILE="$settings" \
	CLARK_EXE_DIR="$fake_exe" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/scripts/classify_metagenome.sh" -P "$input1" "$input2" -R "$result" --light

	grep -Fq "binary=CLARK-l" "$capture" || fail "--light did not select CLARK-l"
	grep -Fq "arg=<$input1>" "$capture" || fail "paired input 1 path was not preserved"
	grep -Fq "arg=<$input2>" "$capture" || fail "paired input 2 path was not preserved"
	grep -Fq "pair1-bytes=16" "$capture" || fail "paired input 1 was not passed to classifier"
	grep -Fq "pair2-bytes=16" "$capture" || fail "paired input 2 was not passed to classifier"
	pass "classify wrapper preserves paired-end paths and selects CLARK-l"
}

test_classify_wrapper_rejects_conflicting_variants() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-conflict-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	input="$tmp/input.fastq"
	result="$tmp/result.csv"
	fake_exe="$tmp/fake-exe"
	capture="$tmp/capture.txt"

	mkdir -p "$dbdir"
	printf 'target.fa 12345\n' > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir/" > "$settings"
	printf '@r1\nACGT\n+\n!!!!\n' > "$input"
	create_fake_exe "$fake_exe"

	if CLARK_SETTINGS_FILE="$settings" \
		CLARK_EXE_DIR="$fake_exe" \
		CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/scripts/classify_metagenome.sh" -O "$input" -R "$result" --light --spaced >/dev/null 2>&1; then
		fail "classify wrapper accepted conflicting --light and --spaced options"
	fi
	pass "classify wrapper rejects conflicting variants"
}

test_scripts_directory_entrypoint() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-scripts-layout-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	input="$tmp/input.fastq"
	result="$tmp/result.csv"
	fake_exe="$tmp/fake-exe"
	capture="$tmp/capture.txt"

	mkdir -p "$dbdir"
	printf 'target.fa 12345\n' > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir/" > "$settings"
	printf '@r1\nACGT\n+\n!!!!\n' > "$input"
	create_fake_exe "$fake_exe"

	CLARK_SETTINGS_FILE="$settings" \
	CLARK_EXE_DIR="$fake_exe" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/scripts/classify_metagenome.sh" -O "$input" -R "$result"

	grep -Fq "binary=CLARK" "$capture" || fail "scripts/ classify entrypoint did not run CLARK"
	grep -Fq "arg=<$input>" "$capture" || fail "scripts/ classify entrypoint did not preserve input path"
	pass "scripts directory entrypoints resolve the repository root"
}

test_set_targets_records_absolute_db_paths() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-set-targets-path-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake home"
	dbdir_input="$tmp/db dir"

	mkdir -p "$fake_home/scripts" "$fake_home/exe" "$dbdir_input/Custom" "$dbdir_input/taxonomy"
	dbdir="$(cd -P "$dbdir_input" >/dev/null 2>&1 && pwd)"
	ref="$dbdir/Custom/ref A.fa"
	ln -s "$REPO_DIR/scripts/make_metadata.sh" "$fake_home/scripts/make_metadata.sh"
	ln -s "$REPO_DIR/scripts/download_taxondata.sh" "$fake_home/scripts/download_taxondata.sh"
	ln -s "$REPO_DIR/scripts/download_RefSeqDB.sh" "$fake_home/scripts/download_RefSeqDB.sh"
	for binary in getTargetsDef getfilesToTaxNodes getAccssnTaxID; do
		ln -s "$REPO_DIR/exe/$binary" "$fake_home/exe/$binary"
	done

	printf '>NC_000001.1 synthetic custom reference\nACGT\n' > "$ref"
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$ref" "111" "111" "222" "333" "444" "555" "666" \
		> "$dbdir/.custom.fileToTaxIDs"
	printf '%s\t%s\t%s\n' "$ref" "NC_000001.1" "111" > "$dbdir/.custom.fileToAccssnTaxID"
	touch "$dbdir/.taxondata"

	(
		cd "$tmp"
		CLARK_HOME="$fake_home" "$REPO_DIR/scripts/set_targets.sh" "db dir" custom --species >/dev/null
	)

	grep -Fq -- "-T $dbdir/targets.txt" "$fake_home/.settings" || fail "set_targets did not store an absolute targets path"
	grep -Fq -- "-D $dbdir/custom_0/" "$fake_home/.settings" || fail "set_targets did not store an absolute database path"
	[ "$(cat "$fake_home/.DBDirectory")" = "$dbdir" ] || fail ".DBDirectory did not record the absolute database directory"
	[ "$(cat "$fake_home/.dbAddress")" = "$dbdir/custom_0" ] || fail ".dbAddress did not record the absolute k-mer database directory"
	pass "set_targets records absolute database paths from another working directory"
}

test_set_targets_passes_refseq_download_options() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-set-targets-download-options-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir_input="$tmp/db"
	capture="$tmp/capture.txt"
	mkdir -p "$fake_home/scripts" "$fake_home/exe" "$dbdir_input"
	dbdir="$(cd -P "$dbdir_input" >/dev/null 2>&1 && pwd)"

	cat > "$fake_home/scripts/make_metadata.sh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
db="$1"
dbdir="$2"
{
	printf 'threads=%s\n' "${CLARK_REFSEQ_THREADS:-}"
	printf 'resume=%s\n' "${CLARK_REFSEQ_RESUME:-}"
	printf 'category=%s\n' "${CLARK_REFSEQ_CATEGORY:-}"
	printf 'assembly=%s\n' "${CLARK_REFSEQ_ASSEMBLY_LEVEL:-}"
} >> "$CLARK_TEST_CAPTURE"
printf '%s\n' "$dbdir/ref.fa" > "$dbdir/.$db"
touch "$dbdir/.taxondata"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$dbdir/ref.fa" "111" "111" "222" "333" "444" "555" "666" > "$dbdir/.$db.fileToTaxIDs"
STUB
	cat > "$fake_home/exe/getTargetsDef" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
cat "$1" >/dev/null
printf 'ref.fa\t111\n'
STUB
	chmod +x "$fake_home/scripts/make_metadata.sh" "$fake_home/exe/getTargetsDef"

	CLARK_HOME="$fake_home" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/scripts/set_targets.sh" "$dbdir_input" bacteria --download-threads 4 --resume-downloads --refseq-category representative --assembly-level "Complete Genome" --genus >/dev/null

	grep -Fq "threads=4" "$capture" || fail "set_targets did not pass download thread count"
	grep -Fq "resume=1" "$capture" || fail "set_targets did not pass resume mode"
	grep -Fq "category=representative" "$capture" || fail "set_targets did not pass RefSeq category"
	grep -Fq "assembly=Complete Genome" "$capture" || fail "set_targets did not pass assembly level"
	grep -Fq -- "-D $dbdir/bacteria_1/" "$fake_home/.settings" || fail "set_targets did not preserve taxonomy rank while parsing download options"
	pass "set_targets passes RefSeq download options to metadata preparation"
}

test_set_targets_defaults_to_parallel_resume_downloads() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-set-targets-download-defaults-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir_input="$tmp/db"
	capture="$tmp/capture.txt"
	mkdir -p "$fake_home/scripts" "$fake_home/exe" "$dbdir_input"

	cat > "$fake_home/scripts/make_metadata.sh" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
db="$1"
dbdir="$2"
{
	printf 'threads=%s\n' "${CLARK_REFSEQ_THREADS:-}"
	printf 'resume=%s\n' "${CLARK_REFSEQ_RESUME:-}"
} >> "$CLARK_TEST_CAPTURE"
printf '%s\n' "$dbdir/ref.fa" > "$dbdir/.$db"
touch "$dbdir/.taxondata"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$dbdir/ref.fa" "111" "111" "222" "333" "444" "555" "666" > "$dbdir/.$db.fileToTaxIDs"
STUB
	cat > "$fake_home/exe/getTargetsDef" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
cat "$1" >/dev/null
printf 'ref.fa\t111\n'
STUB
	chmod +x "$fake_home/scripts/make_metadata.sh" "$fake_home/exe/getTargetsDef"

	CLARK_HOME="$fake_home" \
	CLARK_TEST_CAPTURE="$capture" \
		"$REPO_DIR/scripts/set_targets.sh" "$dbdir_input" bacteria >/dev/null

	grep -Fq "threads=8" "$capture" || fail "set_targets does not default RefSeq download threads to 8"
	grep -Fq "resume=1" "$capture" || fail "set_targets does not default RefSeq resume mode to on"
	pass "set_targets defaults to parallel resumable RefSeq downloads"
}

test_update_taxonomy_requires_db_directory() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-update-taxonomy-state-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	err="$tmp/update-taxonomy.err"

	mkdir -p "$fake_home"

	if CLARK_HOME="$fake_home" "$REPO_DIR/scripts/updateTaxonomy.sh" >/dev/null 2>"$err"; then
		fail "updateTaxonomy accepted a missing .DBDirectory state file"
	fi

	grep -Fq "scripts/set_targets.sh" "$err" || fail "updateTaxonomy did not explain how to configure the database directory"
	pass "updateTaxonomy reports missing database configuration"
}

create_refseq_gzip_fixture() {
	local root="$1"
	local accession="$2"
	local header="${3:-$accession}"
	mkdir -p "$root/$accession"
	printf '>%s synthetic reference\nACGTACGT\n' "$header" | gzip > "$root/$accession/${accession}_genomic.fna.gz"
}

create_refseq_static_gzip_fixture() {
	local root="$1"
	local path="$2"
	mkdir -p "$root/$(dirname "$path")"
	printf '>%s synthetic reference\nACGTACGT\n' "$(basename "$path" .gz)" | gzip > "$root/$path"
}

start_refseq_fixture_server() {
	local root="$1"
	local mode="$2"
	local state="$3"
	local port_file="$4"
	local server_script="$5"
	cat > "$server_script" <<'PY'
#!/usr/bin/env python3
import http.server
import os
import re
import socketserver
import sys

root, mode, state, port_file = sys.argv[1:5]

class Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        return

    def do_GET(self):
        name = os.path.basename(self.path)
        if mode == "transient" and re.match(r"GCF_10000000[3-7]\.1_.*_genomic\.fna\.gz", name):
            marker = os.path.join(state, name + ".failed-once")
            try:
                os.mkdir(marker)
                self.send_error(503, "temporary fixture failure")
                return
            except FileExistsError:
                pass
        if mode == "deferred":
            if name == "GCF_300000001.1_Deferred1_genomic.fna.gz" and not os.path.isdir(os.path.join(state, "first-pass-finished")):
                self.send_error(503, "temporary fixture outage")
                return
            if name == "GCF_300000005.1_Deferred5_genomic.fna.gz":
                os.makedirs(os.path.join(state, "first-pass-finished"), exist_ok=True)
        if mode == "permanent" and name == "GCF_400000002.1_Failure2_genomic.fna.gz":
            self.send_error(503, "permanent fixture failure")
            return

        accession = re.sub(r"_genomic\.fna\.gz$", "", name)
        path = os.path.join(root, name)
        nested_path = os.path.join(root, accession, name)
        if not os.path.isfile(path) and os.path.isfile(nested_path):
            path = nested_path
        if not os.path.isfile(path):
            self.send_error(404, "missing fixture")
            return
        self.send_response(200)
        self.send_header("Content-Length", str(os.path.getsize(path)))
        self.end_headers()
        with open(path, "rb") as handle:
            self.wfile.write(handle.read())

with socketserver.ThreadingTCPServer(("127.0.0.1", 0), Handler) as httpd:
    with open(port_file, "w") as handle:
        handle.write(str(httpd.server_address[1]))
    httpd.serve_forever()
PY
	chmod +x "$server_script"
	python3 "$server_script" "$root" "$mode" "$state" "$port_file" > "$state/server.stdout" 2> "$state/server.stderr" &
	local pid="$!"
	for _ in $(seq 1 50); do
		[ -s "$port_file" ] && break
		sleep 0.1
	done
	[ -s "$port_file" ] || fail "fixture HTTP server did not start"
	printf '%s\n' "$pid"
}

create_refseq_success_wget_stub() {
	local bindir="$1"
	mkdir -p "$bindir"
	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
quiet=0
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet)
			quiet=1
			shift
			;;
		-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
[ -n "$out" ] || {
	echo "wget stub missing -O output for $url" >&2
	exit 1
}
if [ -n "${CLARK_TEST_WGET_LOG:-}" ]; then
	printf 'quiet=%s url=%s\n' "$quiet" "$url" >> "$CLARK_TEST_WGET_LOG"
fi
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			printf 'GCF_900000001.1\tna\tna\tna\trepresentative genome\t10239\t10239\tMatrix virus\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tMatrixVirus\tCLARK\tna\tna\t%s/GCF_900000001.1_MatrixVirus/\n' "$CLARK_TEST_BASE_URL"
		} > "$out"
		;;
	*/protozoa/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			printf 'GCF_900000002.1\tna\tna\tna\trepresentative genome\t5759\t5759\tMatrix protozoan\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-02\tMatrixProtozoa\tCLARK\tna\tna\t%s/GCF_900000002.1_MatrixProtozoa/\n' "$CLARK_TEST_BASE_URL"
		} > "$out"
		;;
	*.fna.gz)
		file="${url##*/}"
		printf '>%s synthetic reference\nACGTACGT\n' "${file%.gz}" | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"
}

test_refseq_downloader_dry_run_manifest() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-dry-run-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	out="$tmp/dry-run.out"

	"$REPO_DIR/scripts/download_RefSeqDB.sh" --dry-run "$dbdir" plasmid > "$out"

	grep -Fq "Dry run" "$out" || fail "RefSeq downloader dry-run did not report dry-run mode"
	require_file "$dbdir/.plasmid.download_manifest.tsv"
	require_file "$dbdir/.plasmid.provenance.tsv"
	grep -Fq "plasmid.1.1.genomic.fna.gz" "$dbdir/.plasmid.download_manifest.tsv" || fail "RefSeq downloader dry-run did not plan plasmid downloads"
	grep -Fq "plasmid.1.1" "$dbdir/.plasmid.provenance.tsv" || fail "RefSeq downloader dry-run did not record static accession provenance"
	[ ! -d "$dbdir/Plasmid" ] || fail "RefSeq downloader dry-run created a sequence directory"
	pass "RefSeq downloader dry-run writes manifest and provenance"
}

test_refseq_downloader_normalizes_ncbi_trailing_slash_paths() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-ncbi-url-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	log="$tmp/wget.log"
	mkdir -p "$bindir"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
while [ "${1:-}" = "-q" ] || [ "${1:-}" = "--quiet" ] || [ "${1:-}" = "-c" ]; do
	shift
done
if [ "$1" = "-O" ]; then
	out="$2"
	url="$3"
	printf '%s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	case "$url" in
		*/bacteria/assembly_summary.txt)
			{
				printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
				printf 'GCF_055383595.1\tna\tna\tna\tna\t111\t111\tMock bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tASM5538359v1\tNCBI\tna\tna\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/055/383/595/GCF_055383595.1_ASM5538359v1/\n'
				printf 'GCF_900128725.1\tna\tna\tna\tna\t222\t222\tMock bacterium 2\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tBCifornacula_v1.0\tNCBI\tna\tna\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/128/725/GCF_900128725.1_BCifornacula_v1.0/\n'
			} > "$out"
			;;
		*/archaea/assembly_summary.txt)
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n' > "$out"
			;;
		*)
			echo "unexpected wget -O URL: $url" >&2
			exit 1
			;;
	esac
else
	echo "unexpected sequence download during dry-run: $*" >&2
	exit 1
fi
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_WGET_LOG="$log" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --dry-run "$dbdir" bacteria > "$tmp/downloader.out"

	grep -Fq "Selected 2 bacteria RefSeq genome(s)" "$tmp/downloader.out" || fail "RefSeq downloader did not select the NCBI-shaped fixture rows"
	grep -Fq "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/055/383/595/GCF_055383595.1_ASM5538359v1/GCF_055383595.1_ASM5538359v1_genomic.fna.gz" "$dbdir/.bacteria.download_manifest.tsv" ||
		fail "RefSeq downloader did not normalize the real NCBI trailing-slash ftp_path shape"
	grep -Fq "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/128/725/GCF_900128725.1_BCifornacula_v1.0/GCF_900128725.1_BCifornacula_v1.0_genomic.fna.gz" "$dbdir/.bacteria.download_manifest.tsv" ||
		fail "RefSeq downloader did not normalize NCBI assembly_summary FTP paths to HTTPS"
	if grep -Fq "ftp://ftp.ncbi.nlm.nih.gov" "$dbdir/.bacteria.download_manifest.tsv"; then
		fail "RefSeq downloader left an NCBI FTP URL in the generated manifest"
	fi
	if grep -Fq "//_genomic.fna.gz" "$dbdir/.bacteria.download_manifest.tsv"; then
		fail "RefSeq downloader emitted the broken double-slash empty-basename URL"
	fi
	pass "RefSeq downloader normalizes NCBI trailing-slash assembly_summary paths"
}

test_refseq_downloader_skips_non_url_ftp_path_rows() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-non-url-ftp-path-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	log="$tmp/wget.log"
	mkdir -p "$bindir"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
while [ "${1:-}" = "-q" ] || [ "${1:-}" = "--quiet" ] || [ "${1:-}" = "-c" ]; do
	shift
done
if [ "$1" = "-O" ]; then
	out="$2"
	url="$3"
	printf '%s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	case "$url" in
		*/bacteria/assembly_summary.txt)
			{
				printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
				printf 'GCF_055383595.1\tna\tna\tna\tna\t111\t111\tMock bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tASM5538359v1\tNCBI\tna\tna\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/055/383/595/GCF_055383595.1_ASM5538359v1/\n'
				printf 'GCF_999999999.1\tna\tna\tna\tna\t333\t333\tMalformed-row bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tASM999v1\tNCBI\tna\tna\tidentical\n'
			} > "$out"
			;;
		*/archaea/assembly_summary.txt)
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n' > "$out"
			;;
		*)
			echo "unexpected wget -O URL: $url" >&2
			exit 1
			;;
	esac
else
	echo "unexpected sequence download during dry-run: $*" >&2
	exit 1
fi
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_WGET_LOG="$log" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --dry-run "$dbdir" bacteria > "$tmp/downloader.out"

	grep -Fq "Selected 1 bacteria RefSeq genome(s)" "$tmp/downloader.out" ||
		fail "RefSeq downloader did not skip the row with a non-URL ftp_path"
	grep -Fq "GCF_055383595.1_ASM5538359v1_genomic.fna.gz" "$dbdir/.bacteria.download_manifest.tsv" ||
		fail "RefSeq downloader lost the valid row alongside the malformed one"
	if grep -Fq "identical" "$dbdir/.bacteria.download_manifest.tsv" "$dbdir/.bacteria.provenance.tsv"; then
		fail "RefSeq downloader generated a URL from a non-URL ftp_path value"
	fi
	pass "RefSeq downloader skips assembly_summary rows whose ftp_path is not a real URL"
}

test_refseq_downloader_rejects_malformed_urls_before_download() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-malformed-url-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	static_urls="$tmp/static-urls.tsv"
	printf 'plasmid\trefseq_static\thttps://example.org/refseq/GCF_055383595.1_ASM5538359v1//_genomic.fna.gz\n' > "$static_urls"

	if CLARK_REFSEQ_STATIC_URLS="$static_urls" "$REPO_DIR/scripts/download_RefSeqDB.sh" --dry-run "$dbdir" plasmid > "$tmp/stdout" 2> "$tmp/stderr"; then
		fail "RefSeq downloader accepted a malformed empty-basename URL"
	fi
	grep -Fq "generated malformed RefSeq download URL" "$tmp/stderr" ||
		fail "RefSeq downloader did not explain malformed URL rejection"
	pass "RefSeq downloader rejects malformed generated URLs before download"
}

test_refseq_downloader_tty_progress_rewrites_line() {
	python3 - "$REPO_DIR" <<'PY'
import importlib.util
import io
from pathlib import Path
import sys

repo = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("refseq_downloader", repo / "scripts" / "refseq_downloader.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class TtyBuffer(io.StringIO):
    def isatty(self):
        return True

buffer = TtyBuffer()
original_stdout = sys.stdout
try:
    sys.stdout = buffer
    module.print_progress("RefSeq download progress", 0, 100)
    module.print_progress("RefSeq download progress", 50, 100)
    module.print_progress("RefSeq download progress", 100, 100)
finally:
    sys.stdout = original_stdout

expected = (
    "\rRefSeq download progress: 0/100 files complete (0%)."
    "\rRefSeq download progress: 50/100 files complete (50%)."
    "\rRefSeq download progress: 100/100 files complete (100%).\n"
)
if buffer.getvalue() != expected:
    raise SystemExit("unexpected TTY progress output: %r" % buffer.getvalue())
PY
	pass "RefSeq downloader rewrites interactive terminal progress in place"
}

test_refseq_downloader_quiet_aggregate_progress() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-quiet-progress-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	fixtures="$tmp/refseq"
	log="$tmp/wget.log"
	mkdir -p "$bindir" "$fixtures"
	for i in 1 2 3 4 5 6 7 8 9 10; do
		create_refseq_gzip_fixture "$fixtures" "$(printf 'GCF_000000%03d.1_Mock%d' "$i" "$i")" "mock-virus-$i"
	done

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
quiet=0
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet)
			quiet=1
			shift
			;;
		-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
printf 'quiet=%s url=%s\n' "$quiet" "$url" >> "$CLARK_TEST_WGET_LOG"
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			for i in 1 2 3 4 5 6 7 8 9 10; do
				printf 'GCF_000000%03d.1\tna\tna\tna\trepresentative genome\t10239\t10239\tMock virus %d\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tMock%d\tCLARK\tna\tna\t%s/GCF_000000%03d.1_Mock%d/\n' "$i" "$i" "$i" "$CLARK_TEST_BASE_URL" "$i" "$i"
			done
		} > "$out"
		;;
	*_genomic.fna.gz)
		printf '>mock-virus\nACGT\n' | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_WGET_LOG="$log" \
	CLARK_TEST_BASE_URL="file://$fixtures" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" "$dbdir" viruses > "$tmp/stdout" 2> "$tmp/stderr"

	if grep -Fq "quiet=0" "$log"; then
		fail "RefSeq downloader invoked wget without quiet mode"
	fi
	grep -Fq "RefSeq download progress:" "$tmp/stdout" || fail "RefSeq downloader did not report aggregate progress"
	grep -Fq "10/10" "$tmp/stdout" || fail "RefSeq downloader did not report final aggregate completion"
	if grep -Eq "Saving to:|HTTP request sent|Resolving |Connecting to " "$tmp/stdout" "$tmp/stderr"; then
		fail "RefSeq downloader leaked per-file network chatter"
	fi
	pass "RefSeq downloader uses quiet network commands and aggregate progress"
}

test_refseq_downloader_retries_parallel_transient_failures() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-parallel-retry-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	state="$tmp/state"
	fixtures="$tmp/refseq"
	mkdir -p "$bindir" "$state" "$fixtures"
	for i in 1 2 3 4 5 6 7 8 9 10; do
		create_refseq_gzip_fixture "$fixtures" "$(printf 'GCF_100000%03d.1_Retry%d' "$i" "$i")" "retry-virus-$i"
	done
	server_pid="$(start_refseq_fixture_server "$fixtures" transient "$state" "$tmp/port" "$tmp/refseq_server.py")"
	trap 'kill "$server_pid" 2>/dev/null || true; rm -rf "$tmp"' RETURN
	base_url="http://127.0.0.1:$(cat "$tmp/port")/refseq"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet|-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			for i in 1 2 3 4 5 6 7 8 9 10; do
				printf 'GCF_100000%03d.1\tna\tna\tna\trepresentative genome\t10239\t10239\tMock virus %d\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tRetry%d\tCLARK\tna\tna\t%s/GCF_100000%03d.1_Retry%d/\n' "$i" "$i" "$i" "$CLARK_TEST_BASE_URL" "$i" "$i"
			done
		} > "$out"
		;;
	*_genomic.fna.gz)
		base="${url##*/}"
		case "$base" in
			GCF_10000000[3-7].1_*)
				marker="$CLARK_TEST_STATE/$base.failed-once"
				if mkdir "$marker" 2>/dev/null; then
					printf 'partial transfer for %s\n' "$base" > "$out"
					exit 1
				fi
				;;
		esac
		printf '>retry-virus\nACGT\n' | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_STATE="$state" \
	CLARK_TEST_BASE_URL="$base_url" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
	CLARK_REFSEQ_RETRY_DELAY=0 \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 8 "$dbdir" viruses > "$tmp/stdout" 2> "$tmp/stderr"

	count=$(find "$dbdir/Viruses" -type f -name '*.fna' | wc -l | tr -d ' ')
	[ "$count" = "10" ] || fail "RefSeq downloader did not recover all transient parallel download failures"
	if find "$dbdir/Viruses" -type f -name '*.part' -print -quit | grep -q .; then
		fail "RefSeq downloader left .part files after successful retries"
	fi
	grep -Fq "GCF_100000003.1_Retry3_genomic.fna" "$dbdir/.viruses" ||
		fail "RefSeq downloader did not list a retried and decompressed FASTA"
	pass "RefSeq downloader retries transient parallel failures and finalizes archives"
}

test_refseq_downloader_defers_retry_until_after_first_pass() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-deferred-retry-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	state="$tmp/state"
	fixtures="$tmp/refseq"
	mkdir -p "$bindir" "$state" "$fixtures"
	for i in 1 2 3 4 5; do
		create_refseq_gzip_fixture "$fixtures" "$(printf 'GCF_300000%03d.1_Deferred%d' "$i" "$i")" "deferred-virus-$i"
	done
	server_pid="$(start_refseq_fixture_server "$fixtures" deferred "$state" "$tmp/port" "$tmp/refseq_server.py")"
	trap 'kill "$server_pid" 2>/dev/null || true; rm -rf "$tmp"' RETURN
	base_url="http://127.0.0.1:$(cat "$tmp/port")/refseq"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet|-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			for i in 1 2 3 4 5; do
				printf 'GCF_300000%03d.1\tna\tna\tna\trepresentative genome\t10239\t10239\tDeferred virus %d\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tDeferred%d\tCLARK\tna\tna\t%s/GCF_300000%03d.1_Deferred%d/\n' "$i" "$i" "$i" "$CLARK_TEST_BASE_URL" "$i" "$i"
			done
		} > "$out"
		;;
	*_genomic.fna.gz)
		base="${url##*/}"
		case "$base" in
			GCF_300000001.1_Deferred1_genomic.fna.gz)
				if [ ! -d "$CLARK_TEST_STATE/first-pass-finished" ]; then
					printf 'temporary NCBI outage for %s\n' "$base" > "$out"
					exit 1
				fi
				;;
			GCF_300000005.1_Deferred5_genomic.fna.gz)
				mkdir "$CLARK_TEST_STATE/first-pass-finished" 2>/dev/null || true
				;;
		esac
		printf '>deferred-virus\nACGT\n' | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_STATE="$state" \
	CLARK_TEST_BASE_URL="$base_url" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
	CLARK_REFSEQ_RETRY_DELAY=0 \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 1 "$dbdir" viruses > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq "Retrying 1 deferred RefSeq download(s)" "$tmp/stdout" ||
		fail "RefSeq downloader did not defer failed URLs for a later retry pass"
	count=$(find "$dbdir/Viruses" -type f -name '*.fna' | wc -l | tr -d ' ')
	[ "$count" = "5" ] || fail "RefSeq downloader did not complete all files after deferred retry"
	grep -Fq "deferred" "$dbdir/.viruses.download_manifest.tsv" ||
		fail "RefSeq downloader did not record the deferred first-pass status"
	grep -Fq "GCF_300000001.1_Deferred1_genomic.fna" "$dbdir/.viruses" ||
		fail "RefSeq downloader did not list the deferred and later downloaded FASTA"
	pass "RefSeq downloader retries deferred URLs after the first pass"
}

test_refseq_downloader_fails_when_deferred_downloads_remain() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-permanent-failure-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	state="$tmp/state"
	fixtures="$tmp/refseq"
	mkdir -p "$bindir" "$state" "$fixtures"
	for i in 1 2 3 4; do
		create_refseq_gzip_fixture "$fixtures" "$(printf 'GCF_400000%03d.1_Failure%d' "$i" "$i")" "failure-virus-$i"
	done
	server_pid="$(start_refseq_fixture_server "$fixtures" permanent "$state" "$tmp/port" "$tmp/refseq_server.py")"
	trap 'kill "$server_pid" 2>/dev/null || true; rm -rf "$tmp"' RETURN
	base_url="http://127.0.0.1:$(cat "$tmp/port")/refseq"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet|-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			for i in 1 2 3 4; do
				printf 'GCF_400000%03d.1\tna\tna\tna\trepresentative genome\t10239\t10239\tFailure virus %d\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tFailure%d\tCLARK\tna\tna\t%s/GCF_400000%03d.1_Failure%d/\n' "$i" "$i" "$i" "$CLARK_TEST_BASE_URL" "$i" "$i"
			done
		} > "$out"
		;;
	*GCF_400000002.1_Failure2_genomic.fna.gz)
		printf 'permanent transfer failure\n' > "$out"
		exit 1
		;;
	*_genomic.fna.gz)
		printf '>failure-virus\nACGT\n' | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"

	if PATH="$bindir:$PATH" \
		CLARK_TEST_BASE_URL="$base_url" \
		CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		CLARK_REFSEQ_RETRY_DELAY=0 \
			"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 4 "$dbdir" viruses > "$tmp/stdout" 2> "$tmp/stderr"; then
		fail "RefSeq downloader reported success despite a permanent deferred failure"
	fi

	grep -Fq "Failed to download 1 RefSeq genome file(s) after deferred retries" "$tmp/stderr" ||
		fail "RefSeq downloader did not summarize permanent deferred failures"
	[ ! -e "$dbdir/.viruses" ] || fail "RefSeq downloader wrote a success marker after an incomplete download"
	grep -Fq "failed" "$dbdir/.viruses.download_manifest.tsv" ||
		fail "RefSeq downloader did not record failed final status in the manifest"
	require_file "$dbdir/.viruses.failed_downloads.tsv"
	pass "RefSeq downloader fails loudly when deferred downloads remain incomplete"
}

test_refseq_downloader_reports_early_parallel_progress() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-progress-refresh-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	fixtures="$tmp/refseq"
	mkdir -p "$bindir" "$fixtures"
	for i in $(seq 1 200); do
		create_refseq_gzip_fixture "$fixtures" "$(printf 'GCF_200000%03d.1_Progress%d' "$i" "$i")" "progress-virus-$i"
	done

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
out=""
url=""
while [ "$#" -gt 0 ]; do
	case "$1" in
		-q|--quiet|-c)
			shift
			;;
		-O)
			out="$2"
			shift 2
			;;
		*)
			url="$1"
			shift
			;;
	esac
done
case "$url" in
	*/viral/assembly_summary.txt)
		{
			printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			for i in $(seq 1 200); do
				printf 'GCF_200000%03d.1\tna\tna\tna\trepresentative genome\t10239\t10239\tMock virus %d\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2026-06-01\tProgress%d\tCLARK\tna\tna\t%s/GCF_200000%03d.1_Progress%d/\n' "$i" "$i" "$i" "$CLARK_TEST_BASE_URL" "$i" "$i"
			done
		} > "$out"
		;;
	*_genomic.fna.gz)
		printf '>progress-virus\nACGT\n' | gzip > "$out"
		;;
	*)
		echo "unexpected URL: $url" >&2
		exit 1
		;;
esac
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_BASE_URL="file://$fixtures" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 8 "$dbdir" viruses > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq "RefSeq download progress:" "$tmp/stdout" ||
		fail "RefSeq downloader did not report progress during parallel downloads"
	grep -Fq "RefSeq download progress: 200/200 files complete (100%)." "$tmp/stdout" ||
		fail "RefSeq downloader did not report final aggregate completion"
	pass "RefSeq downloader refreshes aggregate progress during parallel downloads"
}

test_refseq_downloader_success_matrix() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-success-matrix-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	bindir="$tmp/bin"
	log="$tmp/wget.log"
	fixtures="$tmp/refseq"
	static_urls="$tmp/static-urls.tsv"
	mkdir -p "$fixtures"
	create_refseq_gzip_fixture "$fixtures" "GCF_900000001.1_MatrixVirus" "matrix-virus"
	create_refseq_gzip_fixture "$fixtures" "GCF_900000002.1_MatrixProtozoa" "matrix-protozoa"
	create_refseq_static_gzip_fixture "$fixtures" "static/plasmid.1.1.genomic.fna.gz"
	create_refseq_static_gzip_fixture "$fixtures" "static/plastid.1.1.genomic.fna.gz"
	create_refseq_success_wget_stub "$bindir"
	{
		printf 'plasmid\trefseq_static\tfile://%s/static/plasmid.1.1.genomic.fna.gz\n' "$fixtures"
		printf 'plastid\trefseq_static\tfile://%s/static/plastid.1.1.genomic.fna.gz\n' "$fixtures"
	} > "$static_urls"

	for db in viruses plasmid plastid protozoa; do
		dbdir="$tmp/db-$db"
		PATH="$bindir:$PATH" \
		CLARK_TEST_WGET_LOG="$log" \
		CLARK_TEST_BASE_URL="file://$fixtures" \
		CLARK_REFSEQ_STATIC_URLS="$static_urls" \
			"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 4 "$dbdir" "$db" > "$tmp/$db.out" 2> "$tmp/$db.err"

		require_file "$dbdir/.$db"
		require_file "$dbdir/.$db.download_manifest.tsv"
		require_file "$dbdir/.$db.provenance.tsv"
		grep -Fq "RefSeq download progress:" "$tmp/$db.out" || fail "RefSeq downloader did not report progress for $db"
		grep -Fq "downloaded" "$dbdir/.$db.download_manifest.tsv" || fail "RefSeq downloader did not record completed downloads for $db"
		if find "$dbdir" -type f -name '*.part' -print -quit | grep -q .; then
			fail "RefSeq downloader left partial files for $db"
		fi
		case "$db" in
			viruses)
				data_dir="Viruses"
				grep -Fq "GCF_900000001.1_MatrixVirus_genomic.fna" "$dbdir/.viruses" || fail "RefSeq downloader did not list decompressed virus FASTA"
				grep -Fq "2026-06-01" "$dbdir/.viruses.provenance.tsv" || fail "RefSeq downloader did not preserve virus source date"
				;;
			protozoa)
				data_dir="Protozoa"
				grep -Fq "GCF_900000002.1_MatrixProtozoa_genomic.fna" "$dbdir/.protozoa" || fail "RefSeq downloader did not list decompressed protozoa FASTA"
				grep -Fq "2026-06-02" "$dbdir/.protozoa.provenance.tsv" || fail "RefSeq downloader did not preserve protozoa source date"
				;;
			plasmid)
				data_dir="Plasmid"
				[ "$(find "$dbdir" -type f -name '*.fa' | wc -l | tr -d ' ')" -gt 0 ] || fail "RefSeq downloader did not split $db FASTA records"
				grep -Fq "$db.1.1.genomic" "$dbdir/.$db.provenance.tsv" || fail "RefSeq downloader did not record static $db provenance"
				;;
			plastid)
				data_dir="Plastid"
				[ "$(find "$dbdir" -type f -name '*.fa' | wc -l | tr -d ' ')" -gt 0 ] || fail "RefSeq downloader did not split $db FASTA records"
				grep -Fq "$db.1.1.genomic" "$dbdir/.$db.provenance.tsv" || fail "RefSeq downloader did not record static $db provenance"
				;;
		esac
		[ ! -e "$dbdir/$data_dir/download.sh" ] || fail "RefSeq downloader generated a legacy download.sh script for $db"
	done
	if grep -Fq "quiet=0" "$log"; then
		fail "RefSeq downloader success matrix invoked wget without quiet mode"
	fi
	pass "RefSeq downloader completes viruses, plasmid, plastid, and protozoa downloads"
}

test_refseq_downloader_filters_and_resumes() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-filter-resume-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	fixtures="$tmp/refseq"
	log="$tmp/wget.log"
	mkdir -p "$bindir" "$dbdir/Bacteria" "$fixtures"
	printf '>existing-bacterium\nACGT\n' | gzip > "$dbdir/Bacteria/GCF_111111111.1_Existing_genomic.fna.gz"
	create_refseq_gzip_fixture "$fixtures" "GCF_222222222.1_New" "new-bacterium"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
while [ "${1:-}" = "-q" ] || [ "${1:-}" = "--quiet" ] || [ "${1:-}" = "-c" ]; do
	shift
done
if [ "$1" = "-O" ]; then
	out="$2"
	url="$3"
	printf 'summary %s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	case "$url" in
		*/bacteria/assembly_summary.txt)
			{
				printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
				printf 'GCF_111111111.1\tna\tna\tna\trepresentative genome\t111\t111\tExisting bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2025-01-01\tExisting\tCLARK\tna\tna\t%s/GCF_111111111.1_Existing/\n' "$CLARK_TEST_BASE_URL"
				printf 'GCF_222222222.1\tna\tna\tna\trepresentative genome\t222\t222\tNew bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2025-01-02\tNew\tCLARK\tna\tna\t%s/GCF_222222222.1_New/\n' "$CLARK_TEST_BASE_URL"
				printf 'GCF_333333333.1\tna\tna\tna\tna\t333\t333\tUnselected bacterium\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2025-01-03\tUnselected\tCLARK\tna\tna\t%s/GCF_333333333.1_Unselected\n' "$CLARK_TEST_BASE_URL"
				printf 'GCF_444444444.1\tna\tna\tna\trepresentative genome\t444\t444\tDraft bacterium\tna\tna\tlatest\tScaffold\tMajor\tFull\t2025-01-04\tDraft\tCLARK\tna\tna\t%s/GCF_444444444.1_Draft\n' "$CLARK_TEST_BASE_URL"
			} > "$out"
			;;
		*/archaea/assembly_summary.txt)
			{
				printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
			} > "$out"
			;;
		*_genomic.fna.gz)
			printf 'download %s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
			printf '>new-bacterium\nTGCA\n' | gzip > "$out"
			;;
		*)
			echo "unexpected wget -O URL: $url" >&2
			exit 1
			;;
	esac
else
	url="${@: -1}"
	printf 'download %s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	file="${url##*/}"
	printf '>new-bacterium\nTGCA\n' | gzip > "$file"
fi
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_WGET_LOG="$log" \
	CLARK_TEST_BASE_URL="file://$fixtures" \
	CLARK_REFSEQ_STATIC_URLS="$tmp/missing-static-urls.tsv" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" --threads 2 --resume --refseq-category representative "$dbdir" bacteria > "$tmp/downloader.out"

	grep -Fq "Selected 2 bacteria RefSeq genome(s)" "$tmp/downloader.out" || fail "RefSeq downloader did not report filtered bacteria selection"
	grep -Fq "GCF_111111111.1_Existing_genomic.fna" "$dbdir/.bacteria" || fail "RefSeq downloader lost existing resumed bacteria FASTA"
	grep -Fq "GCF_222222222.1_New_genomic.fna" "$dbdir/.bacteria" || fail "RefSeq downloader did not list newly downloaded bacteria FASTA"
	grep -Fq "GCF_222222222.1" "$dbdir/.bacteria.provenance.tsv" || fail "RefSeq downloader did not keep new representative provenance"
	if grep -Fq "GCF_333333333.1" "$dbdir/.bacteria.provenance.tsv" || grep -Fq "GCF_444444444.1" "$dbdir/.bacteria.provenance.tsv"; then
		fail "RefSeq downloader provenance includes assemblies outside the requested filters"
	fi
	grep -Fq "skipped-existing" "$dbdir/.bacteria.download_manifest.tsv" || fail "RefSeq downloader did not record resumed existing archive"
	if grep -Fq "GCF_111111111.1_Existing_genomic.fna.gz" "$log"; then
		fail "RefSeq downloader re-downloaded an existing archive during resume"
	fi
	pass "RefSeq downloader filters RefSeq assemblies and resumes existing archives"
}

test_make_metadata_uses_refseq_provenance_taxids() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-provenance-taxid-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir="$tmp/db"
	mkdir -p "$fake_home/scripts" "$fake_home/exe" "$dbdir/Bacteria" "$dbdir/taxonomy"
	ln -s "$REPO_DIR/scripts/make_metadata.sh" "$fake_home/scripts/make_metadata.sh"
	ln -s "$REPO_DIR/scripts/download_taxondata.sh" "$fake_home/scripts/download_taxondata.sh"
	ln -s "$REPO_DIR/scripts/download_RefSeqDB.sh" "$fake_home/scripts/download_RefSeqDB.sh"
	for binary in getTargetsDef getfilesToTaxNodes getAccssnTaxID; do
		ln -s "$REPO_DIR/exe/$binary" "$fake_home/exe/$binary"
	done

	ref="$dbdir/Bacteria/GCF_555555555.1_FastPath_genomic.fna"
	printf '>NC_555555.1 mock bacteria\nACGT\n' > "$ref"
	printf '%s\n' "$ref" > "$dbdir/.bacteria"
	{
		printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
		printf 'bacteria\tbacteria\tGCF_555555555.1\t555\t555\t2025-02-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_555555555.1_FastPath/GCF_555555555.1_FastPath_genomic.fna.gz\n'
	} > "$dbdir/.bacteria.provenance.tsv"
	printf '555 | 2 | species |\n2 | 1 | superkingdom |\n' > "$dbdir/taxonomy/nodes.dmp"
	printf '1 | 1 |\n' > "$dbdir/taxonomy/merged.dmp"
	touch "$dbdir/.taxondata"

	CLARK_HOME="$fake_home" "$REPO_DIR/scripts/make_metadata.sh" bacteria "$dbdir" > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq "$ref	GCF_555555555.1	555" "$dbdir/.bacteria.fileToAccssnTaxID" || fail "make_metadata did not use RefSeq provenance taxid mapping"
	grep -Fq "$ref	555	555" "$dbdir/.bacteria.fileToTaxIDs" || fail "make_metadata did not build lineage from provenance taxid mapping"
	if grep -Fq "Re-building bacteria.fileToAccssnTaxID" "$tmp/stdout" "$tmp/stderr"; then
		fail "make_metadata fell back to global accession lookup despite usable RefSeq provenance"
	fi
	pass "make_metadata uses RefSeq provenance taxids without global accession-map lookup"
}

test_documentation_script_paths() {
	if grep -E "\./(buildSpacedDB|classify_metagenome|clean|download_RefSeqDB|download_taxondata|estimate_abundance|evaluate_density_confidence|evaluate_density_gamma|extractSequences|getTargetsKmers_distribution|install|makeSummaryTables|make_metadata|resetCustomDB|set_targets|updateTaxonomy)\.sh|\./scripts/" \
		"$REPO_DIR/README.md" "$REPO_DIR/README_FULL.md" "$REPO_DIR/docs/QUICKSTART.md" "$REPO_DIR/scripts/README.md" >/dev/null; then
		fail "documentation still contains root-relative script examples"
	fi
	pass "README documents scripts/ commands without root-relative script paths"
}

test_make_sample_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-make-sample-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	genome="$tmp/refA.fa"
	{
		printf '>refA synthetic\n'
		python3 -c "import random; random.seed(1); print(''.join(random.choice('ACGT') for _ in range(2000)))"
	} > "$genome"

	targets="$tmp/targets.txt"
	printf '%s\t111\n' "$genome" > "$targets"

	settings="$tmp/.settings"
	printf -- '-T %s\n-D %s/db/\n' "$targets" "$tmp" > "$settings"

	out="$tmp/sample.fa"
	CLARK_SETTINGS_FILE="$settings" "$REPO_DIR/scripts/make_sample.sh" -n 4 -l 50 -o "$out" --seed 42 >/dev/null

	[ -s "$out" ] || fail "make_sample.sh did not write an output file"
	count=$(grep -c '^>' "$out")
	[ "$count" -eq 4 ] || fail "make_sample.sh did not write the requested number of reads (got $count)"
	grep -Fq "taxid=111" "$out" || fail "make_sample.sh did not tag reads with the target taxid"

	if grep -v '^>' "$out" | grep -qvE '^[ACGT]{50}$'; then
		fail "make_sample.sh emitted a read that is not exactly the requested length"
	fi

	pass "make_sample.sh samples reads of the requested count and length from configured targets"
}

test_make_sample_requires_configured_targets() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-make-sample-missing-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	err="$tmp/stderr"
	if CLARK_SETTINGS_FILE="$tmp/missing-settings" "$REPO_DIR/scripts/make_sample.sh" -n 3 -o "$tmp/sample.fa" >/dev/null 2>"$err"; then
		fail "make_sample.sh accepted a missing settings file"
	fi
	grep -Fq "scripts/set_targets.sh" "$err" || fail "make_sample.sh did not explain how to configure targets"
	pass "make_sample.sh requires configured targets before sampling"
}

test_make_benchmark_reads_hiseq_profile() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-make-benchmark-hiseq-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	genome="$tmp/refA.fa"
	{
		printf '>refA synthetic\n'
		python3 -c "import random; random.seed(2); print(''.join(random.choice('ACGT') for _ in range(5000)))"
	} > "$genome"

	targets="$tmp/targets.txt"
	printf '%s\t222\n' "$genome" > "$targets"

	settings="$tmp/.settings"
	printf -- '-T %s\n-D %s/db/\n' "$targets" "$tmp" > "$settings"

	out="$tmp/benchmark.fa"
	truth="$tmp/benchmark.fa.truth.tsv"
	CLARK_SETTINGS_FILE="$settings" "$REPO_DIR/scripts/make_benchmark_reads.sh" -p hiseq -n 15 -o "$out" --seed 5 >/dev/null

	[ -s "$out" ] || fail "make_benchmark_reads.sh did not write an output FASTA"
	count=$(grep -c '^>' "$out")
	[ "$count" -eq 15 ] || fail "make_benchmark_reads.sh did not write the requested number of reads (got $count)"

	if grep -v '^>' "$out" | grep -qvE '^[ACGT]{92}$'; then
		fail "make_benchmark_reads.sh hiseq profile did not emit 92 bp reads"
	fi
	grep -Fq "taxid=222" "$out" || fail "make_benchmark_reads.sh did not tag reads with the source taxid"
	grep -Fq "profile=hiseq" "$out" || fail "make_benchmark_reads.sh did not tag reads with their profile"

	require_file "$truth"
	[ "$(tail -n +2 "$truth" | wc -l)" -eq 15 ] || fail "ground-truth TSV does not have one row per read"
	head -n1 "$truth" | grep -Fq "read_id" || fail "ground-truth TSV is missing its header"
	grep -Fq "222" "$truth" || fail "ground-truth TSV did not record the source taxid"

	pass "make_benchmark_reads.sh hiseq profile emits 92 bp reads with matching ground truth"
}

test_make_benchmark_reads_error_injection() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-make-benchmark-error-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	genome="$tmp/refA.fa"
	seq="$(python3 -c "import random; random.seed(3); print(''.join(random.choice('ACGT') for _ in range(2000)))")"
	printf '>refA synthetic\n%s\n' "$seq" > "$genome"

	targets="$tmp/targets.txt"
	printf '%s\t333\n' "$genome" > "$targets"

	out="$tmp/benchmark.fa"
	truth="$tmp/benchmark.fa.truth.tsv"
	"$REPO_DIR/scripts/make_benchmark_reads.sh" -p custom -n 5 -l 40 -e 1.0 -T "$targets" -o "$out" --seed 9 >/dev/null

	if python3 - "$out" "$seq" <<'PYEOF'
import sys

out_path, seq = sys.argv[1], sys.argv[2]
with open(out_path) as handle:
    lines = [line.rstrip("\n") for line in handle]

i = 0
checked = 0
while i < len(lines):
    header = lines[i]
    read = lines[i + 1]
    i += 2
    pos = header.split("pos=")[1].split("|")[0]
    start, end = (int(x) for x in pos.split("-"))
    original = seq[start:end]
    for a, b in zip(original, read):
        if a == b:
            sys.exit("error_rate=1.0 should mutate every base, but found an unchanged base")
    checked += 1

if checked != 5:
    sys.exit("expected to check 5 reads, checked %d" % checked)
PYEOF
	then :; else fail "make_benchmark_reads.sh error injection did not mutate every base at error_rate=1.0"; fi

	out0="$tmp/benchmark0.fa"
	"$REPO_DIR/scripts/make_benchmark_reads.sh" -p custom -n 5 -l 40 -e 0 -T "$targets" -o "$out0" --seed 9 >/dev/null
	if python3 - "$out0" "$seq" <<'PYEOF'
import sys

out_path, seq = sys.argv[1], sys.argv[2]
with open(out_path) as handle:
    lines = [line.rstrip("\n") for line in handle]

i = 0
while i < len(lines):
    header = lines[i]
    read = lines[i + 1]
    i += 2
    pos = header.split("pos=")[1].split("|")[0]
    start, end = (int(x) for x in pos.split("-"))
    original = seq[start:end]
    if original != read:
        sys.exit("error_rate=0 should leave reads identical to the source genome")
PYEOF
	then :; else fail "make_benchmark_reads.sh emitted a mutated base at error_rate=0"; fi

	pass "make_benchmark_reads.sh error injection respects the requested substitution rate"
}

test_eval_accuracy_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-eval-accuracy-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	truth="$tmp/truth.tsv"
	{
		printf 'read_id\ttrue_taxid\tsource_file\tprofile\n'
		printf 'read1\t100\tgenomeA\tcustom\n'
		printf 'read2\t100\tgenomeA\tcustom\n'
		printf 'read3\t200\tgenomeB\tcustom\n'
		printf 'read4\t200\tgenomeB\tcustom\n'
		printf 'read5\t300\tgenomeC\tcustom\n'
	} > "$truth"

	results="$tmp/results.csv"
	{
		printf 'Object_ID, Length, Assignment\n'
		printf 'read1,40,100\n'
		printf 'read2,40,200\n'
		printf 'read3,40,NA\n'
		printf 'read4,40,200\n'
		printf 'read5,40,100\n'
	} > "$results"

	out="$tmp/stdout"
	python3 "$REPO_DIR/scripts/eval_accuracy.py" --results "$results" --ground-truth "$truth" --per-taxid > "$out"

	grep -Fq "reads_total=5 tp=2 fp=2 fn=3 sensitivity=0.400000 precision=0.500000" "$out" ||
		fail "eval_accuracy.py did not compute the expected micro-averaged sensitivity/precision"
	grep -Fq "taxid=100 tp=1 fp=1 fn=1 sensitivity=0.500000 precision=0.500000" "$out" ||
		fail "eval_accuracy.py did not compute the expected per-taxid breakdown for taxid 100"
	grep -Fq "taxid=200 tp=1 fp=1 fn=1 sensitivity=0.500000 precision=0.500000" "$out" ||
		fail "eval_accuracy.py did not compute the expected per-taxid breakdown for taxid 200"
	grep -Fq "taxid=300 tp=0 fp=0 fn=1 sensitivity=0.000000 precision=0.000000" "$out" ||
		fail "eval_accuracy.py did not compute the expected per-taxid breakdown for taxid 300"

	pass "eval_accuracy.py computes micro-averaged and per-taxid sensitivity/precision"
}

test_rehome_db_fixes_stale_absolute_paths() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rehome-db-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	# Simulate copying a populated db directory from one machine/home
	# directory to another: the genome files live under the *new* prefix,
	# but the cached metadata (.protozoa, .fileToAccssnTaxID,
	# .fileToTaxIDs) still has the *old* machine's absolute paths baked in
	# -- exactly the failure a user hit after rsync-ing ncbi-db/ between
	# /home/h/... and /home/u/... hosts.
	new_dbdir="$tmp/new-home/ncbi-db"
	mkdir -p "$new_dbdir/Protozoa"
	genome="$new_dbdir/Protozoa/GCF_000091205.1_ASM9120v1_genomic.fna"
	{
		printf '>NC_000001.1 synthetic protozoan\n'
		python3 -c "import random; random.seed(31); print(''.join(random.choice('ACGT') for _ in range(500)))"
	} > "$genome"

	old_genome="/home/h/sykim/kmer-matching-project/ncbi-db/Protozoa/GCF_000091205.1_ASM9120v1_genomic.fna"
	printf '%s\n' "$old_genome" > "$new_dbdir/.protozoa"
	printf '%s\tGCF_000091205.1\t5722\n' "$old_genome" > "$new_dbdir/.protozoa.fileToAccssnTaxID"
	printf '%s\t5722\tspecies\n' "$old_genome" > "$new_dbdir/.protozoa.fileToTaxIDs"

	out="$tmp/stdout"
	err="$tmp/stderr"
	"$REPO_DIR/scripts/rehome_db.sh" "$new_dbdir" >"$out" 2>"$err" || fail "rehome_db.sh exited non-zero: $(cat "$err")"

	grep -Fq "$genome" "$new_dbdir/.protozoa" || fail "rehome_db.sh did not repair .protozoa"
	grep -Fq "$genome" "$new_dbdir/.protozoa.fileToAccssnTaxID" || fail "rehome_db.sh did not repair .fileToAccssnTaxID"
	grep -Fq "$genome" "$new_dbdir/.protozoa.fileToTaxIDs" || fail "rehome_db.sh did not repair .fileToTaxIDs"
	grep -Fq "GCF_000091205.1" "$new_dbdir/.protozoa.fileToAccssnTaxID" || fail "rehome_db.sh dropped fields other than the path"
	grep -Fq "5722" "$new_dbdir/.protozoa.fileToTaxIDs" || fail "rehome_db.sh dropped the taxid field"

	if grep -Fq "/home/h/" "$new_dbdir/.protozoa" "$new_dbdir/.protozoa.fileToAccssnTaxID" "$new_dbdir/.protozoa.fileToTaxIDs"; then
		fail "rehome_db.sh left a stale old-host path behind"
	fi

	pass "rehome_db.sh repairs stale absolute paths after copying a db directory to a new host"
}

test_rehome_db_dry_run_leaves_files_untouched() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rehome-db-dry-run-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	new_dbdir="$tmp/new-home/ncbi-db"
	mkdir -p "$new_dbdir/Viruses"
	genome="$new_dbdir/Viruses/GCF_2.fna"
	printf '>seq\nACGTACGT\n' > "$genome"

	old_genome="/home/h/old-project/ncbi-db/Viruses/GCF_2.fna"
	printf '%s\n' "$old_genome" > "$new_dbdir/.viruses"
	before_hash="$(sha256sum "$new_dbdir/.viruses")"

	"$REPO_DIR/scripts/rehome_db.sh" "$new_dbdir" --dry-run >/dev/null 2>&1 || fail "rehome_db.sh --dry-run exited non-zero"

	after_hash="$(sha256sum "$new_dbdir/.viruses")"
	[ "$before_hash" = "$after_hash" ] || fail "rehome_db.sh --dry-run modified a file on disk"

	pass "rehome_db.sh --dry-run reports without modifying files"
}

test_rehome_db_reports_unresolvable_paths() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rehome-db-unresolved-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	new_dbdir="$tmp/new-home/ncbi-db"
	mkdir -p "$new_dbdir/Fungi"
	# No matching genome file exists anywhere under new_dbdir -- this
	# genome was never copied over, a real gap, not a path-prefix issue.
	printf '/home/h/old-project/ncbi-db/Fungi/GCF_missing.fna\n' > "$new_dbdir/.fungi"

	err="$tmp/stderr"
	if "$REPO_DIR/scripts/rehome_db.sh" "$new_dbdir" >/dev/null 2>"$err"; then
		fail "rehome_db.sh should exit non-zero when a path cannot be resolved"
	fi
	grep -Fq "GCF_missing.fna" "$err" || fail "rehome_db.sh did not report the unresolvable path"
	grep -Fq "/home/h/old-project/ncbi-db/Fungi/GCF_missing.fna" "$new_dbdir/.fungi" || fail "rehome_db.sh should leave an unresolvable path unchanged"

	pass "rehome_db.sh reports genuinely missing genomes instead of silently dropping them"
}

test_batch_classify_run_all() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-batch-classify-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir="$tmp/db"
	mkdir -p "$fake_home" "$dbdir/Viruses" "$dbdir/Protozoa" "$dbdir/taxonomy"
	ln -s "$REPO_DIR/scripts" "$fake_home/scripts"
	ln -s "$REPO_DIR/exe" "$fake_home/exe"
	ln -s "$REPO_DIR/batch-classify" "$fake_home/batch-classify"

	genome_v="$dbdir/Viruses/GCF_700000001.1_VirusX_genomic.fna"
	genome_p="$dbdir/Protozoa/GCF_700000002.1_ProtoY_genomic.fna"
	{
		printf '>NC_700000001.1 Virus X\n'
		python3 -c "import random; random.seed(11); print(''.join(random.choice('ACGT') for _ in range(3000)))"
	} > "$genome_v"
	{
		printf '>NC_700000002.1 Protozoan Y\n'
		python3 -c "import random; random.seed(12); print(''.join(random.choice('ACGT') for _ in range(3000)))"
	} > "$genome_p"

	printf '%s\n' "$genome_v" > "$dbdir/.viruses"
	printf '%s\n' "$genome_p" > "$dbdir/.protozoa"
	{
		printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
		printf 'viruses\tviral\tGCF_700000001.1\t700001\t700001\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_700000001.1_VirusX/GCF_700000001.1_VirusX_genomic.fna.gz\n'
	} > "$dbdir/.viruses.provenance.tsv"
	{
		printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
		printf 'protozoa\tprotozoa\tGCF_700000002.1\t700002\t700002\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_700000002.1_ProtoY/GCF_700000002.1_ProtoY_genomic.fna.gz\n'
	} > "$dbdir/.protozoa.provenance.tsv"
	printf '700001 | 2 | species |\n700002 | 2 | species |\n2 | 1 | superkingdom |\n' > "$dbdir/taxonomy/nodes.dmp"
	printf '1 | 1 |\n' > "$dbdir/taxonomy/merged.dmp"
	touch "$dbdir/.taxondata"

	results_dir="$tmp/results"
	CLARK_HOME="$fake_home" \
		"$REPO_DIR/batch-classify/run_all.sh" -d "$dbdir" -t viruses,protozoa \
			-x CLARK-l -n 2 -c 6 -l 100 \
			-o "$results_dir" > "$tmp/stdout" 2> "$tmp/stderr" ||
		fail "batch-classify/run_all.sh exited non-zero: $(cat "$tmp/stderr")"

	v_results="$results_dir/viruses/results.csv"
	p_results="$results_dir/protozoa/results.csv"
	require_file "$v_results"
	require_file "$p_results"

	[ "$(grep -c '^sample_' "$v_results")" -eq 6 ] || fail "batch-classify did not classify 6 viruses reads"
	[ "$(grep -c '^sample_' "$p_results")" -eq 6 ] || fail "batch-classify did not classify 6 protozoa reads"

	# A tiny synthetic genome has few discriminative k-mers, so an occasional
	# read can come back "NA" (no marker overlap); that's expected noise, not
	# a bug. What must never happen is a read assigned to the *other* type's
	# taxid, which would mean cross-contamination between iterations.
	if grep '^sample_' "$v_results" | awk -F',' '{ gsub(/ /, "", $3); if ($3 == "700002") exit 1 }'; then :; else
		fail "batch-classify's viruses run assigned a read to the protozoa taxid"
	fi
	grep -Fq ',700001' "$v_results" || fail "batch-classify did not correctly classify any viruses read"
	if grep '^sample_' "$p_results" | awk -F',' '{ gsub(/ /, "", $3); if ($3 == "700001") exit 1 }'; then :; else
		fail "batch-classify's protozoa run assigned a read to the viruses taxid"
	fi
	grep -Fq ',700002' "$p_results" || fail "batch-classify did not correctly classify any protozoa read"

	CLARK_HOME="$fake_home" \
		"$REPO_DIR/batch-classify/run_all.sh" -d "$dbdir" -t viruses \
			-x CLARK-l -n 2 -c 6 -l 100 \
			-o "$results_dir" -p > "$tmp/stdout-profile" 2> "$tmp/stderr-profile" ||
		fail "batch-classify/run_all.sh -p exited non-zero: $(cat "$tmp/stderr-profile")"
	grep -Fq "CUD_PROFILE" "$tmp/stdout-profile" || fail "batch-classify/run_all.sh -p did not print CLARK's CuD performance breakdown"

	pass "batch-classify/run_all.sh classifies multiple database types sharing one directory"
}

test_rapl_energy_measures_package_domains() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rapl-energy-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	sysfs="$tmp/powercap"
	mkdir -p "$sysfs/intel-rapl:0" "$sysfs/intel-rapl:0:0" "$sysfs/intel-rapl:1"
	printf 'package-0\n' > "$sysfs/intel-rapl:0/name"
	printf 'core\n' > "$sysfs/intel-rapl:0:0/name"
	printf 'package-1\n' > "$sysfs/intel-rapl:1/name"
	printf '1000000\n' > "$sysfs/intel-rapl:0/energy_uj"
	printf '500000\n' > "$sysfs/intel-rapl:0:0/energy_uj"
	printf '2000000\n' > "$sysfs/intel-rapl:1/energy_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:0/max_energy_range_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:1/max_energy_range_uj"

	# The wrapped command bumps both package counters (+300000 and +400000
	# uJ) but also the "core" sub-domain (+900000 uJ), which must NOT be
	# double-counted since it is already included in package-0's own total.
	out="$tmp/stdout"
	err="$tmp/stderr"
	python3 "$REPO_DIR/scripts/rapl_energy.py" --sysfs-dir "$sysfs" -- bash -c "
		printf '1300000\n' > '$sysfs/intel-rapl:0/energy_uj'
		printf '1400000\n' > '$sysfs/intel-rapl:0:0/energy_uj'
		printf '2400000\n' > '$sysfs/intel-rapl:1/energy_uj'
	" > "$out" 2> "$err" || fail "rapl_energy.py exited non-zero: $(cat "$err")"

	grep -Fq "pkg_joules=0.700000" "$out" || fail "rapl_energy.py did not sum only the package domains (expected 0.7 J): $(cat "$out")"
	grep -Fq "pkg_joules_scaled=0.490000" "$out" || fail "rapl_energy.py did not apply the default 0.70 scale factor: $(cat "$out")"
	grep -Fq "scale=0.70" "$out" || fail "rapl_energy.py did not report the scale factor used"
	grep -Fq "pkg_domains=2" "$out" || fail "rapl_energy.py did not report exactly 2 package domains (core sub-domain should be excluded)"
	grep -Fq "dram_joules=0.000000" "$out" || fail "rapl_energy.py should report zero dram energy when no dram domain exists"
	grep -Fq "dram_domains=0" "$out" || fail "rapl_energy.py should report zero dram domains when none exist"
	grep -Fq "total_joules=0.700000" "$out" || fail "rapl_energy.py's total_joules should equal pkg_joules when there is no dram domain"

	pass "rapl_energy.py sums package-only RAPL domains and reports raw + scaled-down energy"
}

test_rapl_energy_measures_dram_domain() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rapl-energy-dram-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	sysfs="$tmp/powercap"
	mkdir -p "$sysfs/intel-rapl:0" "$sysfs/intel-rapl:0:0" "$sysfs/intel-rapl:0:2"
	printf 'package-0\n' > "$sysfs/intel-rapl:0/name"
	printf 'core\n' > "$sysfs/intel-rapl:0:0/name"
	printf 'dram\n' > "$sysfs/intel-rapl:0:2/name"
	printf '1000000\n' > "$sysfs/intel-rapl:0/energy_uj"
	printf '500000\n' > "$sysfs/intel-rapl:0:0/energy_uj"
	printf '300000\n' > "$sysfs/intel-rapl:0:2/energy_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:0/max_energy_range_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:0:2/max_energy_range_uj"

	out="$tmp/stdout"
	python3 "$REPO_DIR/scripts/rapl_energy.py" --sysfs-dir "$sysfs" -- bash -c "
		printf '1300000\n' > '$sysfs/intel-rapl:0/energy_uj'
		printf '350000\n' > '$sysfs/intel-rapl:0:2/energy_uj'
	" > "$out" 2>/dev/null || fail "rapl_energy.py exited non-zero with a dram domain present"

	grep -Fq "pkg_joules=0.300000" "$out" || fail "rapl_energy.py did not correctly measure package energy alongside dram: $(cat "$out")"
	grep -Fq "dram_joules=0.050000" "$out" || fail "rapl_energy.py did not correctly measure the dram domain: $(cat "$out")"
	grep -Fq "dram_domains=1" "$out" || fail "rapl_energy.py did not count the dram domain"
	grep -Fq "total_joules=0.350000" "$out" || fail "rapl_energy.py's total_joules should be pkg_joules + dram_joules"
	grep -Fq "total_joules_scaled=0.245000" "$out" || fail "rapl_energy.py did not scale total_joules by the default 0.70 factor"

	pass "rapl_energy.py separately measures a dram RAPL domain when the platform exposes one"
}

test_rapl_energy_handles_counter_wraparound() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rapl-energy-wrap-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	sysfs="$tmp/powercap"
	mkdir -p "$sysfs/intel-rapl:0"
	printf 'package-0\n' > "$sysfs/intel-rapl:0/name"
	printf '99000000\n' > "$sysfs/intel-rapl:0/energy_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:0/max_energy_range_uj"

	out="$tmp/stdout"
	python3 "$REPO_DIR/scripts/rapl_energy.py" --sysfs-dir "$sysfs" -- bash -c "
		printf '500000\n' > '$sysfs/intel-rapl:0/energy_uj'
	" > "$out" 2>/dev/null || fail "rapl_energy.py exited non-zero on a wrapping counter"

	# (100000000 - 99000000) + 500000 = 1500000 uJ = 1.5 J
	grep -Fq "pkg_joules=1.500000" "$out" || fail "rapl_energy.py did not correctly account for counter wraparound: $(cat "$out")"

	pass "rapl_energy.py accounts for RAPL counter wraparound using max_energy_range_uj"
}

test_rapl_energy_unavailable_still_runs_command() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-rapl-energy-unavailable-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	out="$tmp/stdout"
	err="$tmp/stderr"
	python3 "$REPO_DIR/scripts/rapl_energy.py" --sysfs-dir "$tmp/does-not-exist" -- echo hello > "$out" 2> "$err"
	rc=$?

	[ "$rc" -eq 0 ] || fail "rapl_energy.py should propagate the wrapped command's exit code (echo exits 0)"
	grep -Fq "hello" "$out" || fail "rapl_energy.py did not run the wrapped command when RAPL is unavailable"
	grep -Fq "ENERGY_PROFILE" "$out" && fail "rapl_energy.py printed an ENERGY_PROFILE line despite RAPL being unavailable"
	grep -Fiq "unavailable" "$err" || fail "rapl_energy.py did not warn that RAPL energy accounting is unavailable"

	pass "rapl_energy.py falls back to running the command when RAPL is unreadable, with a warning"
}

test_batch_classify_run_all_energy_flag() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-batch-classify-energy-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir="$tmp/db"
	mkdir -p "$fake_home" "$dbdir/Viruses" "$dbdir/taxonomy" "$fake_home/exe"
	ln -s "$REPO_DIR/scripts" "$fake_home/scripts"
	ln -s "$REPO_DIR/batch-classify" "$fake_home/batch-classify"
	for binary in getTargetsDef getfilesToTaxNodes getAccssnTaxID; do
		ln -s "$REPO_DIR/exe/$binary" "$fake_home/exe/$binary"
	done
	ln -s "$REPO_DIR/exe/CLARK-l" "$fake_home/exe/CLARK-l"

	genome="$dbdir/Viruses/GCF_700000006.1_VirusR_genomic.fna"
	{
		printf '>NC_700000006.1 Virus R\n'
		python3 -c "import random; random.seed(51); print(''.join(random.choice('ACGT') for _ in range(3000)))"
	} > "$genome"
	printf '%s\n' "$genome" > "$dbdir/.viruses"
	{
		printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
		printf 'viruses\tviral\tGCF_700000006.1\t700006\t700006\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_700000006.1_VirusR/GCF_700000006.1_VirusR_genomic.fna.gz\n'
	} > "$dbdir/.viruses.provenance.tsv"
	printf '700006 | 2 | species |\n2 | 1 | superkingdom |\n' > "$dbdir/taxonomy/nodes.dmp"
	printf '1 | 1 |\n' > "$dbdir/taxonomy/merged.dmp"
	touch "$dbdir/.taxondata"

	sysfs="$tmp/powercap"
	mkdir -p "$sysfs/intel-rapl:0"
	printf 'package-0\n' > "$sysfs/intel-rapl:0/name"
	printf '0\n' > "$sysfs/intel-rapl:0/energy_uj"
	printf '100000000\n' > "$sysfs/intel-rapl:0/max_energy_range_uj"

	results_dir="$tmp/results"
	CLARK_HOME="$fake_home" RAPL_SYSFS_DIR="$sysfs" \
		"$REPO_DIR/batch-classify/run_all.sh" -d "$dbdir" -t viruses -x CLARK-l -n 2 -c 6 -l 100 \
			-o "$results_dir" -e > "$tmp/stdout" 2> "$tmp/stderr" ||
		fail "batch-classify/run_all.sh -e exited non-zero: $(cat "$tmp/stderr")"

	grep -Fq "ENERGY_PROFILE" "$tmp/stdout" || fail "batch-classify/run_all.sh -e did not print an ENERGY_PROFILE line: $(cat "$tmp/stdout")"
	require_file "$results_dir/viruses/results.csv"

	pass "batch-classify/run_all.sh -e wraps the classify step with RAPL energy measurement"
}

test_batch_classify_run_all_gpu_flag() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-batch-classify-gpu-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/fake-home"
	dbdir="$tmp/db"
	mkdir -p "$fake_home" "$dbdir/Viruses" "$dbdir/taxonomy" "$fake_home/exe"
	ln -s "$REPO_DIR/scripts" "$fake_home/scripts"
	ln -s "$REPO_DIR/batch-classify" "$fake_home/batch-classify"
	for binary in getTargetsDef getfilesToTaxNodes getAccssnTaxID; do
		ln -s "$REPO_DIR/exe/$binary" "$fake_home/exe/$binary"
	done

	genome="$dbdir/Viruses/GCF_700000005.1_VirusQ_genomic.fna"
	{
		printf '>NC_700000005.1 Virus Q\n'
		python3 -c "import random; random.seed(41); print(''.join(random.choice('ACGT') for _ in range(3000)))"
	} > "$genome"
	printf '%s\n' "$genome" > "$dbdir/.viruses"
	{
		printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
		printf 'viruses\tviral\tGCF_700000005.1\t700005\t700005\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_700000005.1_VirusQ/GCF_700000005.1_VirusQ_genomic.fna.gz\n'
	} > "$dbdir/.viruses.provenance.tsv"
	printf '700005 | 2 | species |\n2 | 1 | superkingdom |\n' > "$dbdir/taxonomy/nodes.dmp"
	printf '1 | 1 |\n' > "$dbdir/taxonomy/merged.dmp"
	touch "$dbdir/.taxondata"

	# Stand in for cuCLARK: record how it was invoked instead of actually
	# running a GPU classifier (none is available in this test environment).
	fake_cuclark="$tmp/cuCLARK"
	capture="$tmp/cuclark-capture.txt"
	cat > "$fake_cuclark" <<EOF
#!/usr/bin/env bash
echo "invoked with: \$*" > "$capture"
results=""
while [ "\$#" -gt 0 ]; do
	if [ "\$1" = "-R" ]; then
		results="\$2"
	fi
	shift
done
[ -n "\$results" ] && printf 'Object_ID, Length, Assignment\n' > "\$results.csv"
EOF
	chmod +x "$fake_cuclark"

	results_dir="$tmp/results"
	CLARK_HOME="$fake_home" CUCLARK_EXE="$fake_cuclark" \
		"$REPO_DIR/batch-classify/run_all.sh" -d "$dbdir" -t viruses -g -n 2 -c 3 -l 80 \
			-o "$results_dir" > "$tmp/stdout" 2> "$tmp/stderr" ||
		fail "batch-classify/run_all.sh -g exited non-zero: $(cat "$tmp/stderr")"

	require_file "$capture"
	grep -Fq -- "-k 31" "$capture" || fail "batch-classify/run_all.sh -g did not pass -k through to cuCLARK"
	grep -Fq -- "-n 2" "$capture" || fail "batch-classify/run_all.sh -g did not pass -n through to cuCLARK"

	pass "batch-classify/run_all.sh -g dispatches to the configured cuCLARK executable"
}

test_clark_cud_profile_opt_in() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-cud-profile-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	genome="$tmp/refA.fa"
	{
		printf '>refA\n'
		python3 -c "import random; random.seed(3); print(''.join(random.choice('ACGT') for _ in range(3000)))"
	} > "$genome"

	targets="$tmp/targets.txt"
	printf '%s\t111\n' "$genome" > "$targets"

	dbd="$tmp/dbd/"
	mkdir -p "$dbd"

	objects="$tmp/objects.fa"
	python3 -c "
seq = ''.join(l.strip() for l in open('$genome') if not l.startswith('>'))
print('>read1')
print(seq[100:250])
" > "$objects"

	"$REPO_DIR/exe/CLARK-l" -k 31 -T "$targets" -D "$dbd" -O "$objects" -R "$tmp/results1" -n 1 > "$tmp/out1" 2>&1
	if grep -Fq "CUD_PROFILE" "$tmp/out1"; then
		fail "CLARK-l printed CUD_PROFILE without CLARK_CUD_PROFILE=1 being set"
	fi

	CLARK_CUD_PROFILE=1 "$REPO_DIR/exe/CLARK-l" -k 31 -T "$targets" -D "$dbd" -O "$objects" -R "$tmp/results2" -n 1 > "$tmp/out2" 2>&1
	line="$(grep "CUD_PROFILE" "$tmp/out2" || true)"
	[ -n "$line" ] || fail "CLARK-l did not print a CUD_PROFILE line with CLARK_CUD_PROFILE=1"

	# CLARK-l ignores -k and always uses its own light k-mer size (27).
	printf '%s\n' "$line" | grep -Eq 'build_ns=[0-9.e+-]+ load_ns=[0-9.e+-]+ match_ns=[0-9.e+-]+ write_ns=[0-9.e+-]+ kmer=27 nbObjects=1$' \
		|| fail "CUD_PROFILE line has an unexpected format: $line"

	pass "CLARK_CUD_PROFILE opt-in prints a build/load/match/write breakdown only when requested"
}

test_get_targets_def_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-targets-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	lineage="$tmp/fileToTaxIDs.txt"
	out="$tmp/targets.txt"

	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
		"$tmp/refA.fa" "111" "111" "222" "333" "444" "555" "666" \
		"$tmp/refB.fa" "112" "112" "222" "333" "444" "555" "666" \
		> "$lineage"

	(
		cd "$tmp"
		"$REPO_DIR/exe/getTargetsDef" "$lineage" 0 > "$out"
	)

	grep -Fq "$tmp/refA.fa	111" "$out" || fail "getTargetsDef did not emit species target for refA"
	grep -Fq "$tmp/refB.fa	112" "$out" || fail "getTargetsDef did not emit species target for refB"
	pass "getTargetsDef helper emits expected target definitions"
}

test_get_targets_def_exit_code_ignores_excluded_count() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-targets-excluded-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	lineage="$tmp/fileToTaxIDs.txt"
	out="$tmp/targets.txt"

	{
		printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tmp/refA.fa" "111" "111" "222" "333" "444" "555" "666"
		printf '%s\t-1\n' "$tmp/refB.fa"
		printf '%s\t-1\n' "$tmp/refC.fa"
	} > "$lineage"

	set +e
	(
		cd "$tmp"
		"$REPO_DIR/exe/getTargetsDef" "$lineage" 0 > "$out"
	)
	rc=$?
	set -e

	[ "$rc" -eq 0 ] || fail "getTargetsDef exited non-zero ($rc) when some files were excluded; this makes scripts/set_targets.sh (set -e) abort silently with no error message"
	grep -Fq "$tmp/refA.fa	111" "$out" || fail "getTargetsDef did not emit the target for the one resolvable file"
	grep -Fq "$tmp/refB.fa" "$tmp/files_excluded.txt" || fail "getTargetsDef did not record refB as excluded"
	grep -Fq "$tmp/refC.fa" "$tmp/files_excluded.txt" || fail "getTargetsDef did not record refC as excluded"
	pass "getTargetsDef exits 0 regardless of how many files were excluded"
}

test_get_accssn_taxid_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-accession-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	files="$tmp/files.txt"
	accessions="$tmp/nucl_accession2taxid"
	merged="$tmp/merged.dmp"
	ref_a="$tmp/refA.fa"
	ref_b="$tmp/refB.fa"
	out="$tmp/file_taxid.txt"

	printf '>NC_000001.1 synthetic reference A\nACGT\n' > "$ref_a"
	printf '>NC_000002.1 synthetic reference B\nTGCA\n' > "$ref_b"
	printf '%s\n%s\n' "$ref_a" "$ref_b" > "$files"
	printf 'NC_000001\tNC_000001.1\t111\t111\n' > "$accessions"
	printf '111 | 211 |\n' > "$merged"

	"$REPO_DIR/exe/getAccssnTaxID" "$files" "$accessions" "$merged" > "$out" 2> "$tmp/stderr"

	grep -Fq "$ref_a	NC_000001	211" "$out" || fail "getAccssnTaxID did not map and merge refA taxid"
	grep -Fq "$ref_b	NC_000002	-1" "$out" || fail "getAccssnTaxID did not report unmapped refB"
	pass "getAccssnTaxID maps accession IDs for tiny FASTA inputs"
}

test_getfiles_to_taxnodes_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-lineage-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	nodes="$tmp/nodes.dmp"
	file_taxid="$tmp/file_taxid.txt"
	out="$tmp/fileToTaxIDs.txt"

	cat > "$nodes" <<'EOF'
1 | 1 | root |
60 | 1 | phylum |
50 | 60 | class |
40 | 50 | order |
30 | 40 | family |
20 | 30 | genus |
10 | 20 | species |
EOF
	printf '%s\t%s\t%s\n' "$tmp/refA.fa" "NC_000001" "10" > "$file_taxid"

	"$REPO_DIR/exe/getfilesToTaxNodes" "$nodes" "$file_taxid" > "$out" 2> "$tmp/stderr"

	grep -Fq "$tmp/refA.fa	10	10	20	30	40	50	UNKNOWN" "$out" || fail "getfilesToTaxNodes did not emit expected lineage"
	pass "getfilesToTaxNodes expands a synthetic taxonomy lineage"
}

test_exe_seq_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-exeseq-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	multifasta="$tmp/multi.fa"
	outdir="$tmp/split"
	mkdir -p "$outdir"
	cat > "$multifasta" <<'EOF'
>seqA
ACGT
>seqB
TGCA
EOF

	"$REPO_DIR/exe/exeSeq" "$multifasta" "$outdir" > "$tmp/stdout" 2> "$tmp/stderr"

	[ "$(find "$outdir" -type f -name '*.fa' | wc -l | tr -d ' ')" = "2" ] || fail "exeSeq did not split two FASTA records"
	grep -R -Fq "ACGT" "$outdir" || fail "exeSeq output does not contain seqA bases"
	grep -R -Fq "TGCA" "$outdir" || fail "exeSeq output does not contain seqB bases"
	pass "exeSeq splits a multi-FASTA file"
}

test_dscript_maker_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-dscript-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	assembly="$tmp/assembly_summary.txt"
	out="$tmp/download.sh"
	printf 'https://example.org/refseq/GCF_000001405.40\n' > "$assembly"

	"$REPO_DIR/exe/dscriptMaker" "$assembly" > "$out"

	grep -Fq "wget https://example.org/refseq/GCF_000001405.40/GCF_000001405.40_genomic.fna.gz" "$out" || fail "dscriptMaker did not emit expected download command"
	pass "dscriptMaker emits a deterministic download command"
}

test_density_helpers_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-density-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	result="$tmp/results.csv"
	gamma_out="$tmp/gamma.txt"
	conf_out="$tmp/confidence.txt"
	cat > "$result" <<'EOF'
Object_ID,Length,Gamma,1st_assignment,score1,2nd_assignment,score2,confidence
read1,100,0.50,111,9,222,3,0.80
read2,100,1.20,111,8,222,4,1.00
EOF

	"$REPO_DIR/exe/getGammaDensity" "$result" > "$gamma_out" 2> "$tmp/gamma.err"
	"$REPO_DIR/exe/getConfidenceDensity" "$result" > "$conf_out" 2> "$tmp/conf.err"

	grep -Fq "assignments with Gamma score found" "$tmp/gamma.err" || fail "getGammaDensity did not process gamma scores"
	grep -Fq "[0.50,0.52[" "$gamma_out" || fail "getGammaDensity did not report expected 0.50 gamma bucket"
	grep -Fq "[>=1]" "$gamma_out" || fail "getGammaDensity did not report the >=1 gamma bucket"
	grep -Fq "assignments with confidence score found" "$tmp/conf.err" || fail "getConfidenceDensity did not process confidence scores"
	grep -Fq "[0.80,0.82[" "$conf_out" || fail "getConfidenceDensity did not report expected 0.80 confidence bucket"
	grep -Fq "[1]" "$conf_out" || fail "getConfidenceDensity did not report the exact confidence score 1 bucket"
	grep -Fq "[4.00,4.02[" "$conf_out" && fail "getConfidenceDensity emitted an impossible interval"
	grep -Fq "Interval" "$conf_out" || fail "getConfidenceDensity did not emit a density table"
	pass "density helpers summarize synthetic CLARK scores"
}

test_extract_seqs_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-extract-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	reads="$tmp/reads.fq"
	results="$tmp/results.csv"
	out_prefix="$tmp/extracted"
	cat > "$reads" <<'EOF'
@read1
ACGT
+
!!!!
@read2
TGCA
+
!!!!
EOF
	cat > "$results" <<'EOF'
Object_ID,Length,Assignment
read1,4,111
read2,4,222
EOF

	"$REPO_DIR/exe/extractSeqs" 111 "$reads" "$results" "$out_prefix" > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq "@read1" "$out_prefix.fq" || fail "extractSeqs did not extract the matching read"
	if grep -Fq "@read2" "$out_prefix.fq"; then
		fail "extractSeqs extracted a read assigned to a different taxid"
	fi
	pass "extractSeqs extracts matching FASTQ records"
}

test_extract_seqs_long_report_thresholds() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-extract-long-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	reads="$tmp/reads.fa"
	results="$tmp/results.csv"
	out_prefix="$tmp/extracted"
	cat > "$reads" <<'EOF'
>read1 description
ACGT
>read2
TGCA
>read3
CCCC
EOF
	cat > "$results" <<'EOF'
Object_ID,Length,Gamma,1st_assignment,score1,2nd_assignment,score2,confidence
read1,4,0.50,111,9,222,3,0.80
read2,4,0.01,111,9,222,3,0.90
read3,4,0.50,222,9,111,3,0.90
EOF

	"$REPO_DIR/exe/extractSeqs" 111 "$reads" "$results" "$out_prefix" 0.03 0.75 > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq ">read1 description" "$out_prefix.fa" || fail "extractSeqs did not keep FASTA header for passing long-report assignment"
	grep -Fq "ACGT" "$out_prefix.fa" || fail "extractSeqs did not keep FASTA sequence for passing long-report assignment"
	if grep -Fq ">read2" "$out_prefix.fa" || grep -Fq ">read3" "$out_prefix.fa"; then
		fail "extractSeqs kept filtered or wrong-taxid FASTA records"
	fi
	pass "extractSeqs applies gamma/confidence thresholds for long reports"
}

test_get_abundance_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-abundance-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	results="$tmp/results.csv"
	out="$tmp/abundance.csv"
	cat > "$results" <<'EOF'
Object_ID,Length,Assignment
read1,4,111
read2,4,111
read3,4,NA
EOF

	(
		cd "$tmp"
		"$REPO_DIR/exe/getAbundance" -F "$results" > "$out"
	)

	grep -Fq "Name,TargetID,Count,Proportion_All(%),Proportion_Classified(%)" "$out" || fail "getAbundance did not emit expected header"
	grep -Fq "111,111,2," "$out" || fail "getAbundance did not count assigned reads"
	grep -Fq "UNKNOWN,UNKNOWN,1," "$out" || fail "getAbundance did not count unassigned reads"
	pass "getAbundance summarizes a tiny assignment file"
}

test_get_abundance_filters_extended_scores() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-abundance-filter-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	results="$tmp/results.csv"
	out="$tmp/abundance.csv"
	cat > "$results" <<'EOF'
Object_ID,Length,Gamma,1st_assignment,score1,2nd_assignment,score2,confidence
read1,100,0.50,111,9,222,3,0.80
read2,100,0.01,222,9,111,3,0.90
read3,100,0.50,333,9,111,3,0.60
read4,100,0.50,NA,0,NA,0,1.00
EOF

	(
		cd "$tmp"
		"$REPO_DIR/exe/getAbundance" -F "$results" -c 0.75 -g 0.03 > "$out"
	)

	grep -Fq "111,111,1,25,100" "$out" || fail "getAbundance did not keep the high-confidence assignment"
	grep -Fq "UNKNOWN,UNKNOWN,3,75,-" "$out" || fail "getAbundance did not move filtered assignments to UNKNOWN"
	if "$REPO_DIR/exe/getAbundance" -c 0.80 > "$tmp/no_file.out" 2> "$tmp/no_file.err"; then
		fail "getAbundance accepted missing -F results"
	fi
	grep -Fq "Please provide one (or several) CLARK output file" "$tmp/no_file.err" ||
		fail "getAbundance did not explain missing -F results"
	pass "getAbundance filters extended scores and validates required input"
}

test_make_summary_tables_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-summary-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	report_a="$tmp/sampleA.csv"
	report_b="$tmp/sampleB.csv"
	cat > "$report_a" <<'EOF'
Name,TaxID,Count,Proportion_All(%),Proportion_Classified(%)
TaxonA,111,3,75,100
UNKNOWN,UNKNOWN,1,25,-
EOF
	cat > "$report_b" <<'EOF'
Name,TaxID,Count,Proportion_All(%),Proportion_Classified(%)
TaxonB,222,2,50,100
UNKNOWN,UNKNOWN,2,50,-
EOF

	(
		cd "$tmp"
		"$REPO_DIR/exe/makeSummaryTables" 2 0 "$report_a" "$report_b" > stdout 2> stderr
	)

	grep -Fq "sampleA" "$tmp/TableSummary_per_Report.csv" || fail "makeSummaryTables did not include sampleA"
	grep -Fq "TaxonA" "$tmp/TableSummary_HitCount.csv" || fail "makeSummaryTables did not include TaxonA"
	grep -Fq "TaxonB" "$tmp/TableSummary_HitCount.csv" || fail "makeSummaryTables did not include TaxonB"
	pass "makeSummaryTables writes summary tables for tiny reports"
}

test_make_summary_tables_all_unknown() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-summary-empty-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	report="$tmp/all_unknown.csv"
	cat > "$report" <<'EOF'
Name,TaxID,Count,Proportion_All(%),Proportion_Classified(%)
UNKNOWN,UNKNOWN,4,100,-
EOF

	(
		cd "$tmp"
		"$REPO_DIR/exe/makeSummaryTables" 3 0 "$report" > stdout 2> stderr
	)

	grep -Fq "all_unknown,4,0," "$tmp/TableSummary_per_Report.csv" || fail "makeSummaryTables did not summarize an all-UNKNOWN report"
	grep -Fq "#TotalReads,4," "$tmp/TableSummary_HitCount.csv" || fail "makeSummaryTables did not write total reads for all-UNKNOWN report"
	grep -Fq "#TotalReadsMapped,0," "$tmp/TableSummary_HitCount.csv" || fail "makeSummaryTables did not write mapped reads for all-UNKNOWN report"
	pass "makeSummaryTables handles reports with no classified taxa"
}

test_target_specific_kmers_stat_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-kmer-stat-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	fake_home="$tmp/home"
	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	label_file="$dbdir/db_central_k2_t2_s1610612741_m0.tsk.lb"
	mkdir -p "$dbdir" "$fake_home/exe"
	ln -s "$REPO_DIR/exe/getTargetSpecificKmersStat" "$fake_home/exe/getTargetSpecificKmersStat"
	printf '%s\t%s\n%s\t%s\n' "$tmp/refA.fa" "111" "$tmp/refB.fa" "222" > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir" > "$settings"
	cp "$settings" "$fake_home/.settings"
	printf '\000\000\001\000\001\000' > "$label_file"

	(
		cd "$tmp"
		"$REPO_DIR/exe/getTargetSpecificKmersStat" "$settings" 2 0 > stdout 2> stderr
	)

	grep -Fq "111,1," "$tmp/targets.distribution.csv" || fail "getTargetSpecificKmersStat did not count target 111"
	grep -Fq "222,2," "$tmp/targets.distribution.csv" || fail "getTargetSpecificKmersStat did not count target 222"
	rm "$tmp/targets.distribution.csv"
	(
		cd "$tmp"
		CLARK_HOME="$fake_home" "$REPO_DIR/scripts/getTargetsKmers_distribution.sh" 2 > stdout.script 2> stderr.script
	)
	grep -Fq "111,1," "$tmp/targets.distribution.csv" || fail "getTargetsKmers_distribution.sh did not default min frequency to 0"
	pass "getTargetSpecificKmersStat counts labels in a tiny boundary-k database"
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
	if grep -R "ftp://ftp.ncbi.nih.gov\|ftp://ftp.ncbi.nlm.nih.gov" "$REPO_DIR/scripts/download_RefSeqDB.sh" "$REPO_DIR/scripts/download_taxondata.sh" >/dev/null; then
		fail "legacy or misspelled NCBI FTP URL remains"
	fi
	grep -Fq "https://ftp.ncbi.nlm.nih.gov" "$REPO_DIR/scripts/download_taxondata.sh" || fail "taxonomy downloader does not use HTTPS NCBI URL"
	pass "NCBI download URLs use HTTPS host"
}

test_portable_script_paths() {
	if grep -R "[r]eadlink -f" "$REPO_DIR" --include='*.sh' --exclude-dir='.git' --exclude-dir='build' --exclude-dir='exe' >/dev/null; then
		fail "non-portable shell path resolver remains"
	fi
	pass "shell scripts use portable path resolution"
}

test_shell_syntax
test_no_root_shell_scripts
test_required_executables
test_version_binaries
test_file_hash_unit_binary
test_kso_requires_spectrum_mode
test_classify_wrapper_quotes_paths
test_classify_wrapper_gzip
test_classify_wrapper_paired_light_variant
test_classify_wrapper_rejects_conflicting_variants
test_scripts_directory_entrypoint
test_set_targets_records_absolute_db_paths
test_set_targets_passes_refseq_download_options
test_set_targets_defaults_to_parallel_resume_downloads
test_update_taxonomy_requires_db_directory
test_refseq_downloader_dry_run_manifest
test_refseq_downloader_normalizes_ncbi_trailing_slash_paths
test_refseq_downloader_skips_non_url_ftp_path_rows
test_refseq_downloader_rejects_malformed_urls_before_download
test_refseq_downloader_tty_progress_rewrites_line
test_refseq_downloader_quiet_aggregate_progress
test_refseq_downloader_retries_parallel_transient_failures
test_refseq_downloader_defers_retry_until_after_first_pass
test_refseq_downloader_fails_when_deferred_downloads_remain
test_refseq_downloader_reports_early_parallel_progress
test_refseq_downloader_success_matrix
test_refseq_downloader_filters_and_resumes
test_make_metadata_uses_refseq_provenance_taxids
test_documentation_script_paths
test_make_sample_smoke
test_make_sample_requires_configured_targets
test_make_benchmark_reads_hiseq_profile
test_make_benchmark_reads_error_injection
test_eval_accuracy_smoke
test_rehome_db_fixes_stale_absolute_paths
test_rehome_db_dry_run_leaves_files_untouched
test_rehome_db_reports_unresolvable_paths
test_batch_classify_run_all
test_batch_classify_run_all_gpu_flag
test_rapl_energy_measures_package_domains
test_rapl_energy_measures_dram_domain
test_rapl_energy_handles_counter_wraparound
test_rapl_energy_unavailable_still_runs_command
test_batch_classify_run_all_energy_flag
test_clark_cud_profile_opt_in
test_get_targets_def_smoke
test_get_targets_def_exit_code_ignores_excluded_count
test_get_accssn_taxid_smoke
test_getfiles_to_taxnodes_smoke
test_exe_seq_smoke
test_dscript_maker_smoke
test_density_helpers_smoke
test_extract_seqs_smoke
test_extract_seqs_long_report_thresholds
test_get_abundance_smoke
test_get_abundance_filters_extended_scores
test_make_summary_tables_smoke
test_make_summary_tables_all_unknown
test_target_specific_kmers_stat_smoke
test_clark_l_label_bug_regression
test_ncbi_urls
test_portable_script_paths
