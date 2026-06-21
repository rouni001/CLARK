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
			*"Version: 1.4.1.0-a"*) ;;
			*) fail "unexpected version output from $binary: $version" ;;
		esac
	done
	pass "CLARK variant version executables run"
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

test_refseq_downloader_mocked_assembly_summary() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-refseq-mocked-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	dbdir="$tmp/db"
	bindir="$tmp/bin"
	log="$tmp/wget.log"
	mkdir -p "$bindir"

	cat > "$bindir/wget" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
if [ "$1" = "-O" ]; then
	out="$2"
	url="$3"
	printf '%s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	case "$url" in
		*/assembly_summary.txt)
			{
				printf '# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\tsubmitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\n'
				printf 'GCF_999999999.1\tna\tna\tna\trepresentative genome\t10239\t10239\tMock virus\tna\tna\tlatest\tComplete Genome\tMajor\tFull\t2024-01-02\tMockVirus1\tCLARK\tna\tna\thttps://example.org/refseq/GCF_999999999.1_MockVirus1\n'
				printf 'GCF_000000000.1\tna\tna\tna\trepresentative genome\t10239\t10239\tOld virus\tna\tna\treplaced\tComplete Genome\tMajor\tFull\t2020-01-02\tOldVirus\tCLARK\tna\tna\thttps://example.org/refseq/GCF_000000000.1_OldVirus\n'
			} > "$out"
			;;
		*)
			echo "unexpected wget -O URL: $url" >&2
			exit 1
			;;
	esac
else
	url="$1"
	printf '%s\n' "$url" >> "$CLARK_TEST_WGET_LOG"
	file="${url##*/}"
	printf '>mock-virus\nACGT\n' | gzip > "$file"
fi
STUB
	chmod +x "$bindir/wget"

	PATH="$bindir:$PATH" \
	CLARK_TEST_WGET_LOG="$log" \
		"$REPO_DIR/scripts/download_RefSeqDB.sh" "$dbdir" viruses > "$tmp/downloader.out"

	require_file "$dbdir/.viruses"
	grep -Fq "GCF_999999999.1_MockVirus1_genomic.fna" "$dbdir/.viruses" || fail "RefSeq downloader did not list decompressed virus FASTA"
	grep -Fq "GCF_999999999.1" "$dbdir/.viruses.provenance.tsv" || fail "RefSeq downloader did not record assembly accession provenance"
	grep -Fq "2024-01-02" "$dbdir/.viruses.provenance.tsv" || fail "RefSeq downloader did not preserve assembly source date"
	grep -Fq "downloaded" "$dbdir/.viruses.download_manifest.tsv" || fail "RefSeq downloader did not record completed downloads"
	[ ! -e "$dbdir/Viruses/download.sh" ] || fail "RefSeq downloader generated a legacy download.sh script"
	pass "RefSeq downloader uses mocked assembly summaries without generated scripts"
}

test_documentation_script_paths() {
	if grep -E "\./(buildSpacedDB|classify_metagenome|clean|download_RefSeqDB|download_taxondata|estimate_abundance|evaluate_density_confidence|evaluate_density_gamma|extractSequences|getTargetsKmers_distribution|install|makeSummaryTables|make_metadata|resetCustomDB|set_targets|updateTaxonomy)\.sh|\./scripts/" \
		"$REPO_DIR/README.md" "$REPO_DIR/README_FULL.md" "$REPO_DIR/docs/QUICKSTART.md" "$REPO_DIR/scripts/README.md" >/dev/null; then
		fail "documentation still contains root-relative script examples"
	fi
	pass "README documents scripts/ commands without root-relative script paths"
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
	grep -Fq "[>=1]" "$gamma_out" || fail "getGammaDensity did not report the >=1 gamma bucket"
	grep -Fq "assignments with confidence score found" "$tmp/conf.err" || fail "getConfidenceDensity did not process confidence scores"
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

	"$REPO_DIR/exe/extractSeqs" 111 "$reads" "$results" "$out_prefix" 0 0 > "$tmp/stdout" 2> "$tmp/stderr"

	grep -Fq "@read1" "$out_prefix.fq" || fail "extractSeqs did not extract the matching read"
	if grep -Fq "@read2" "$out_prefix.fq"; then
		fail "extractSeqs extracted a read assigned to a different taxid"
	fi
	pass "extractSeqs extracts matching FASTQ records"
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

test_target_specific_kmers_stat_smoke() {
	tmp="$(mktemp -d "${TMPDIR:-/tmp}/clark-kmer-stat-test.XXXXXX")"
	trap 'rm -rf "$tmp"' RETURN

	settings="$tmp/settings"
	targets="$tmp/targets.txt"
	dbdir="$tmp/db"
	label_file="$dbdir/db_central_k3_t2_s1610612741_m0.tsk.lb"
	mkdir -p "$dbdir"
	printf '%s\t%s\n%s\t%s\n' "$tmp/refA.fa" "111" "$tmp/refB.fa" "222" > "$targets"
	printf -- '-T %s\n-D %s\n' "$targets" "$dbdir" > "$settings"
	printf '\000\000\001\000\001\000' > "$label_file"

	(
		cd "$tmp"
		"$REPO_DIR/exe/getTargetSpecificKmersStat" "$settings" 3 0 > stdout 2> stderr
	)

	grep -Fq "111,1," "$tmp/targets.distribution.csv" || fail "getTargetSpecificKmersStat did not count target 111"
	grep -Fq "222,2," "$tmp/targets.distribution.csv" || fail "getTargetSpecificKmersStat did not count target 222"
	pass "getTargetSpecificKmersStat counts labels in a tiny database"
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
test_classify_wrapper_quotes_paths
test_classify_wrapper_gzip
test_classify_wrapper_paired_light_variant
test_classify_wrapper_rejects_conflicting_variants
test_scripts_directory_entrypoint
test_set_targets_records_absolute_db_paths
test_update_taxonomy_requires_db_directory
test_refseq_downloader_dry_run_manifest
test_refseq_downloader_mocked_assembly_summary
test_documentation_script_paths
test_get_targets_def_smoke
test_get_accssn_taxid_smoke
test_getfiles_to_taxnodes_smoke
test_exe_seq_smoke
test_dscript_maker_smoke
test_density_helpers_smoke
test_extract_seqs_smoke
test_get_abundance_smoke
test_make_summary_tables_smoke
test_target_specific_kmers_stat_smoke
test_clark_l_label_bug_regression
test_ncbi_urls
test_portable_script_paths
