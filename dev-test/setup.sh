#!/usr/bin/env bash
# Prepares dev-test/db so scripts/set_targets.sh can build targets against it
# without any network access. Run this once before set_targets.sh.

set -euo pipefail

DEV_TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="$DEV_TEST_DIR/db"

genomeA="$DB_DIR/Bacteria/GCF_900000001.1_DemoOrgA_genomic.fna"
genomeB="$DB_DIR/Bacteria/GCF_900000002.1_DemoOrgB_genomic.fna"

[ -f "$genomeA" ] || { echo "Error: missing $genomeA" >&2; exit 1; }
[ -f "$genomeB" ] || { echo "Error: missing $genomeB" >&2; exit 1; }

printf '%s\n%s\n' "$genomeA" "$genomeB" > "$DB_DIR/.bacteria"

{
	printf 'database\tsource\taccession\ttaxid\tspecies_taxid\tseq_rel_date\tassembly_level\tversion_status\turl\n'
	printf 'bacteria\tbacteria\tGCF_900000001.1\t900001\t900001\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_900000001.1_DemoOrgA/GCF_900000001.1_DemoOrgA_genomic.fna.gz\n'
	printf 'bacteria\tbacteria\tGCF_900000002.1\t900002\t900002\t2026-01-01\tComplete Genome\tlatest\thttps://example.org/refseq/GCF_900000002.1_DemoOrgB/GCF_900000002.1_DemoOrgB_genomic.fna.gz\n'
} > "$DB_DIR/.bacteria.provenance.tsv"

touch "$DB_DIR/.taxondata"

echo "dev-test database metadata ready in $DB_DIR"
echo "Next steps:"
echo "  CLARK_HOME=\"\$(pwd)\" scripts/set_targets.sh $DB_DIR bacteria --species"
echo "  CLARK_HOME=\"\$(pwd)\" scripts/classify_metagenome.sh -O $DEV_TEST_DIR/sample.fa -R $DEV_TEST_DIR/result --light"
