# CLARK Quickstart for Metagenomics and Genomics Researchers

This guide is for researchers who want to classify reads or contigs but do not
want to debug compiler flags, shell quoting, or hidden CLARK state files.

## 1. Install

```sh
cd /path/to/CLARK
make all
make test
```

`make all` builds the CLARK helper tools plus the three classifiers:

- `exe/CLARK`
- `exe/CLARK-l`
- `exe/CLARK-S`

If your compiler does not support OpenMP, the build still succeeds but prints
that CLARK was built in single-threaded mode. On macOS with Apple clang, install
a compiler with OpenMP support and rerun with `CXX=g++ make all` if you need
parallel execution.

## 2. Prepare a Database

Choose a database directory with enough disk space. RefSeq bacterial databases
can be very large.

```sh
scripts/set_targets.sh /path/to/clark-db bacteria viruses --species
```

This writes the target configuration into CLARK's local `.settings` file and
stores database-specific files under `/path/to/clark-db`. CLARK records this
database directory as an absolute path so later classification and maintenance
commands are not dependent on your current working directory.

For custom references, place FASTA files in:

```text
/path/to/clark-db/Custom/
```

Then run:

```sh
scripts/set_targets.sh /path/to/clark-db custom --species
```

## 3. Classify Reads or Contigs

Single input:

```sh
scripts/classify_metagenome.sh -O sample.fastq -R sample.results.csv -m 2 -n 8
```

Paired-end input:

```sh
scripts/classify_metagenome.sh -P sample_R1.fastq sample_R2.fastq -R sample.results.csv -m 2 -n 8
```

Gzipped input:

```sh
scripts/classify_metagenome.sh -O sample.fastq.gz -R sample.results.csv --gzipped -m 2 -n 8
```

Use `--light` for CLARK-l and `--spaced` for CLARK-S.

## 4. Keep Runs Reproducible

Record these details with each analysis:

- CLARK version: `exe/CLARK --version`
- Git commit or release archive used
- Database directory path
- Date when taxonomy/reference data were downloaded
- Database choices and taxonomy rank passed to `scripts/set_targets.sh`
- CLARK command line used for classification
- RefSeq download manifest/provenance files written as
  `.<database>.download_manifest.tsv` and `.<database>.provenance.tsv`

## 5. Script Layout

All repository shell scripts live under `scripts/`. Build CLARK with `make all`,
then run commands such as `scripts/set_targets.sh` and
`scripts/classify_metagenome.sh` directly from the repository root. From another
directory, use an absolute script path such as
`/path/to/CLARK/scripts/classify_metagenome.sh`.

## 6. Troubleshooting

- If `scripts/classify_metagenome.sh` says targets are not configured, run
  `scripts/set_targets.sh` first.
- If a script says an executable is missing, run `make all`.
- If `-n` does not speed up classification, check the build output for
  `OpenMP: disabled`.
- Paths containing spaces are supported by the modernized scripts. Avoid moving
  the database directory after running `scripts/set_targets.sh`; rerun
  `scripts/set_targets.sh` if the database path changes.
- To preview a large RefSeq download before starting it, run
  `scripts/download_RefSeqDB.sh --dry-run <DIR_DB/> <database>`.
