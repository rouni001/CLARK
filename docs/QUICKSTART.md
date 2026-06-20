# CLARK Quickstart for Metagenomics and Genomics Researchers

This guide is for researchers who want to classify reads or contigs but do not
want to debug compiler flags, shell quoting, or hidden CLARK state files.

## 1. Install

```sh
./install.sh
make test
```

`install.sh` builds the CLARK helper tools plus the three classifiers:

- `exe/CLARK`
- `exe/CLARK-l`
- `exe/CLARK-S`

If your compiler does not support OpenMP, the build still succeeds but prints
that CLARK was built in single-threaded mode. On macOS with Apple clang, install
a compiler with OpenMP support and rerun with `CXX=g++ ./install.sh` if you need
parallel execution.

## 2. Prepare a Database

Choose a database directory with enough disk space. RefSeq bacterial databases
can be very large.

```sh
./set_targets.sh /path/to/clark-db bacteria viruses --species
```

This writes the target configuration into CLARK's local `.settings` file and
stores database-specific files under `/path/to/clark-db`.

For custom references, place FASTA files in:

```text
/path/to/clark-db/Custom/
```

Then run:

```sh
./set_targets.sh /path/to/clark-db custom --species
```

## 3. Classify Reads or Contigs

Single input:

```sh
./classify_metagenome.sh -O sample.fastq -R sample.results.csv -m 2 -n 8
```

Paired-end input:

```sh
./classify_metagenome.sh -P sample_R1.fastq sample_R2.fastq -R sample.results.csv -m 2 -n 8
```

Gzipped input:

```sh
./classify_metagenome.sh -O sample.fastq.gz -R sample.results.csv --gzipped -m 2 -n 8
```

Use `--light` for CLARK-l and `--spaced` for CLARK-S.

## 4. Keep Runs Reproducible

Record these details with each analysis:

- CLARK version: `exe/CLARK --version`
- Git commit or release archive used
- Database directory path
- Date when taxonomy/reference data were downloaded
- Database choices and taxonomy rank passed to `set_targets.sh`
- CLARK command line used for classification

## 5. Troubleshooting

- If `classify_metagenome.sh` says targets are not configured, run
  `set_targets.sh` first.
- If a wrapper says an executable is missing, run `./install.sh`.
- If `-n` does not speed up classification, check the install output for
  `OpenMP: disabled`.
- Paths containing spaces are supported by the modernized wrappers, but avoid
  moving the database directory after running `set_targets.sh`; rerun
  `set_targets.sh` if the path changes.
