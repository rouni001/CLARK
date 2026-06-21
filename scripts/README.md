# CLARK Scripts

This directory contains the implementation scripts used to prepare databases,
classify reads, estimate abundance, and run maintenance tasks.

The repository root intentionally does not keep duplicate script launchers. Use
these scripts directly from this directory so the top-level project stays
focused on source, docs, tests, and build metadata.

All scripts resolve the CLARK repository root through `CLARK_HOME` when it is
set, or by using this directory's parent. Commands written as `scripts/...`
assume the current directory is the repository root. From another directory, use
the absolute script path. For example:

```sh
scripts/classify_metagenome.sh ...
/path/to/CLARK/scripts/classify_metagenome.sh ...
```

Before downloading large RefSeq datasets, you can inspect the planned work:

```sh
scripts/download_RefSeqDB.sh --dry-run /path/to/clark-db bacteria
```

RefSeq downloads use 8 parallel workers and resume mode by default. To override
the worker count, pass `--download-threads` to the user-facing setup script:

```sh
scripts/set_targets.sh /path/to/clark-db bacteria --download-threads 4
```

For faster exploratory databases, restrict RefSeq assembly summaries to
representative or reference genomes:

```sh
scripts/set_targets.sh /path/to/clark-db bacteria --refseq-category representative
```

The downloader writes `.<database>.download_manifest.tsv` and
`.<database>.provenance.tsv` in the database directory. These files record the
planned or completed URLs, source accessions, and RefSeq source dates when they
are available from NCBI assembly summaries. Standard RefSeq libraries with
assembly-summary taxids can use this provenance directly, avoiding the large
global accession-map lookup when it is not needed.
