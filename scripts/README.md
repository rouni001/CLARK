# CLARK Scripts

This directory contains the implementation scripts used to install CLARK,
prepare databases, classify reads, estimate abundance, and run maintenance
tasks.

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

The downloader writes `.<database>.download_manifest.tsv` and
`.<database>.provenance.tsv` in the database directory. These files record the
planned or completed URLs, source accessions, and RefSeq source dates when they
are available from NCBI assembly summaries.
