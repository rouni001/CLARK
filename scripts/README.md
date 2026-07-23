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

## Moving a populated database directory to another machine

`set_targets.sh` caches several hidden metadata files under the database
directory (`.protozoa`, `.protozoa.fileToAccssnTaxID`,
`.protozoa.fileToTaxIDs`, one triple per database type). The first field on
every line of these files is an **absolute path** to a genome FASTA file,
recorded when the sequences were first downloaded/scanned. If you `rsync`,
`scp`, or otherwise copy the whole database directory to a different
machine (or a different path on the same machine), those cached files still
point at the *old* absolute paths -- `set_targets.sh` sees the cache is
non-empty and skips re-downloading/re-scanning, and CLARK then looks for
genomes at a path that no longer exists there.

Symptom: `set_targets.sh`/`batch-classify/run_all.sh` prints "Collecting
metadata... done." almost instantly, then a stream of `Warning: failed to
read /old/host/path/...: [Errno 2] No such file or directory` from
`make_sample.py`, followed by `no target genome is long enough for a
N-base read` (or `getTargetsDef` excluding every file).

Fix: after copying the directory, run this once (no network access
needed):

```sh
scripts/rehome_db.sh /path/to/copied/ncbi-db
```

It rewrites each cached file's path field in place by locating the same
relative path (type subdirectory + filename) under the new database
directory, leaving already-correct paths untouched and warning (with a
non-zero exit code) about any genome that genuinely isn't there. Use
`--dry-run` to preview the changes first. Afterwards, re-run
`scripts/set_targets.sh` as usual to regenerate a clean `targets.txt` from
the repaired cache.
