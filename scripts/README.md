# CLARK Scripts

This directory contains the implementation scripts used to install CLARK,
prepare databases, classify reads, estimate abundance, and run maintenance
tasks.

The repository root intentionally does not keep duplicate script launchers. Use
these scripts directly from this directory so the top-level project stays
focused on source, docs, tests, and build metadata.

All scripts resolve the CLARK repository root through `CLARK_HOME` when it is
set, or by using this directory's parent. For example:

```sh
scripts/classify_metagenome.sh ...
```
