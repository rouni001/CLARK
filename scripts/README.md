# CLARK Scripts

This directory contains the implementation scripts used to install CLARK,
prepare databases, classify reads, estimate abundance, and run maintenance
tasks.

The repository root keeps small compatibility launchers such as
`install.sh`, `set_targets.sh`, and `classify_metagenome.sh` so existing
documentation and user workflows continue to work. New maintenance scripts
should be added here instead of expanding the root directory.

All scripts resolve the CLARK repository root through `CLARK_HOME` when it is
set, or by using this directory's parent. This allows both of these forms:

```sh
./classify_metagenome.sh ...
scripts/classify_metagenome.sh ...
```
