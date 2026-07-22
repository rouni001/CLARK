# dev-test: offline pipeline smoke test

This folder is a small, self-contained example that exercises the whole
CLARK_CuD pipeline (`set_targets.sh` -> `classify_metagenome.sh`) without
downloading anything from NCBI. It exists purely to let you verify the
scripts/executables work in your environment before pointing them at a
real, large database.

It is **not** biological data — `db/Bacteria/*.fna` are randomly
generated sequences standing in for two "organisms", and `sample.fa` is
a few reads taken directly from one of them.

## Layout

- `db/Bacteria/` - two synthetic genomes (taxids `900001` and `900002`)
- `db/taxonomy/` - a minimal `nodes.dmp`/`merged.dmp` covering just those
  two taxids (no real NCBI taxonomy download needed)
- `setup.sh` - (re)generates the path-dependent metadata files
  (`db/.bacteria`, `db/.bacteria.provenance.tsv`, `db/.taxondata`) so
  `make_metadata.sh` sees the sequences as already downloaded, wherever
  you cloned this repo
- `sample.fa` - 3 reads taken from `GCF_900000001.1_DemoOrgA`, expected
  to classify as taxid `900001`

## Usage

Run these from the repository root (after `make all`):

```sh
dev-test/setup.sh
CLARK_HOME="$(pwd)" scripts/set_targets.sh dev-test/db bacteria --species
CLARK_HOME="$(pwd)" scripts/classify_metagenome.sh -O dev-test/sample.fa -R dev-test/result --light
cat dev-test/result.csv
```

Expected output: all 3 reads in `dev-test/result.csv` assigned to taxid
`900001`.

### Why `--light`?

CLARK's default/full mode allocates a fixed ~1.6 billion-bucket hash
table (`HTSIZE` in `src/parameters.hh`) regardless of how small the
target database is, so even this 2-genome demo can throw
`std::bad_alloc` on a machine with less RAM than CLARK's baseline
footprint (the main README documents ~58 GB to load / ~156 GB to build
a real bacterial database). `--light` (CLARK-l) is designed for
constrained RAM and is enough to run this demo on an ordinary machine.
If you have enough RAM, you can drop `--light` and use the default mode
instead.

## Cleaning up

`setup.sh`'s generated files and the build output (`db/targets.txt`,
`db/bacteria_*/`, `result*`) are gitignored (`dev-test/.gitignore`) —
re-running `setup.sh` and the pipeline again is always safe.
