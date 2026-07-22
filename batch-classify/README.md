# batch-classify

Runs a full classification pass (sample -> CLARK) against each of several
already-downloaded databases that share one database directory (e.g. one
`ncbi-db/` populated by running `scripts/set_targets.sh` once per type:
`viruses`, `plasmid`, `plastid`, `fungi`, `human`, ...).

## Why re-run `set_targets.sh` per type?

`set_targets.sh <db-dir> <type>` writes `<db-dir>/targets.txt` and points
`.settings` at `<db-dir>/<type>_0/`. If you ran it once per type against the
same `<db-dir>`, only the **last** type's `targets.txt` survives -- each
type's own `<type>_0/` database folder is untouched, but the shared
`targets.txt` gets overwritten every time.

`run_all.sh` handles this by re-running `set_targets.sh` for each type right
before classifying it. Since the sequences and taxonomy are already
downloaded, this makes no network calls -- it only regenerates
`targets.txt`/`.settings` for that type and points `-D` at its existing
`<type>_0/` folder (already-built k-mer databases are reused, not rebuilt).

## Usage

```sh
batch-classify/run_all.sh /path/to/ncbi-db viruses plasmid plastid fungi human
```

Omit the type list to use that same default set:

```sh
batch-classify/run_all.sh /path/to/ncbi-db
```

For each type this:
1. Re-runs `scripts/set_targets.sh <db-dir> <type> --species`.
2. Samples a fresh `objects.fa` from that type's own targets via
   `scripts/make_sample.sh`.
3. Classifies it with `exe/CLARK -k <kmer> -T <targets> -D <DBD/> -O
   <objects.fa> -R <results> -n <threads>`.

Results land in `batch-classify/results/<type>/objects.fa` and
`.../results.csv` (gitignored -- generated output, not tracked).

## Environment overrides

| Variable              | Default    | Meaning                                  |
|------------------------|-----------|-------------------------------------------|
| `CLARK_HOME`           | repo root | where `scripts/` and `exe/` live           |
| `CLARK_VARIANT_EXE`    | `CLARK`   | executable to run (`CLARK`, `CLARK-l`, `CLARK-S`) |
| `CLARK_KMER`           | `31`      | k-mer size                                 |
| `CLARK_THREADS`        | `8`       | thread count                               |
| `CLARK_RANK`           | `--species` | taxonomy rank flag for `set_targets.sh`  |
| `CLARK_SAMPLE_COUNT`   | `20`      | reads sampled per type                     |
| `CLARK_SAMPLE_LEN`     | `150`     | sampled read length                        |
| `CLARK_BATCH_RESULTS_DIR` | `batch-classify/results` | where per-type results land |

Note: CLARK's default/full mode has a large fixed baseline memory
footprint regardless of database size (see main `README.md`'s Memory
Requirements section). If a type's build/classify step runs out of
memory, set `CLARK_VARIANT_EXE=CLARK-l` for a lower-RAM run.

Note: CLARK-l crashes with a "Bus error" if the thread count
(`CLARK_THREADS`) is greater than or equal to the number of objects
being classified (`CLARK_SAMPLE_COUNT`). This is a pre-existing issue
in CLARK-l itself, not specific to this script. Keep `CLARK_SAMPLE_COUNT`
comfortably above `CLARK_THREADS` (the defaults, 20 vs. 8, already do).
