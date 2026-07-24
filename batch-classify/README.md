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
batch-classify/run_all.sh --db-dir /path/to/ncbi-db --types viruses,plasmid,plastid,fungi,human
```

Omit `--types` to use that same default set:

```sh
batch-classify/run_all.sh --db-dir /path/to/ncbi-db
```

For each type this:
1. Re-runs `scripts/set_targets.sh <db-dir> <type> <rank>`.
2. Samples a fresh `objects.fa` from that type's own targets via
   `scripts/make_sample.sh`.
3. Classifies it with `exe/<variant> -k <kmer> -T <targets> -D <DBD/> -O
   <objects.fa> -R <results> -n <threads>`.

Results land in `batch-classify/results/<type>/objects.fa` and
`.../results.csv` (gitignored -- generated output, not tracked).

## Options

| Flag                | Default    | Meaning                                  |
|----------------------|-----------|-------------------------------------------|
| `--db-dir <path>`   | (required) | database directory passed to `set_targets.sh` |
| `--types <list>`    | `viruses,plasmid,plastid,fungi,human` | comma- or space-separated database types |
| `--threads <n>`     | `8`       | thread count (CLARK's `-n`)                |
| `--sample-count <n>`| `20`      | reads sampled per type                     |
| `--sample-len <n>`  | `150`     | sampled read length                        |
| `--kmer <n>`        | `31`      | k-mer size (CLARK's `-k`)                  |
| `--profile`         | off       | set `CLARK_CUD_PROFILE=1` and `CLARK_CUD_PROFILE_FINE=1` for the CLARK run, printing its opt-in performance breakdown (build/load/match/write time, plus match/hits-update/classify sub-timings) |
| `--variant-exe <name>` | `CLARK` | executable to run (`CLARK`, `CLARK-l`, `CLARK-S`) |
| `--rank <flag>`     | `--species` | taxonomy rank flag for `set_targets.sh`  |
| `--results-dir <path>` | `batch-classify/results` | where per-type results land |

Run `batch-classify/run_all.sh --help` for the full usage message.

`CLARK_HOME` (env var) overrides the repository root if you're not running
from a checkout where this script's parent directory is the repo root.

Note: CLARK's default/full mode has a large fixed baseline memory
footprint regardless of database size (see main `README.md`'s Memory
Requirements section). If a type's build/classify step runs out of
memory, pass `--variant-exe CLARK-l` for a lower-RAM run.

Note: CLARK-l crashes with a "Bus error" if the thread count
(`--threads`) is greater than or equal to the number of objects being
classified (`--sample-count`). This is a pre-existing issue in CLARK-l
itself, not specific to this script. Keep `--sample-count` comfortably
above `--threads` (the defaults, 20 vs. 8, already do).
