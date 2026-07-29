# batch-classify

Runs a full classification pass (sample -> classify) against each of several
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
batch-classify/run_all.sh -d /path/to/ncbi-db -t viruses,plasmid,plastid,fungi,human
```

Omit `-t` to use that same default set:

```sh
batch-classify/run_all.sh -d /path/to/ncbi-db
```

For each type this:
1. Re-runs `scripts/set_targets.sh <db-dir> <type> <rank>`.
2. Samples a fresh `objects.fa` from that type's own targets via
   `scripts/make_sample.sh`.
3. Classifies it with `-k <kmer> -T <targets> -D <DBD/> -O <objects.fa>
   -R <results> -n <threads>`.

Results land in `batch-classify/results/<type>/objects.fa` and
`.../results.csv` (gitignored -- generated output, not tracked).

## Options

| Flag         | Default    | Meaning                                  |
|--------------|-----------|--------------------------------------------|
| `-d <path>`  | (required) | database directory passed to `set_targets.sh` |
| `-t <list>`  | `viruses,plasmid,plastid,fungi,human` | comma- or space-separated database types |
| `-n <n>`     | `8`       | thread count (`-n`)                        |
| `-c <n>`     | `20`      | reads sampled per type                     |
| `-l <n>`     | `150`     | sampled read length                        |
| `-k <n>`     | `31`      | k-mer size (`-k`)                          |
| `-p`         | off       | set `CLARK_CUD_PROFILE=1` and `CLARK_CUD_PROFILE_FINE=1` for the run, printing CLARK's opt-in performance breakdown (build/load/match/write time, plus match/hits-update/classify sub-timings). CPU only; ignored with `-g`. |
| `-e`         | off       | measure CPU package energy (RAPL) for the run via `scripts/rapl_energy.py`, printing an `ENERGY_PROFILE` line. |
| `-g`         | off       | run on GPU with cuCLARK instead of CPU CLARK. Overrides `-x`. |
| `-x <name>`  | `CLARK`   | CPU executable variant (`CLARK`, `CLARK-l`, `CLARK-S`). Ignored if `-g` is set. |
| `-r <flag>`  | `--species` | taxonomy rank flag for `set_targets.sh`  |
| `-o <path>`  | `batch-classify/results` | where per-type results land |

Run `batch-classify/run_all.sh -h` for the full usage message.

`CLARK_HOME` (env var) overrides the repository root if you're not running
from a checkout where this script's parent directory is the repo root.

### Energy (`-e`, RAPL)

`-e` wraps the classify step with `scripts/rapl_energy.py`, which reads
Intel RAPL's package-energy counters (`/sys/class/powercap/intel-rapl:*`,
summing one `package-*` domain per CPU socket) before and after the run
and prints:

```
ENERGY_PROFILE unit=joules pkg_joules=12.345678 pkg_joules_scaled=8.641975 scale=0.70 elapsed_s=1.234567 domains=1
```

This is the same MSR data Intel's PCM tool reports, just read directly
from sysfs instead of via PCM/perf (no root or `msr` kernel module
needed, as long as `energy_uj` is readable). RAPL package energy is
**not** whole-system wall power -- it excludes the PSU, fans, disks, and
(on most server/Xeon platforms) DRAM. There's no universal correction
factor for that gap, so alongside the raw `pkg_joules` figure,
`pkg_joules_scaled` reports the same value scaled down by a fixed 30%
(`scale=0.70`) as a rough, clearly-labeled second estimate -- treat both
as component-level numbers, not a measurement of true wall power.

If RAPL is unavailable or unreadable, `rapl_energy.py` prints a warning
and still runs the classifier normally, just without an `ENERGY_PROFILE`
line -- `-e` never blocks a run. `RAPL_SYSFS_DIR` (env var) overrides the
`/sys/class/powercap` path if yours is mounted elsewhere.

### GPU (cuCLARK)

`-g` swaps the classifier from the CPU `exe/<variant>` binary to cuCLARK.
Everything else about the run (targets, database directory, objects,
results, `-k`/`-n`) is passed to it exactly as it would be to CLARK, since
cuCLARK accepts the same CLI.

By default `-g` looks for a binary named `cuCLARK` sitting right next to
`run_all.sh`, i.e. `batch-classify/cuCLARK` -- drop your executable there
under that name and `-g` works with no further setup. If yours lives
elsewhere or is named differently, either set the `CUCLARK_EXE`
environment variable to its path, or edit the `CUCLARK_EXE` line at the
top of `run_all.sh`:

```sh
CUCLARK_EXE=/path/to/your-cuCLARK-binary batch-classify/run_all.sh -d /path/to/ncbi-db -g
```

`-p`'s CuD profiling env vars are CLARK-specific and are not set when
`-g` is used.

Note: CLARK's default/full mode has a large fixed baseline memory
footprint regardless of database size (see main `README.md`'s Memory
Requirements section). If a type's build/classify step runs out of
memory, pass `-x CLARK-l` for a lower-RAM run.

Note: CLARK-l crashes with a "Bus error" if the thread count (`-n`) is
greater than or equal to the number of objects being classified (`-c`).
This is a pre-existing issue in CLARK-l itself, not specific to this
script. Keep `-c` comfortably above `-n` (the defaults, 20 vs. 8, already
do).
