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

## `run_benchmark.sh`: a citable benchmark input, with accuracy scoring

`run_all.sh` samples uniformly-random substrings via `scripts/make_sample.sh`
-- fine for a smoke test, but not something a paper can cite as a benchmark.
`run_benchmark.sh` instead uses `scripts/make_benchmark_reads.sh`, which
simulates reads matching the *published characteristics* of datasets other
metagenomic-classifier papers have used:

| Profile  | Read length | Source genomes         | Error rate | Modeled after |
|----------|--------------|-------------------------|------------|----------------|
| `hiseq`  | 92 bp        | 10, equal proportion     | baseline   | Kraken's "HiSeq" set |
| `miseq`  | 156 bp       | 10, equal proportion     | baseline   | Kraken's "MiSeq" set |
| `simba5` | 100 bp       | broad draw, all available | 5x baseline | Kraken's "simBA-5" set |
| `simhc`  | 800 bp       | up to 113, power-law abundance | elevated | FAMeS's "simHC" set |
| `custom` | you choose   | you choose               | you choose | -- |

CLARK's own BMC Genomics (2015) paper benchmarked against the Kraken
project's "HiSeq"/"MiSeq"/"simBA-5" read sets and the FAMeS project's
"simHC" mock community; the ISCA 2021 "Sieve" paper (hardware acceleration
for k-mer classification) reuses the same lineage. The *original* files
are hosted on decade-old project pages (`ccb.jhu.edu`, `fames.jgi-psf.org`)
that may no longer serve the exact original bytes -- `run_benchmark.sh`
simulates reads matching each profile's spec from whatever genomes
`scripts/set_targets.sh` already downloaded, rather than re-fetching those
files. Treat this as a documented approximation in anything you write up,
not a claim of using the literal original dataset.

```sh
batch-classify/run_benchmark.sh /path/to/ncbi-db hiseq viruses plasmid plastid fungi human
```

For each type this samples benchmark reads (with ground truth), classifies
them, and scores sensitivity/precision (overall and per-taxid) with
`scripts/eval_accuracy.py`. Results land in `batch-classify/results/<type>/`:
`objects.fa` (reads), `objects.fa.truth.tsv` (ground truth), `results.csv`
(CLARK's output), and `accuracy.txt` (sensitivity/precision).

Additional environment overrides (on top of the ones above):

| Variable                     | Meaning                                         |
|-------------------------------|-------------------------------------------------|
| `CLARK_BENCHMARK_COUNT`       | override the profile's default read count        |
| `CLARK_BENCHMARK_LEN`         | override the profile's default read length        |
| `CLARK_BENCHMARK_ERROR_RATE`  | override the profile's default error rate         |
| `CLARK_BENCHMARK_NGENOMES`    | override the profile's default source-genome count |
| `CLARK_BENCHMARK_SEED`        | random seed for reproducible sampling             |

`-p custom` requires count/length/error-rate (`CLARK_BENCHMARK_COUNT` /
`CLARK_BENCHMARK_LEN` / `CLARK_BENCHMARK_ERROR_RATE`) to be set explicitly.
