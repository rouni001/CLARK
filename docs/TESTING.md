# Testing and Coverage Notes

The current regression suite is intentionally lightweight. It is designed to
catch install, wrapper, and portability regressions quickly on both Linux and
macOS.

Run:

```sh
make test
```

To run a clean coverage-instrumented build and enforce the minimum coverage
gate:

```sh
make coverage
```

## Current Coverage

The suite currently contains 20 regression tests. It verifies:

- all required executables are created by the build
- `CLARK`, `CLARK-l`, and `CLARK-S` respond to `--version`
- modernized shell entrypoints parse successfully
- `classify_metagenome.sh` preserves paths with spaces
- gzipped inputs are decompressed and temporary files are cleaned up
- paired-end inputs are passed through correctly
- `--light` selects `CLARK-l`
- conflicting `--light` and `--spaced` options are rejected
- `getTargetsDef` emits expected target definitions for a tiny synthetic input
- `getAccssnTaxID` maps accession IDs and handles unmapped FASTA records
- `getfilesToTaxNodes` expands a tiny taxonomy lineage
- `exeSeq` splits a multi-FASTA file
- `dscriptMaker` emits deterministic download commands
- `getGammaDensity` and `getConfidenceDensity` summarize CLARK score columns
- `extractSeqs` extracts matching FASTQ records from a CLARK result file
- `getAbundance` summarizes a tiny assignment file
- `makeSummaryTables` writes summary tables for tiny abundance reports
- `getTargetSpecificKmersStat` counts target labels in a tiny database fixture
- the CLARK-l target filename regression stays fixed
- NCBI download URLs use the HTTPS host
- shell scripts avoid non-portable `readlink -f`

`make coverage` measures line coverage for compiled CLARK C++ sources with
`gcov`. It excludes system headers and merges the generated `build/default`,
`build/light`, and `build/spaced` copies back to their original `src/` paths so
shared lines are counted once.

GitHub Actions runs `make test` on Linux and macOS, and runs `make coverage` as
a dedicated coverage gate. The PR should not be merged if the coverage job drops
below the required minimum.

Current coverage from `make coverage`:

- covered executable C++ lines: 908
- total executable C++ lines: 5,657
- line coverage: 16.05%
- required minimum: 10.00%

## Coverage Limitations

The suite still does not measure shell line coverage; shell entrypoints are
validated through behavioral regression tests instead. More importantly, it does
not yet exercise the core k-mer database build and classification paths against
biological fixtures. The current C++ coverage is stronger for helper utilities
than for the core classifier algorithm.

## Recommended Next Tests

- small custom-database end-to-end tests for `CLARK`, `CLARK-l`, and `CLARK-S`
- golden FASTA/FASTQ fixtures with expected assignments
- tests for abundance estimation with taxonomy names, Krona output, and MPA
  output
- downloader tests using mocked NCBI manifests instead of live network calls
- converter tests for tiny contiguous-to-spaced k-mer databases
- continued replacement of fixed-size buffers and `sprintf` usage, followed by
  negative-path tests for long filenames and malformed input
