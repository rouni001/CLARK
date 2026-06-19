# Testing and Coverage Notes

The current regression suite is intentionally lightweight. It is designed to
catch install, wrapper, and portability regressions quickly on both Linux and
macOS.

Run:

```sh
make test
```

## Current Coverage

The suite currently verifies:

- all required executables are created by the build
- `CLARK`, `CLARK-l`, and `CLARK-S` respond to `--version`
- modernized shell entrypoints parse successfully
- `classify_metagenome.sh` preserves paths with spaces
- gzipped inputs are decompressed and temporary files are cleaned up
- paired-end inputs are passed through correctly
- `--light` selects `CLARK-l`
- conflicting `--light` and `--spaced` options are rejected
- `getTargetsDef` emits expected target definitions for a tiny synthetic input
- the CLARK-l target filename regression stays fixed
- NCBI download URLs use the HTTPS host
- shell scripts avoid non-portable `readlink -f`

## Coverage Limitations

The suite does not yet measure line or branch coverage. More importantly, it
does not yet exercise the core k-mer database build and classification paths
against biological fixtures. The current C++ coverage is smoke-level: CLI
startup/version paths and one helper utility are tested, but the classifier
algorithm is not validated end-to-end.

## Recommended Next Tests

- small custom-database end-to-end tests for `CLARK`, `CLARK-l`, and `CLARK-S`
- golden FASTA/FASTQ fixtures with expected assignments
- abundance-estimation tests with known CSV inputs
- downloader tests using mocked NCBI manifests instead of live network calls
- optional instrumented coverage target once compiler/tooling support is
  standardized across Linux and macOS
