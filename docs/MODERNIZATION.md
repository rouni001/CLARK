# Modernization Notes

This version focuses on making CLARK easier to install, safer to run, and more
reliable for researchers who are not bioinformatics specialists.

## Critical Issues Addressed

1. **Unreliable installation and hidden build behavior**
   - Added a `Makefile` used by both `make` and `scripts/install.sh`.
   - Added explicit OpenMP detection and clear single-threaded fallback.
   - Generated binaries now live under ignored `exe/` output instead of being
     treated as source artifacts.

2. **CLARK-l target filename correctness bug**
   - Fixed CLARK-l target hash-table path generation so the regular target-label
     loop uses `m_labels[t]`, not `m_labels_c[t]`.

3. **Unsafe shell scripts for common research paths**
   - Rewrote `scripts/classify_metagenome.sh`, `scripts/set_targets.sh`, and
     `scripts/make_metadata.sh` with portable script-directory resolution,
     quoted arguments, clear errors, and safer gzipped input handling.
   - Consolidated shell scripts under `scripts/` so the repository root does
     not keep duplicate launchers.

4. **Fragile database/taxonomy download setup**
   - Updated taxonomy download to use HTTPS NCBI URLs, `curl`/`wget` fallback,
     and safer target-directory handling.
   - Replaced misspelled/legacy NCBI FTP hosts in the RefSeq downloader with
     `https://ftp.ncbi.nlm.nih.gov`.

5. **No automated regression safety net**
   - Added focused tests for script quoting, gzipped input behavior, source-level
     CLARK-l regression, URL sanity, and version execution.
   - Added a CI workflow that builds and runs tests on Linux and macOS.

## Remaining Important Work

- Replace the large shell RefSeq downloader with a manifest-driven downloader
  that records accessions, checksums, dates, and taxonomy release metadata.
- Add small biological golden datasets with expected assignments for CLARK,
  CLARK-l, and CLARK-S.
- Add package recipes, such as Conda/Bioconda and container images.
- Continue replacing fixed-size C buffers and legacy parsing with typed C++
  interfaces.
