# CLARK: Fast and Accurate Sequence Classifier

## Authors

- **Rachid Ounit** (1)
- **Steve Wanamaker** (2)
- **Timothy J. Close** (2)
- **Stefano Lonardi** (1)

1. Computer Science and Engineering Department, University of California, Riverside, CA 92521
2. Botany and Plant Sciences Department, University of California, Riverside, CA 92521

## Publications

Please cite the following peer-reviewed publications when using CLARK:

1. Ounit, R., Wanamaker, S., Close, T. J., & Lonardi, S. (2015). CLARK: fast and accurate classification of metagenomic and genomic sequences using discriminative k-mers. BMC Genomics, 16(1), 236.
2. Ounit, R., & Lonardi, S. (2016). Higher classification sensitivity of short metagenomic reads with CLARK-S. Bioinformatics, 32(24), 3823-3825.

## About

CLARK (CLAssifier based on Reduced K-mers) addresses the computational challenges of DNA sequence classification in molecular biology, genomics, metagenomics, and genetics. It uses discriminative k-mers for supervised sequence classification, providing a confidence score for each assignment. CLARK supports three classifiers:

- **CLARK (default):** Designed for powerful workstations and large databases, using exact k-mer matching.
- **CLARK-l:** A memory-efficient version for small metagenomes and limited-memory workstations, also using exact k-mer matching.
- **CLARK-S:** Uses spaced k-mers for higher sensitivity, requiring more RAM.

## General Notes

- CLARK was introduced at the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (ACM-BCB) in 2014.
- CLARK-S was presented at the 14th Workshop on Algorithms in Bioinformatics in 2015.

## License

CLARK is distributed under the GNU General Public License (GPL) v3. It is free software, available without any warranty. See the GNU General Public License for more details: [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)

## Releases

### Version 1.3.0 (May 16, 2024)
- Fixed error when processing bacterial/viral genomes during database construction.
- Update reference for human genome.
- Extend list of references to download for plasmid/plastids.
- Code improvements and refactoring included.

### Previous Versions
- **Version 1.2.6.1 (May 11, 2019):** Various improvements and bug fixes.
- **Version 1.2.6 (February 21, 2019):** New features and performance enhancements.
- **Version 1.2.5.1 (September 20, 2018):** Added script for k-mer distribution and support for long reads.
- **Version 1.2.5 (February 13, 2018):** Extended fungi database and added plasmid, plastid, and protozoa databases.
- **Version 1.2.4 (November 25, 2017):** Introduced fungi database for hospital-acquired infections.
- **Version 1.2.3.2 (August 21, 2017):** Improved summary script and support for paired-end reads in Fasta format.
- **Version 1.2.3.1 (February 27, 2017):** Added summary script and updated contact information.
- **Version 1.2.3 (April 7, 2016):** Finalized CLARK-S and added scripts for custom databases.
- **Version 1.2.2-b (December 11, 2015):** Improved speed and RAM usage for full and spaced modes.
- **Version 1.2.1-b (June 15, 2015):** Introduced discriminative spaced k-mers.
- **Version 1.1.3 (June 3, 2015):** Extended features and simplified output.
- **Version 1.1.2 (April 22, 2015):** Added scripts for abundance estimation and gzip support.
- **Version 1.1.1 (April 15, 2015):** Improved efficiency for loading databases.
- **Version 1.1 (February 20, 2015):** Improved storage efficiency and multithreaded database loading.
- **Version 1.0 (September 1, 2014):** Initial release.

## Compatibility Between Versions

- Databases created with version 1.0 are not compatible with v1.1 or newer. Rebuild the database using v1.1 or newer.

## Software & System Requirements

### External Tools & C++ Compiler Version
- CLARK requires a 64-bit OS (Linux or Mac) and GNU GCC 4.4 or higher.
- Multithreading operations use the OpenMP libraries.

### Memory Requirements
- CLARK can be RAM-intensive. Detailed memory usage for different versions and configurations is provided below:

  - **CLARK (default):** 58 GB RAM for loading the bacterial database, 156 GB for building it.
  - **CLARK-l:** Designed for 4 GB RAM laptops, suitable for small metagenomes.
  - **CLARK-S:** Up to 101 GB RAM for spaced k-mer classification.

## Installation

1. Download the latest version from the CLARK webpage: [http://clark.cs.ucr.edu](http://clark.cs.ucr.edu)
2. Uncompress the tar.gz file: `tar -xvf CLARKV1.3.0.tar.gz`
3. Navigate to the CLARK directory and run the installation script: `./install.sh`

## Usage

### Scripts

- **set_targets.sh** and **classify_metagenome.sh**: Define your database and classify metagenomes.
- **estimate_abundance.sh**: Compute abundance estimation (count/proportion of objects assigned to targets).
- Additional scripts are provided for various tasks such as building spaced k-mer databases, resetting custom databases, and updating taxonomy data.

### Basic Example

1. Define targets: `./set_targets.sh <DIR_DB/> bacteria`
2. Classify metagenome: `./classify_metagenome.sh -O ./sample.fa -R ./result`

## How to Choose the K-mer Length

- For high sensitivity, use k = 20 or 21.
- For a balance between speed, precision, and RAM usage, use k = 31.

## Manual & Options

A detailed description of all command-line options and parameters is provided below.

### Definitions of Parameters

- **-k \<kmerSize\>**: K-mer length (default: 31, 27 for CLARK-l, 31 for CLARK-S with max mismatches set to 9).
- **-T \<fileTargets\>**: Targets definition file (mandatory).
- **-t \<minFreqTarget\>**: Minimum k-mer frequency in targets (default: 0).
- **-D \<directoryDB\>**: Directory for the database (mandatory).
- **-O \<fileObjects\>**: File containing objects to be classified (mandatory).
- **-P \<file1\> \<file2\>**: Paired-end fastq files.
- **-o \<minFreqObject\>**: Minimum k-mer frequency in objects (spectrum mode only).
- **-R \<fileResults\>**: File to store results (mandatory).
- **-m \<mode\>**: Mode of execution (0: full, 1: default, 2: express, 3: spectrum).
- **-n \<numberofthreads\>**: Number of threads (default: 1).
- **-g \<iterations\>**: Gap for non-overlapping k-mers (CLARK-l only, default: 4).
- **-s \<factor\>**: Sampling factor for discriminative k-mers (default: 2).
- **--long**: For very long sequences (e.g., long contigs, Nanopore/Pacbio reads).
- **--tsk**: Detailed creation of the database (no longer supported).
- **--ldm**: Load database file by memory-mapped file.
- **--kso**: Preliminary k-spectrum analysis (spectrum mode only).
- **--extended**: Extended output (full mode only).

### Examples

- Default mode with 31-mers: `./CLARK -k 31 -T ./targets_addresses.txt -D ./DBD/ -O ./objects.fa -R ./results`
- Full mode with 20-mers: `./CLARK -k 20 -T ./targets_addresses.txt -D ./DBD/ -O ./objects.fa -R ./results -m 0`
- Express mode with 12 threads: `./CLARK -k 20 -T ./targets_addresses.txt -D ./DBD/ -O ./objects.fa -R ./results -m 2 -n 12`
- Paired-end reads: `./CLARK -k 20 -T ./targets_addresses.txt -D ./DBD/ -P ./sample1.fastq ./sample2.fastq -R paired.results`

## Classification of Metagenomic Samples

Scripts provided for metagenomic classification:

### Setting Targets

1. Create a directory to store reference sequences (e.g., `<DIR_DB/>`).
2. Define targets:
   - Only bacteria: `./set_targets.sh <DIR_DB/> bacteria`
   - Bacteria, viruses, and human: `./set_targets.sh <DIR_DB/> bacteria viruses human`
   - Bacteria and custom: `./set_targets.sh <DIR_DB/> bacteria

 custom`

### Running the Classification

Example commands:

- Classify metagenome: `./classify_metagenome.sh -O ./sample.fa -R ./result`
- Use 20-mers: `./classify_metagenome.sh -O ./sample.fa -R ./result -k 20`
- Full mode: `./classify_metagenome.sh -O ./sample.fa -R ./result -m 0`
- Multiple sample files: `./classify_metagenome.sh -O ./samples.txt -R ./samples.txt -m 0`
- Paired-end reads: `./classify_metagenome.sh -O ./samples.R.txt ./samples.L.txt -R ./samples.R.txt -m 0`
- 8 threads: `./classify_metagenome.sh -O ./sample.fa -R ./result -m 2 -n 8`
- Gzipped objects file: `./classify_metagenome.sh -O ./sample.fa.gz -R ./result -m 0 -n 8 --gzipped`
- Use CLARK-l: `./classify_metagenome.sh -P ./sample1.fastq ./sample2.fastq -R ./result --light`
- Use CLARK-S: `./classify_metagenome.sh -P ./sample1.fastq ./sample2.fastq -R ./result --spaced`
- CLARK-S with full mode and 8 threads: `./classify_metagenome.sh -O ./sample.fa -R ./result -m 0 -n 8 --spaced`
- CLARK-S with express mode on gzipped file: `./classify_metagenome.sh -O ./sample.fa.gz -R ./result -m 2 -n 8 --spaced`
- Lower RAM usage for CLARK-S: `./classify_metagenome.sh -O ./sample.fa -R ./result --spaced -s 2`

### Abundance Estimation

Example commands:

- Basic usage: `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/>`
- High confidence assignments: `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/> --highconfidence`
- Filter by confidence score (e.g., 0.8): `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/> -c 0.80`
- Filter by gamma score (e.g., 0.03): `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/> -g 0.03`
- Output in MetaPhlAn format: `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/> -g 0.03 --mpa`
- Filter by abundance (e.g., >2%): `./estimate_abundance.sh -F ./result.csv -D <DIR_DB/> -a 2`

## Results Format

- **Full mode:** `<Object_ID>,<hit count in target 1>,...,<hit count in target N>,<Length of object>,<Gamma>,<first assignment>,<hit count of first>,<second assignment>,<hit count of second>,<confidence score>`
- **Default/Express mode:** `<Object_ID>,<Length of object>,<1st_assignment>`

## COMMUNITY OF CLARK USERS

For discussions, suggestions, or questions please join the Google groups of the CLARK users at:
 https://groups.google.com/g/clarkusers

---

This README provides a comprehensive guide to installing, using, and understanding CLARK. It ensures clarity and ease of use for all users. For more detailed information and examples, please refer to the full README file included with the CLARK distribution.