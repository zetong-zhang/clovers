# Change Log
All notable changes to this project will be documented in this file.

## [1.2.3] - 2026-11-10

### Changed
- Metagenome runs are about 30% faster than v1.2.2: start-codon candidates are now pre-filtered with the init score, so TISA-CNN scores fewer ORFs (4.6 Mb fragmented input: 6.0 s -> 4.0 s on 8 threads). Predictions are unchanged.
- Single-genome prediction also skips ORFs with a very low init score and is about 20% faster than v1.2.2. Predictions are unchanged.
- GFF3 output now reports the gene start score as an additional `r_score` attribute.

### Fixed
- ORF coordinates no longer overflow when writing results for very short contigs.

## [1.2.2] - 2026-10-28

### Changed
- The scalar TIS-CNN backend, used when the build lacks AVX2/FMA, now scores start-codon candidates in parallel with OpenMP. On 8 threads the bundled `example.fa` drops from 6.7 s to 2.0 s, putting baseline builds within about 1.3x of the AVX2 backend. Predictions are unchanged.

## [1.2.1] - 2026-10-20

### Added
- `make install` target, honouring `PREFIX` and `DESTDIR`, and a `native` target for a static AVX-enabled build.

### Changed
- The build honours `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS` and `LDLIBS` from the environment.
- The default build is dynamically linked and targets baseline x86-64; `-static` and `-mavx -mfma -mavx2` moved behind the `STATIC=` and `NATIVE=` switches.
- OpenMP is switchable with `OPENMP=0`, and `OPENMP_FLAGS`/`OPENMP_LIBS` can be overridden for toolchains that keep it outside the default search path.
- `--version` writes to stdout.
- SVM scoring standardises the feature matrix in place instead of duplicating it, halving the temporary memory of that step.
- Reduced the SVM prediction cache from 200 to 50, lowering peak memory without affecting scores.
- `scripts/requirements.txt` records the Python version (3.9-3.11) and minimum versions of all dependencies.

### Fixed
- Corrected one codon's x-coordinate in the Z-curve feature table (`+1` -> `-1`) in both the C++ table (`include/Common.hpp`) and the Python reference table (`scripts/Zcurve.pyx`). Scores on genomes using that codon change slightly.
- The helper script directory is consistently named `scripts/`; it had been committed as `Scripts/`, which broke builds from the source archive on case-sensitive filesystems.

## [1.2.0] - 2026-09-25

### Added
- A pre-trained TIS model and a reinforcement-learning-trained CDS model for the metagenomic procedure.
- `-M/--model` option to force the TIS reviser (`hmm`, `cnn` or `none`). By default the reviser is chosen by the procedure and the genome size.

### Changed
- Option groups renamed: `CLOVERS` -> `CDS-Finder`, `TriTISA` -> `TIS-Reviser`.
- `--maxiter` no longer has the `-M` short form (now used by `--model`).

### Removed
- `-n/--bypass`; use `-M none` instead.

## [1.1.1] - 2026-09-22

### Changed
- ORF coordinates now use signed integers and the G+C counters were widened, removing signed/unsigned conversion and overflow hazards. Prediction results are unchanged from v1.1.0.

## [1.1.0] - 2026-05-20

### Changed
In this version, the code has been completely restructured. Unnecessary logic has been removed, and all functions have been encapsulated to minimize memory errors (e.g., array out-of-bounds and unreleased memory) and improve the running speed.

### Added
- Compilation support for 32-bit systems.
- Output the complete GenBank file containing the sequence.
- Specify the topology for each sequence in one assembly.
- Complete RNA genome processing capability.
- Replaced the `med` output format with standard BED6 (`bed`).

## [1.0.6] - 2026-04-22

### Changed
- Disable the TriTISA trainer in the metagenomic mode.
- Force the output of partial-5'-end ORFs with scores >= threshold for linear genomes. 

## [1.0.5] - 2026-04-03

### Added
- Added RNA sequence handling support for TriTISA and protein translation.
- Added translation table #15 (Blepharisma Nuclear Code).
- Incorporating the RBS score into the determination of the final score.
- Added version information display in the help menu.

### Changed
- Optimized the processing procedure for the phage genomes.
- Fixed some errors present in the Python scripts.
- Optimized the overlapping gene reporter.
- Fixed the calculation error of the Bayesian factor in TriTISA.
- For the protein output, replace non-AGCT chars with "X"s (unknown amino acids) instead of reporting an error.
- Fixed the issue where the external SVM model could not be loaded in certain cases.

## [1.0.4] - 2026-03-23

### Added
- Added auto-detection of translation table before training the model.
- Added more translation tables, and now there are 5 tables in total (1, 4, 11, 16, 25).

## [1.0.3] - 2026-03-20

### Added
- Added EMBL format support for input files.
- Added third-party terms for TriTISA to the LICENSE file.

### Changed
- Fixed a possible array out-of-bounds issue in the function `mm_train`.
- TIS revision has now been moved to take place before the CDS training.

## [1.0.2] - 2026-03-10

### Added
- Added the parameter `minolen` to filter overprinted genes with a minimum overprinted length between two ORFs;
- Added the length/topology of each sequence and date information to the output (GFF3, GBK and BED).

### Changed
- Optimized the exception handling of all the modules to improve performance;
- Fix the issue where GBK output cannot be parsed by Biopython.

## [1.0.1] - 2026-03-09

### Added
- Pretrained markov model for TriTISA to complete TIS revision on sequences where the samples are too few to initiate self-training.

### Changed
- Fixed the bug that overprinted genes cannot be detected at the boundaries of the circular genomes;
- Fixed the incorrect example format in the README file;
- Fixed the issue where svm.h could not be compiled on machines without AVX support.
- Fixed compilation issues on machines without OpenMP support.