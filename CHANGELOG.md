# Change Log
All notable changes to this project will be documented in this file.

## [1.0.1] - 2026-03-09

### Added
- Pretrained markov model for TriTISA+ to complete TIS revision on sequences where the samples are too few to initiate self-training.

### Changed
- Fixed the bug that overprinted genes cannot be detected at the boundaries of the circular genomes;
- Fixed the incorrect example format in the README file;
- Fixed the issue where svm.h could not be compiled on machines without AVX support.
- Fixed compilation issues on machines without OpenMP support.