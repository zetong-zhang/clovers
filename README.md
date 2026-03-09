# CLOVERS
<br/>

![LOGO](logo.png)

<br/>

![C++](https://img.shields.io/badge/language-C%2B%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-brightgreen)
![Genomics](https://img.shields.io/badge/domain-genomics-purple)
![AVX](https://img.shields.io/badge/optimized-AVX%2FFMA-red)
![Static](https://img.shields.io/badge/build-static-lightgrey)
![Dependencies](https://img.shields.io/badge/dependencies-minimal-green)
[![GPLv3 License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![x86_64](https://img.shields.io/badge/arch-x86__64-green)
![GitHub release](https://img.shields.io/github/v/release/zetong-zhang/CLOVERS)
![GitHub downloads](https://img.shields.io/github/downloads/zetong-zhang/CLOVERS/total)

Ab initio prediction of overprinted genes using the Z-curve method

## Contents
- **[Overview](#overview)** - Project introduction and features
- **[Setup](#setup)** - Installation guide for different operating systems
- **[Usage](#usage)** - Command-line interface and options
- **[Examples](#examples)** - Usage examples with command-line options
- **[Input](#input)** - Description of the input files and formats
- **[Output](#output)** - Description of the output files and formats
- **[Structure](#structure)** - Project file structure and organization
- **[Server](#server)** - Free available server for running tasks online
- **[Citation](#citation)** - Citation guide for academic use
- **[Contact](#contact)** - Contact information for questions or support
- **[License](#license)** - GNU General Public License v3.0

## Overview
CLOVERS is an ab initio gene finder that utilizes the Z-curve method for robust gene prediction on prokaryotic and phage genomes. Its novel technical approach gives it an advantage in detecting more smaller overprinted genes, which are often missed by other tools.

## Setup
### Windows/Linux (x86_64)
Download the latest precompiled binary file from the [release page](https://github.com/zetong-zhang/CLOVERS/releases), and decompress it to the directory of your choice.
### Other Operating Systems
Download and compile the source code yourself.

```bash
git clone https://github.com/zetong-zhang/clovers.git
cd clovers
make  # use mingw32-make on windows (MinGW64)
```
*Note:*  
(1) If your system doesn't have zlib installed, or you just don't want to use zlib, just remove the `-DZLIB` and `-lz` option in `CXXFLAGS`;   
(2) If you don't want to use AVX and FMA function, just remove the `-mavx` and `-mfma` option in `CXXFLAGS`;  
(3) If you don't want to use OpenMP function, just remove the `-fopenmp` option in `CXXFLAGS`.

## Usage
### Quick Start
We recommend configure the environment variable `PATH` to include the directory of the executable binary file, such that you can run `clovers` directly in the terminal.
```bash
clovers -i example.fasta -o example.gff -c -f gff
```
### Options

```bash
clovers.exe [OPTION...]
```

#### General Options
* `-h, --help`  
Print help menu and exit. If no option is specified, print the help menu and wait.  

* `-q, --quiet`  
Run quietly with no stderr output.  

* `-T, --threads`  
Number of threads to use. Please set as a positive integer. (default: all)

#### Input/Output Options
* `-i, --input`   
Specify FASTA/Genbank input file or their compressed versions (gzip). (default: stdin)  

* `-o, --output`  
Specify output file. (default: stdout)

* `-f, --format`  
Select output format (gff, gbk, med). (default: gff)

* `-a, --faa`  
Write protein translations of genes to the selected file.

* `-d, --fna`  
Write nucleotide sequences of genes to the selected file.

#### CLOVERS Options  
* `-g, --table`  
Specify a genetic codon table to use. (default: 11)

* `-l, --minlen`  
Specify the mininum length (nt) of ORFs. (default: 90)

* `-c, --circ`  
Treat all the sequences's topology as circular.

* `-p, --proc`  
Select prediction procedure (single or meta).

* `-s, --thres`  
Specify putative gene probability score threshold. (default: 0.5)

* `-t, --train`  
Write (if none exists) or use the specified training file.

#### TriTISA+ Options  
* `-n, --bypass`  
Bypass TriTISA+ and output longest ORFs.

* `-M,--maxiter`  
Max iteration times for TIS revision. (default: 20)

* `-I,--tis`  
Write (if none exists) or use the specified TIS model file.

#### GOP-Reporter options:  
* `-O, --overprint`  
Write overprinted genes to the selected file (GFF3 format).

* `-A, --amino`  
Write protein translations of overprinted genes to the selected file.

* `-D, --nucl`  
Write nucleotide sequences of overprinted genes to the selected file.

* `-M, --ratio`  
Specify minimum overlap ratio for overprinted genes. (default: 0.6)

## Examples
### Basic Usage
```bash
# Predict genes in a FASTA file and output to GFF3
clovers -i example.fasta -o example.gff

# Predict genes in a circular genome
clovers -i circular.fasta -o output.gff -c

# Use a custom translation table (e.g., table 4 for Mycoplasma)
clovers -i input.fasta -o output.gff -g 4
```

### Advanced Usage
```bash
# Generate protein and nucleotide sequences
clovers -i input.fasta -o output.gff -a proteins.faa -d genes.fna

# Detect overprinted genes and their sequences
clovers -i input.fasta -o output.gff -O overprint.gff -A overprint.faa -D overprint.fna

# Use multiple threads for faster processing
clovers -i large.fasta -o output.gff -T 8

# Predict in metagenomes and use longest ORFs
clovers -i input.fasta -o output.gff -p meta -n

# Generate training file and TIS model file
clovers -i input.fasta -o output.gff -t training.dat -I tis_model.bin
```

## Input
CLOVERS supports FASTA format input (stdin or file) or their compressed versions (gzip) with content like follow:
```
>NC_000913.3
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC
. . .
```
Or minimal GenBank format input (stdout or file) or their compressed versions (gzip) with content like follow:
```
LOCUS       NC_000913.3
DEFINITION  NC_000913.3
FEATURES             Location/Qualifiers
ORIGIN   
        1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
       61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
      121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
          . . .
//
```

## Output
- **GFF/GenBank/MED files**: Primary annotation results containing gene locations, scores and other attributes

    **GFF File Example:**
    ```
    ##gff-version 3
    NC_000913.3	CLOVERS_v1.0.0	CDS	337	2799	0.974	+	0	ID=orf000001
    NC_000913.3	CLOVERS_v1.0.0	CDS	2801	3733	0.968	+	0	ID=orf000002
    . . .
    ```
    **GenBank File Example:**
    ```
    LOCUS       NC_000913.3
    DEFINITION  NC_000913.3
    FEATURES             Location/Qualifiers
        CDS             337..2799
                        /note="version=CLOVERS_v1.0.0;ID=orf000001;score=0.974"
        CDS             2801..3733
                        /note="version=CLOVERS_v1.0.0;ID=orf000002;score=0.968"
    . . .
    ORIGIN
    //
    ```
    **MED File Example:**
    ```
    # NC_000913.3
    337 2799	+
    2801 3733	+
    . . .
    ```
- **Protein sequence files (.faa)**: Translated amino acid sequences of predicted genes
    ```
    >NC_000913.3:2801..3733(+) score=0.968
    MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQI
    PVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDI
    ISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLP
    GFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*
    . . .
    ```
- **Nucleotide sequence files (.fna)**: DNA sequences of predicted genes
    ```
    >NC_000913.3:2801..3733(+) score=0.968
    ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCTGT
    TGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCGGCAGAGACATTCAGTCTCAACAACCTCGGACGCTTTGCCG
    ATAAGCTGCCGTCAGAACCACGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCAGGAACTGGGTAAGCAAATT
    CCAGTGGCGATGACCCTGGAAAAGAATATGCCGATCGGTTCGGGCTTAGGCTCCAGTGCCTGTTCGGTGGTCGCGGCGCT
    GATGGCGATGAATGAACACTGCGGCAAGCCGCTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTGGAAGGCCGTA
    TCTCCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTCGGTGGTATGCAGTTGATGATCGAAGAAAACGACATC
    ATCAGCCAGCAAGTGCCAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGGATTAAAGTCTCGACGGCAGAAGC
    CAGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATCTGGCAGGCTTCATTCACGCCT
    GCTATTCCCGTCAGCCTGAGCTTGCCGCGAAGCTGATGAAAGATGTTATCGCTGAACCCTACCGTGAACGGTTACTGCCA
    GGCTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCCGGCTCCGGCCCGACCTTGTT
    CGCTCTGTGTGACAAGCCGGAAACCGCCCAGCGCGTTGCCGACTGGTTGGGTAAGAACTACCTGCAAAATCAGGAAGGTT
    TTGTTCATATTTGCCGGCTGGATACGGCGGGCGCACGAGTACTGGAAAACTAA
    . . .
    ```

## Structure
```
clovers/
├── bin/               # Binary data files
│   └── meta.bin           # Pre-designed heuristic model params
├── build/             # Build directory
├── include/           # Header files
│   ├── BioIO.hpp           # Bioinformatics I/O operation
│   ├── BioStruct.hpp       # Bioinformatics data structure
│   ├── BioUtil.hpp         # Bioinformatics utility
│   ├── cxxopt.hpp          # Command line option parser
│   ├── Encoding.hpp        # Sequence Z-curve encoding
│   ├── Model.hpp           # Model function (CDS, TIS)
│   └── svm.hpp             # LibSVM wrapper
├── Scripts/           # Python scripts
│   ├── Label.py            # Ab initio labeling ORFs
│   ├── README.md           # Scripts instruction
│   ├── Train.py            # Training heuristic models
│   ├── requirements.txt    # Python dependencies
│   ├── setup.py            # CPython extension setup
│   ├── Zcurve.cpp          # Z-curve encoding extension
│   └── Zcurve.pyi
├── src/               # Source files
│   ├── Main.cpp            # Program entry point
│   ├── BioIO.cpp
│   ├── BioStruct.cpp
│   ├── BioUtil.cpp
│   ├── Encoding.cpp
│   ├── Model.cpp
│   └── svm.cpp
├── example.fa         # Example FASTA file
├── LICENSE            # License file
├── Makefile           # Makefile for building
└── README.md          # Readme file
```

## Citation
Zetong Zhang, Yan Lin, Feng Gao. Ab initio prediction of overprinted genes using the Z-curve method.

## Server
- **CLOVERS Server**: Free available at [https://tubic.tju.edu.cn/clovers](https://tubic.tju.edu.cn/clovers/).

## Contact
- **Author**: Zetong Zhang, Yan Lin*, Feng Gao*
- **Homepage**: [https://tubic.tju.edu.cn/clovers](https://tubic.tju.edu.cn/clovers)
- **GitHub**: [https://github.com/zetong-zhang/clovers](https://github.com/zetong-zhang/clovers)
- **Issues**: Report bugs or feature requests on the [GitHub Issues](https://github.com/zetong-zhang/clovers/issues) page

## License
CLOVERS is distributed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

