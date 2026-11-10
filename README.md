# CLOVERS
<br/>

![LOGO](./logo.png)

<br/>

![C++](https://img.shields.io/badge/language-C%2B%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-brightgreen)
![Genomics](https://img.shields.io/badge/domain-genomics-purple)
![AVX](https://img.shields.io/badge/optimized-AVX%2FFMA-red)
![Static](https://img.shields.io/badge/build-static-lightgrey)
![Dependencies](https://img.shields.io/badge/dependencies-minimal-green)
[![GPLv3 License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![x86_64](https://img.shields.io/badge/arch-x86__64-green)

High-performance *ab initio* prokaryotic and phagic gene prediction and with CLOVERS !

## Contents
- **[Overview](#overview)** - Project introduction and features
- **[Setup](#setup)** - Installation guide for different operating systems
- **[Usage](#usage)** - Command-line interface and options
- **[Examples](#examples)** - Usage examples with command-line options
- **[Input](#input)** - Description of the input files and formats
- **[Output](#output)** - Description of the output files and formats
- **[Structure](#structure)** - Project file structure and organization
- **[Web Server](#web-server)** - Free available server for running tasks online
- **[Contact](#contact)** - Contact information for questions or support
- **[License](#license)** - GNU General Public License v3.0

## Overview
CLOVERS is a high-performance *ab initio* gene finder that revisits and reformulats geometric model-based paradigm pioneered by the ZCURVE family ([10.1093/nar/gkg254](https://doi.org/10.1093/nar/gkg254)) with a novel technical route. The tool delivers a state-of-the-art performance, enables the *de novo* discovery of overlapping genes, and runs fast with minimal memory. 

## Setup

Install from [Bioconda](https://bioconda.github.io) (Linux/macOS, 64-bit):

```bash
conda install -c conda-forge -c bioconda clovers
# or with micromamba/mamba
```

> **Tip:** the Bioconda build targets baseline x86-64 (no AVX), so it also runs on pre-Haswell CPUs but is somewhat slower than the AVX precompiled binaries.

Download the latest precompiled binary executable file (Linux/Windows/MacOS, 64-bit) from the [release page](https://tubic.tju.edu.cn/clovers/downloads), and decompress it to the directory of your choice. 

```bash
wget https://tubic.tju.edu.cn/clovers/static/pkg/clovers-1.2.3-linux64.tar.gz
# wget https://tubic.tju.edu.cn/clovers/static/pkg/clovers-1.2.3-win64.zip
# wget https://tubic.tju.edu.cn/clovers/static/pkg/clovers-1.2.3-macos64.tar.gz
tar -xvf clovers-1.2.3-linux64.tar.gz
# tar -xvf clovers-1.2.3-win64.zip
# tar -xvf clovers-1.2.3-macos64.tar.gz
```

If there is no distribution suitable for your operating system, you can modify the code yourself and compile and install it.   

```bash
git clone https://github.com/zetong-zhang/clovers.git
cd clovers
# do some modification if needed
make NATIVE=1 # use mingw32-make (MinGW) on Windows
```

## Usage
### Configuration
We recommend configure the environment variable `PATH` to include the directory of the executable binary file, such that you can run `clovers` directly in the terminal.

```bash
echo 'export PATH=$PATH:/path/to/clovers' >> ~/.bashrc
source ~/.bashrc
```

### Options

```bash
clovers.exe [OPTION...]
```

#### General Options
* `-h, --help`  
Print help menu and exit.

* `-q, --quiet`  
Run quietly with no stderr output (note that fatal error messages will still be printed under quiet mode).

* `-v, --version`  
Print version info and exit (useful when to quickly check whether it has been configured in the environment).

* `-T, --threads`  
Number of threads to use. Please set as a positive integer. (default: half)

#### Input/Output Options
* `-i, --input`   
Specify FASTA/GenBank/EMBL input file or their compressed versions (gzip). (default: stdin)  

* `-o, --output`  
Specify output file or '-' as standard output (stdout).

* `-f, --format`  
Select output format (gff, gbk, bed, gbk-full). (default: gff)

  **Note:** The bed format follows BED6 convention: `chrom  chromStart  chromEnd  name  score  strand`, with 0-based half-open coordinates (`chromStart = start - 1`, `chromEnd = end`) and score scaled from the predicted probability (0-1000).

* `-a, --faa`  
Write protein translations of genes to the selected file or '-' as stdout.

* `-d, --fna`  
Write nucleotide sequences of genes to the selected file or '-' as stdout.

#### CDS-Finder Options  
* `-g, --table`  
Specify a genetic codon table to use (1, 4, 11, 15, 16, 25 or auto). (default: 11)  

  **Note:** Auto-selection of codon table function can only determine whether TGA or TAG have been reprogrammed into sense codons.

* `-l, --minlen`  
Specify the mininum length (nt) of ORFs. (default: 90)

* `-c, --circ`  
Treat the default topology as circular. Note that the topology for each sequence can be set by the words "circular" or "linear" appear in the header line (FASTA/GenBank/EMBL), and this option will only take effect when the words "circular" or "linear" do not exist.

* `-p, --proc`  
Select prediction procedure (single or meta). In meta mode, ORFs are revised by TISA-CNN and scored with the built-in unified CDS-MLP; no genome‑specific model is built, to avoid false positives from contaminated sequences. Contig-edge ORFs down to 60 nt are kept in this mode.

* `-s, --thres`  
Specify putative gene probability score threshold. For chromosomes, the recommended setting is 0.5. For plasmids, viruses, and bacteriophages, the recommended setting is 0.4. (default: 0.5)

* `-t, --train`  
Write (if none exists) or use the specified training file.

#### TIS-Reviser Options  
* `-M, --model`  
Force the TIS reviser: `hmm` (TriTISA), `cnn` (TISA-CNN), or `none` (output the longest ORFs, most left 5'-end). Default: chosen by procedure and genome size — meta mode always uses TISA-CNN; single mode uses TISA-CNN for genomes up to 20 kb and TriTISA above 20 kb.

* `--maxiter`  
Max iteration times for RBS revision (TriTISA only). (default: 20)

* `-R,--rbs`  
Write (if none exists) or use the specified RBS model file (TriTISA only).

#### GOP-Reporter options:  
* `-L, --min-olen`  
Specify the mininum overlapped length between two ORFs. (default: 120)

* `-O, --overlap`  
Write overlapped genes to the selected file (GFF3 format).

* `-A, --amino`  
Write protein translations of overlapped genes to the selected file.

* `-D, --nucl`  
Write nucleotide sequences of overlapped genes to the selected file.

## Examples
```bash
# Predict genes in a FASTA file and output to GFF3
clovers -i example.fasta -o example.gff

# Use auto-detection of translation table
clovers -i input.fasta -o output.gff -g auto

# Generate protein and nucleotide sequences
clovers -i input.fasta -o output.gff -a proteins.faa -d genes.fna

# Detect overlapped genes and their sequences
clovers -i input.fasta -o output.gff -O overlap.gff -A overlap.faa -D overlap.fna

# Predict genes in metagenomes
clovers -i input.fasta -o output.gff -p meta

# Generate training file and RBS model file
clovers -i input.fasta -o output.gff -t training.dat -I rbs_model.bin
```

## Input
CLOVERS supports FASTA format input (stdin or file) or their compressed versions (gzip) with content like follow:
```
>NC_000913.3 circular
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
Or minimal EMBL format input (stdout or file) or their compressed versions (gzip) with content like follow:
```
ID   NC_000913.3
XX
FH   Key             Location/Qualifiers
XX
SQ   
     agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc        60
     tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg       120
     tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac       180
     . . .
//
```
## Output
- **GFF/GenBank/BED files**: Primary annotation results containing gene locations, scores and other attributes

    **GFF File Example:**
    ```
    ##gff-version 3
    # trans_tbl: 11
    # NC_000913.3	4641652 bp	circular	UNA	19-MAY-2026
    NC_000913.3	CLOVERS_v1.2.3	CDS	337	2799	0.974	+	0	ID=ORF000001;partial=00;r_score=0.991
    NC_000913.3	CLOVERS_v1.2.3	CDS	2801	3733	0.968	+	0	ID=ORF000002;partial=00;r_score=0.985
    . . .
    ```
    **GenBank File Example:**
    ```
    LOCUS       NC_000913.3          4641652 bp    DNA     circular UNA 19-MAY-2026
    DEFINITION  NC_000913.3
    FEATURES             Location/Qualifiers
         CDS             337..2799
                         /locus_tag=ORF000001
                         /transl_table=11
                         /codon_start=1
                         /note="Derived by protein-coding gene prediction method: CLOVERS_v1.2.3"
                         /translation="MRVLKFGG...LGV*"
         CDS             2801..3733
                         /locus_tag=ORF000002
                         /transl_table=11
                         /codon_start=1
                         /note="Derived by protein-coding gene prediction method: CLOVERS_v1.2.3"
                         /translation="MVKVYAPA...LEN*"
    . . .
    ORIGIN
    //
    ```
    **BED File Example (BED6, 0-based half-open coordinates, score 0-1000):**
    ```
    # NC_000913.3	4641652 bp	circular	UNA	19-MAY-2026
    NC_000913.3	336	2799	ORF000001	974	+
    NC_000913.3	2800	3733	ORF000002	968	+
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

## Structure

The repository is laid out as follows (`.exe` extensions appear on Windows builds; files and directories marked `*` are generated):

```text
CLOVERS/
├── Makefile                  # builds clovers(.exe) and scripts/labeler(.exe)
│                             #   (mingw32-make on Windows)
├── include/                  # C++ headers shared by all binaries
│   ├── Common.hpp            # core library: genome/ORF containers, Z-curve
│   │                         #   encoding, heuristic models (MLP / PMM /
│   │                         #   RBF-SVM wrappers), I/O utilities
│   ├── svm.hpp               # libsvm C interface (types, training, scoring)
│   ├── cnn.hpp               # embedded TISA-CNN and unified CDS-MLP weights
│   └── cxxopts.hpp           # vendored command-line option parser
├── src/
│   ├── Clovers.cpp           # main gene-finder program (single & meta
│   │                         #   procedures, output & GOL-Reporter)
│   ├── Meta.cpp              # embedded pre-trained heuristic (MLP) weights;
│   │                         #   a numeric-only array included by Clovers.cpp
│   ├── Labeler.cpp           # ab initio ORF Labeler — developer tool that
│   │                         #   auto-labels coding/non-coding ORFs to build
│   │                         #   the heuristic-model training data
│   ├── Eval.cpp              # standalone evaluator: compares predicted CDS
│   │                         #   features with a reference set (5'/3'-end match)
│   └── svm.cpp               # libsvm implementation used by the RBF-SVM model
├── scripts/                  # Python training / evaluation toolchain
│   ├── labeler(.exe) *       # compiled Labeler binary (training-data generator)
│   ├── 01.GC_Partition.py    # GC partitioning + TIS-CNN window dataset
│   ├── 02.Build_Samples.py   # labeled CDS-model training samples
│   ├── 03.Train_MLP.py       # trains the per-GC-bin heuristic MLP models
│   ├── 04.Train_CNN.py       # trains the TISA-CNN start-site model
│   ├── 05.Eval.py            # evaluates predicted vs. reference CDS features
│   ├── 06.Test_Speed.py      # runtime benchmark of clovers vs. other tools
│   ├── 07.Neg_Ctrl.py        # decoy (negative-control) genome generator
│   ├── Zcurve.pyx            # Cython Z-curve encoder used by the scripts
│   ├── setup.py              # builds the Zcurve extension in place
│   ├── requirements.txt
│   └── README.md             # per-script instructions
├── example.fa                # demo input genome (E. coli)
├── clovers(.exe) *           # compiled main binary
├── build/ *                  # intermediate objects produced by make
├── CHANGELOG.md
├── LICENSE
├── logo.png
└── README.md
```

The two executables are assembled by the root `Makefile` as:

| Binary | Translation units |
|--------|-------------------|
| `clovers(.exe)` | `Common.hpp` + `Meta.cpp` + `Clovers.cpp` (object) + `svm.cpp` (object) |
| `scripts/labeler(.exe)` | `Common.hpp` + `Labeler.cpp` (object) + `svm.cpp` (object) |

`Meta.cpp` contains no code of its own — it is the single embedded weight array of the heuristic parameters read by the MLP scoring in `Common.hpp`, and is compiled into `clovers(.exe)` only. `Eval.cpp` is a self-contained utility and is not part of either binary.

## Web Server
- Free available at [https://tubic.tju.edu.cn/clovers/](https://tubic.tju.edu.cn/clovers/).

## Contact
- **Authors**: Zetong Zhang, Zhisong You, Yan Lin*, Feng Gao*
- **Address**: No. 92 Weijin Road Nankai District, Tianjin, China, 300072
- **Emails**: ylin@tju.edu.cn | fgao@tju.edu.cn
- **Telephone**: +86-22-27402697

## License
CLOVERS is distributed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

