# CLOVERS
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

## Project Overview
CLOVERS is an ab initio gene finder that utilizes the Z-curve method for robust gene prediction on prokaryotic and phage genomes. Its novel technical approach gives it an advantage in detecting more smaller overprinted genes, which are often missed by other tools.

## Setup
### Windows/Linux (x86_64)
Download the latest precompiled binary file from the [release page](https://github.com/zetong-zhang/CLOVERS/releases), and decompress it to the directory of your choice.
### Other Operating Systems
Download and compile the source code yourself.

Compilation Example: 
```bash
objcopy --input binary --output elf64-x86-64 --binary-architecture i386:x86-64 meta.bin meta.bin.o
g++ -DZLIB -o clovers -fopenmp -mavx -mfma meta.bin.o Main.cpp BioIO.cpp BioUtil.cpp Encoding.cpp Model.cpp svm.cpp -static -lz -O3
```
*Note:*  
(1) If your system doesn't have zlib installed, or you just don't want to use zlib, just remove the `-DZLIB` and `-lz` option;   
(2) If you don't want to use AVX and FMA function, just remove the `-mavx` and `-mfma` option.

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
Print help menu and exit. If no option is specified, print the help menu.  

* `-q, --quiet`  
Run quietly with no stderr output.  

* `-T, --threads`  
Number of threads to use. (default: all)

#### Input/Output Options
* `-i, --input`   
Specify FASTA/Genbank input file. (default: stdin)  

* `-o, --output`  
Specify output file. (default: stdout)

* `-f, --format`  
Select output format (gff, gbk). (default: gff)

* `-a, --faa`  
Write protein translations to the selected file.

* `-d, --fna`  
Write nucleotide sequences of genes to the selected file.

#### CLOVERS Options  
* `-g, --table`  
Specify a translation table to use. (default: 11)

* `-l, --minlen`  
Specify the mininum length of ORFs. (default: 90)

* `-c, --circ`  
Treat topology as circular.

* `-p, --proc`  
Select procedure (single or meta).

* `-s, --thres`  
Specify putative gene score threshold. (default: 0.5)

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

## Output Files
- **GFF3/GenBank files**: Primary annotation results containing gene positions, orientations, IDs, and other attributes
- **Protein sequence files (.faa)**: Translated amino acid sequences of predicted genes
- **Nucleotide sequence files (.fna)**: DNA sequences of predicted genes

## Usage Examples
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

## Contact
- **Author**: Zetong Zhang, Yan Lin, Feng Gao
- **Homepage**: [https://tubic.tju.edu.cn/clovers](https://tubic.tju.edu.cn/clovers)
- **GitHub**: [https://github.com/zetong-zhang/clovers](https://github.com/zetong-zhang/clovers)
- **Issues**: Report bugs or feature requests on the [GitHub Issues](https://github.com/zetong-zhang/clovers/issues) page

## License
CLOVERS is distributed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

