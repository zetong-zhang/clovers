# Instruction of Scripts
Python scripts for generating the heuristic model of CLOVERS. These scripts belong to the **developer tools** and are not related to the use of the final built program that is released.
## Step I: Installation of third-party packages
| Package              | Description                         |
|----------------------|-------------------------------------|
| biopython            | Python libraries for bioinformatics |
| scikit-learn         | Machine learning library for Python |
| matplotlib           | Python plotting library             |
| orfipy               | Python library for ORF extraction   |
| torch                | PyTorch library for deep learning   |
| scikit-learn-intelex | Intel extension for scikit-learn    |
```bash
cd Scripts
pip install -r requirements.txt
```
## Step II: Setup of Zcurve C++ Extension
Build the C++ extension to accelerate the Z-curve calculation of Python scripts.
```bash
python setup.py build_ext --inplace
```
## Step III: Ab initio labelling on genomes
Use ab initio labelling on anomymous genomes to generate the training data (coding and non-coding ORF samples) for the heuristic model of CLOVERS:
```bash
python Label.py -i genome.fasta -o output_dir # [other options]
```
### Input
For the detailed info of the input options, please refer to the help message:
```bash
python Label.py -h
```
### Output
The output directory will contain the following files:
- `positive_{filename}.fasta`: Coding ORF sequences
- `negative_{filename}.fasta`: Non-coding ORF sequences
- `positive_{filename}.gff`: Coding ORF annotations
- `negative_{filename}.gff`: Non-coding ORF annotations
- `{filename}_clustering.csv`: Clustering results of ORFs

## Step IV: Training of heuristic model
Use the training data to train the heuristic model of CLOVERS:
```bash
python Train.py -p positive.fasta -n negative.fasta -o model.pth # [other options]
```
### Input
For the detailed info of the input options, please refer to the help message:
```bash
python Train.py -h
```
### Output
The output file `model.pth` will contain the trained model.
```json
model = {
    "state_dict": model_state,  # State dictionary of the trained model (torch.nn.Module.state_dict())
    "scaler": scaler,           # StandardScaler object for feature scaling (sklearn.preprocessing.StandardScaler)
    "threshold": threshold,     # Threshold value for classification (float)
}
```
