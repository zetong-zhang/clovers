import os
from Zcurve import encode
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import orfipy_core
from typing import Dict, List
from sklearn.decomposition import PCA
from sklearn.preprocessing import minmax_scale, StandardScaler
from sklearn.svm import SVC
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    print("Warning: scikit-learn-intelex is not available.")
    pass

def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1 / (1 + np.exp(-x))

def extract_orfs(
    records: Dict[str, str],
    starts: list = ["ATG", "GTG", "TTG"],
    stops: list = ["TAG", "TGA", "TAA"],
    minlen: int = 90
) -> List[dict]:
    """
    Extract ORFs from a DNA sequence using Orfipy.

    Parameters
    ----------
        records: Dict[str, str]
            DNA sequences, where keys are record IDs and values are sequences.
        starts: list, optional
            List of start codons. Default is ["ATG", "GTG", "TTG"].
        stops: list, optional
            List of stop codons. Default is ["TAG", "TGA", "TAA"].
        minlen: int, optional
            The minimum length of ORFs to extract. Default is 90.
    Returns
    -------
        List[dict]: List of ORF dictionaries, where each dictionary contains:
            - "id": ORF ID (format: {record_id}_{start}_{end}_{strand})
            - "seq": ORF sequence
            - "start": Start position (0-based)
            - "end": End position (0-based)
            - "strand": Strand ("+" or "-")
    """
    def refine_start(longest, patience: int = 75):
        """
        Refine the start position of the longest ORF.

        If the longest ORF starts with a standard start codon (ATG), return the ORF as is.
        Otherwise, search for the standard start codon within the patience range.
        If found, return the ORF starting from the found start codon.
        Otherwise, return the ORF starting from the first codon.

        Parameters
        ----------
            longest: str
                The longest ORF sequence.
            patience: int, optional
                The range to search for the start codon. Default is 75.
        Returns
        -------
            str: The ORF sequence starting from the refined start codon.
            int: The offset of the start codon in the ORF sequence.
        """
        if longest[:3].upper() == starts[0]:
            return longest, 0
        length = len(longest)
        offset = 0
        for idx in range(3, min(length-2, patience), 3):
            if longest[idx:idx+3].upper() == starts[0] and length - idx >= minlen:
                offset = idx
                break
        suborf = longest[offset:]
        return suborf, offset
    
    orfs: List[dict] = []
    for record_id, sequence in records.items():
        for start, end, strand, _ in orfipy_core.orfs(
            seq=sequence,
            starts=starts,
            stops=stops,
            minlen=minlen-3,
            include_stop=True
        ):
            longest = sequence[start:end]
            if strand == '-':
                longest = str(Seq(longest).reverse_complement())
            suborf, offset = refine_start(longest)
            if strand == '+':
                start = start + offset
            else:
                end = end - offset
            orfs.append({
                "id": record_id,
                "seq": suborf,
                "start": start,
                "end": end,
                "length": len(suborf),
                "strand": strand,
            })
    return orfs

def iter_train(
    features: np.ndarray, 
    positive_orfs: np.ndarray, 
    negative_orfs: np.ndarray,
    quiet: bool = False,
    max_iter: int = 20,
    up_thres: float = 0.7,
    down_thres: float = 0.31,
    ratio_thres: float = 8.0
):
    """
    Iteratively train an SVM classifier to label ORFs as positive or negative.

    Parameters
    ----------
        features: np.ndarray
            Feature matrix, where each row is a feature vector.
        positive_orfs: np.ndarray
            Indices of positive ORFs in the feature matrix.
        negative_orfs: np.ndarray
            Indices of negative ORFs in the feature matrix.
        quiet: bool, optional
            Whether to suppress print messages. Default is False.
        max_iter: int, optional
            Maximum number of iterations. Default is 20.
        up_thres: float, optional
            Threshold for positive ORFs. Default is 0.7.
        down_thres: float, optional
            Threshold for negative ORFs. Default is 0.3.
        ratio_thres: float, optional
            Threshold for the ratio of negative to positive ORFs. Default is 8.0.
    Returns
    -------
        pos_in: np.ndarray
            Indices of positive ORFs after convergence.
        neg_in: np.ndarray
            Indices of negative ORFs after convergence.

    """
    if not quiet:
        print(f"Iteration # 0: ")
    pos_features = features[positive_orfs]
    neg_features = features[negative_orfs]
    data = np.concatenate([pos_features, neg_features], axis=0)
    labels = np.concatenate([np.ones(len(pos_features)), 
                             np.zeros(len(neg_features))], axis=0)
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    svm = SVC(kernel='rbf', C=10.0).fit(data_scaled, labels)
    scores = sigmoid(svm.decision_function(scaler.transform(features)))
    pos_in = np.where(scores > up_thres)[0]
    neg_in = np.where(scores <= down_thres)[0]
    n_pos, n_neg = len(pos_in), len(neg_in)
    if not quiet:
        print(f"\t{n_pos} positive ORFs, {n_neg} negative ORFs")
    for i in range(max_iter):
        if not quiet:
            print(f"Iteration # {i+1}: ")
        should_break = False
        pos_features = features[pos_in]
        neg_features = features[neg_in]
        data = np.concatenate([pos_features, neg_features], axis=0)
        labels = np.concatenate([np.ones(len(pos_features)), 
                                 np.zeros(len(neg_features))], axis=0)
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        svm = SVC(kernel='rbf', C=1.0, class_weight='balanced').fit(data_scaled, labels)
        scores = sigmoid(svm.decision_function(scaler.transform(features)))
        pos_in = np.where(scores > up_thres)[0]
        neg_in = np.where(scores <= down_thres)[0]
        n_pos, n_neg = len(pos_in), len(neg_in)
        if n_neg / n_pos >= ratio_thres:
            should_break = True
        if not quiet:
            print(f"\t{n_pos} positive ORFs, {n_neg} negative ORFs")
        if should_break:
            print(f"Converged after {i+1} iterations.")
            break
    return pos_in, neg_in, scaler, svm

def to_gff(
    orfs: List[dict], 
    filename: str, 
    output_dir: str,
    label: str
) -> None:
    """
    Convert ORF dictionaries to GFF format.

    Parameters
    ----------
        orfs: List[dict]
            List of ORF dictionaries.
        filename: str
            Filename to save the GFF file.
        output_dir: str
            Directory to save the GFF file.
    """
    gff_output = os.path.join(output_dir, f"{label}_{filename}.gff")
    with open(gff_output, "w") as gff_file:
        gff_file.write("##gff-version 3\n")
        for idx, orf in enumerate(orfs):
            host_id, start, end, strand = orf["id"], orf["start"], orf["end"], orf["strand"]
            start = start + 1
            score = orf["score"] if "score" in orf else "."
            gff_file.write(f"{host_id}\tCLOVERS\tCDS\t{start}\t{end}\t{score}\t{strand}\t.\tID=orf{idx+1:06d}\n")

def to_fasta(
    orfs: List[dict], 
    filename: str, 
    output_dir: str,
    label: str
) -> None:
    """
    Convert ORF dictionaries to FASTA format.

    Parameters
    ----------
        orfs: List[dict]
            List of ORF dictionaries.
        filename: str
            Filename to save the FASTA file.
        output_dir: str
            Directory to save the FASTA file.
    """
    fasta_output = os.path.join(output_dir, f"{label}_{filename}.fasta")
    with open(fasta_output, "w") as fasta_file:
        for idx, orf in enumerate(orfs):
            host_id, seq = orf["id"], orf["seq"]
            fasta_file.write(f">{host_id}_{idx+1:06d}\n{seq}\n")


def visual_orf_flower(
    features: np.ndarray, 
    labels: np.ndarray, 
    output_dir: str, 
    filename: str
):
    """
    Visualize ORFs in a 2D flower plot using PCA and K-means clustering.

    Parameters
    ----------
        features: np.ndarray
            Z-curve features for ORFs.
        labels: np.ndarray
            Cluster labels for ORFs.
        pca_features: np.ndarray
            PCA-transformed features for ORFs.
        output_dir: str
            Directory to save the output plot.
        filename: str
            Filename to save the plot.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(features)
    pos_features = pca_features[labels == 1]
    neg_features = pca_features[labels == 0]
    ax.scatter(pos_features[:, 0], pos_features[:, 1], c='red', label='Positive ORFs', s=10)
    ax.scatter(neg_features[:, 0], neg_features[:, 1], c='blue', label='Negative ORFs', s=10)
    ax.set_title('K-means Clustering of Long ORFs')
    ax.set_xlabel('PCA Component 1')
    ax.set_ylabel('PCA Component 2')
    plt.savefig(os.path.join(output_dir, f"{filename}_clustering.png"))

def main():
    parser = argparse.ArgumentParser(
        description="Ab initio labelling ORFs using Z-curve method",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", help="Input FASTA file containing genome sequences", required=True)
    parser.add_argument("-o", "--output", help="Output directory for detailed results", required=True)
    parser.add_argument("-g", "--table", help="Genetic code table number", type=int, default=11)
    parser.add_argument("-m", "--minlen", help="Minimum length of ORFs to extract", type=int, default=90)
    parser.add_argument("-l", "--longorf", help="Minimum length of ORFs to be considered as long ORFs", type=int, default=300)
    parser.add_argument("-c", "--cluster", help="Cluster number for K-means clustering", type=int, default=3)
    parser.add_argument("-n", "--n_jobs", help="Number of threads to use", type=int, default=0)
    parser.add_argument("-q", "--quiet", help="Suppress verbose output", type=bool, default=False)
    args = parser.parse_args()
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    # Load input FASTA file
    filename = os.path.basename(args.input)
    if '.' in filename:
        filename = filename.rsplit('.', 1)[0]
    with open(args.input, "r") as fasta_file:
        records: Dict[str, str] = dict()
        overall_gc_count, overall_len = 0, 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            strseq = str(record.seq)
            gc_count = strseq.count('G') + strseq.count('C')
            overall_gc_count += gc_count
            overall_len += len(strseq)
            records[record.id] = strseq
        gc_content = overall_gc_count / overall_len
        if not args.quiet:
            print(f"Loaded {len(records)} sequences.\n\tTotal length: {overall_len} bp\n\tG+C content: {gc_content:.4f}")
    # Extract ORFs
    orfs = extract_orfs(records, stops=['TAG', 'TAA'] if args.table == 4 else ['TAA', 'TAG', 'TGA'], minlen=args.minlen)
    if not args.quiet:
        print(f"Extracted {len(orfs)} ORFs.")
    # Filter long ORFs
    orf_seqs = [orf["seq"] for orf in orfs]
    long_orfs = np.array([idx for idx, seq in enumerate(orf_seqs) if len(seq) >= args.longorf])
    if len(long_orfs) == 0:
        raise ValueError(f"No long ORFs found with length >= {args.longorf} (>= {args.minlen} bp)")
    elif not args.quiet:
        print(f"Found {len(long_orfs)} long ORFs.")
    # Encode ORFs using Z-curve method
    if not args.quiet:
        print(f"Encoding ORFs ...")
    features = encode(orf_seqs, n_jobs=args.n_jobs)
    long_features = features[long_orfs]
    long_features_scaled = minmax_scale(long_features)
    # Cluster long ORFs using K-means
    if not args.quiet:
        print(f"Clustering long ORFs ...")
    kmeans = KMeans(n_clusters=args.cluster, random_state=42)
    labels = kmeans.fit_predict(long_features_scaled)
    max_avg_len, max_cluster, positive_orfs = 0, None, None
    for cluster in range(args.cluster):
        cluster_orfs = long_orfs[np.where(labels == cluster)[0]]
        cluster_len = sum([orfs[idx]["length"] for idx in cluster_orfs]) / len(cluster_orfs)
        if cluster_len > max_avg_len:
            max_avg_len = cluster_len
            max_cluster = cluster
            positive_orfs = cluster_orfs
        if not args.quiet:
            print(f"\tCluster {cluster}: {len(cluster_orfs)} ORFs, average length: {int(cluster_len)} bp")
    negative_orfs = long_orfs[np.where(labels != max_cluster)[0]]
    if not args.quiet:
        print(f"Initial Positive ORFs: {len(positive_orfs)}")
        print(f"Initial Negative ORFs: {len(negative_orfs)}")
    pos_in, neg_in, scaler, svm = iter_train(features, positive_orfs, negative_orfs, args.quiet)
    # Visualize ORFs in a 2D flower plot
    long_labels = svm.predict(scaler.transform(long_features))
    visual_orf_flower(long_features_scaled, long_labels, args.output, filename)
    pos_orfs = [orfs[idx] for idx in pos_in]
    neg_orfs = [orfs[idx] for idx in neg_in]
    to_gff(pos_orfs, filename, args.output, "positive")
    to_fasta(pos_orfs, filename, args.output, "positive")
    to_gff(neg_orfs, filename, args.output, "negative")
    to_fasta(neg_orfs, filename, args.output, "negative")

if __name__ == "__main__":
    main()