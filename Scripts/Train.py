import argparse
import torch
import joblib
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau
from Bio import SeqIO
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from torch.utils.data import TensorDataset, DataLoader
from typing import Tuple

from Zcurve import encode

# Input dimension
INPUT_DIM = 189
# Hidden layer dimension
HIDDEN_DIM = 100

class HeuristicMLP(nn.Module):
    """
    MLP model for heuristic gene prediction.
    
    Parameters:
    -----------
    input_dim: int, optional (default=189)
        Dimension of the input Z-curve features.
    """
    def __init__(self, input_dim: int = INPUT_DIM) -> None:
        super(HeuristicMLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, HIDDEN_DIM)
        self.fc2 = nn.Linear(HIDDEN_DIM, 1)
        self.relu = nn.ReLU()
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.relu(self.fc1(x))
        x = self.fc2(x)
        return x
    
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """
        Predict the probability.

        Parameters
        ----------
        x: torch.Tensor
            Input tensor of shape (batch_size, input_dim).
        
        Returns:
        --------
        torch.Tensor
            Tensor of shape (batch_size,) containing the predicted probabilities.
        """
        return torch.sigmoid(self.forward(x))
    
    def predict(self, x: torch.Tensor, threshold: float = 0.5) -> torch.Tensor:
        """
        Predict the class labels.
        
        Parameters:
        -----------
        x: torch.Tensor
            Input tensor of shape (batch_size, input_dim).
        threshold: float, optional (default=0.5)
            Threshold for class prediction.
        
        Returns:
        --------
        torch.Tensor
            Tensor of shape (batch_size,) containing the predicted class labels (0 or 1).
        """
        return (self.predict_proba(x) >= threshold).float()

def train_model(
    model, train_loader, 
    val_loader=None, 
    num_epochs=20, 
    learning_rate=0.001,
    threshold=0.5,
    quiet: bool = False
):
    """
    Train the HeuristicMLP model.
    
    Parameters:
    -----------
    model: HeuristicMLP
        The HeuristicMLP model to be trained.
    train_loader: DataLoader
        DataLoader for the training dataset.
    val_loader: DataLoader, optional (default=None)
        DataLoader for the validation dataset.
    num_epochs: int, optional (default=100)
        Number of epochs for training.
    learning_rate: float, optional (default=0.001)
        Learning rate for the optimizer.
    threshold: float, optional (default=0.5)
        Threshold for class prediction.
    """
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(DEVICE)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)
    best_val_f1 = 0.0
    best_model_state = None
    for epoch in range(num_epochs):
        ### Training Stage ###
        model.train()
        train_loss = 0.0
        correct, tp, fp, fn, total = 0, 0, 0, 0, 0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(DEVICE), batch_y.to(DEVICE)
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y.unsqueeze(1))
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * batch_X.size(0)
            predicted = (outputs >= threshold).float()
            correct += (predicted == batch_y.unsqueeze(1)).sum().item()
            tp += ((predicted == 1) & (batch_y.unsqueeze(1) == 1)).sum().item()
            fp += ((predicted == 1) & (batch_y.unsqueeze(1) == 0)).sum().item()
            fn += ((predicted == 0) & (batch_y.unsqueeze(1) == 1)).sum().item()
            total += batch_y.size(0)
        
        train_loss /= len(train_loader.dataset)
        train_acc = correct / total
        train_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        train_recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        train_f1 = 2 * train_precision * train_recall / (train_precision + train_recall) if (train_precision + train_recall) > 0 else 0
        
        ### Validation Stage ###
        model.eval()
        val_loss = 0.0
        correct, tp, fp, fn, total = 0, 0, 0, 0, 0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(DEVICE), batch_y.to(DEVICE)
                outputs = model(batch_X)
                loss = criterion(outputs, batch_y.unsqueeze(1))
                val_loss += loss.item() * batch_X.size(0)
                predicted = (outputs >= threshold).float()
                correct += (predicted == batch_y.unsqueeze(1)).sum().item()
                tp += ((predicted == 1) & (batch_y.unsqueeze(1) == 1)).sum().item()
                fp += ((predicted == 1) & (batch_y.unsqueeze(1) == 0)).sum().item()
                fn += ((predicted == 0) & (batch_y.unsqueeze(1) == 1)).sum().item()
                total += batch_y.size(0)
        
        val_loss /= len(val_loader.dataset)
        val_acc = correct / total
        val_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        val_recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        val_f1 = 2 * val_precision * val_recall / (val_precision + val_recall) if (val_precision + val_recall) > 0 else 0
        scheduler.step(val_loss)

        if not quiet:
            print(f"Epoch {epoch+1}/{num_epochs}, "
                  f"Train Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f}, Train Precision: {train_precision:.4f}, Train Recall: {train_recall:.4f}, Train F1: {train_f1:.4f}, "
                  f"Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}, Val Precision: {val_precision:.4f}, Val Recall: {val_recall:.4f}, Val F1: {val_f1:.4f}")
        
        if val_f1 > best_val_f1:
            best_val_f1 = val_f1
            best_model_state = model.state_dict()
    return best_model_state

def load_pair_datasets(
    pos_file: str, 
    neg_file: str, 
    test_size: float = 0.1,
    quiet: bool = False
) -> Tuple[DataLoader, DataLoader, StandardScaler]:
    """
    Load and preprocess paired datasets of CDS and ncORF.
    
    Parameters:
    -----------
    pos_file: str
        Path to the FASTA file containing CDS sequences.
    neg_file: str
        Path to the FASTA file containing ncORF sequences.
    quiet: bool, optional (default=False)
        If True, suppress verbose output.
    
    Returns:
    --------
    train_loader: DataLoader
        DataLoader for the training dataset.
    val_loader: DataLoader
        DataLoader for the validation dataset.
    scaler: StandardScaler
        StandardScaler fitted on the training features.
    """
    labels = []
    # Load positive sequences
    pos_sequences = [str(record.seq) for record in SeqIO.parse(pos_file, 'fasta')]
    pos_features = encode(pos_sequences)
    labels.extend([1] * len(pos_features))
    if not quiet:
        print(f"Loaded {len(pos_features)} positive sequences.")
    # Load negative sequences
    neg_sequences = [str(record.seq) for record in SeqIO.parse(neg_file, 'fasta')]
    neg_features = encode(neg_sequences)
    labels.extend([0] * len(neg_features))
    if not quiet:
        print(f"Loaded {len(neg_features)} negative sequences.")
    # Split the dataset into training and validation sets
    features = np.concatenate([pos_features, neg_features], axis=0)
    labels = np.array(labels)
    X_train, X_val, y_train, y_val = train_test_split(features, labels, test_size=test_size, random_state=42)
    # Standardize the features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    # Convert to PyTorch tensors
    X_train = torch.FloatTensor(X_train)
    X_val = torch.FloatTensor(X_val)
    y_train = torch.FloatTensor(y_train)
    y_val = torch.FloatTensor(y_val)
    # Create DataLoaders
    train_dataset = TensorDataset(X_train, y_train)
    val_dataset = TensorDataset(X_val, y_val)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)

    return train_loader, val_loader, scaler

def main():
    parser = argparse.ArgumentParser(description="Train Heuristic MLP model for CLOVERS")
    parser.add_argument("-p", "--pos_file", type=str, required=True, help="Path to FASTA file of CDS")
    parser.add_argument("-n", "--neg_file", type=str, required=True, help="Path to FASTA file of ncORF")
    parser.add_argument("-o", "--output", type=str, default="model.pth", help="Path to save the trained model")
    parser.add_argument("-e", "--num_epochs", type=int, default=20, help="Number of training epochs")
    parser.add_argument("-l", "--learning_rate", type=float, default=0.001, help="Learning rate for optimizer")
    parser.add_argument("-t", "--threshold", type=float, default=0.5, help="Threshold for binary classification")
    parser.add_argument("-s", "--test_size", type=float, default=0.1, help="Proportion of dataset to use for validation")
    parser.add_argument("-q", "--quiet", help="Suppress verbose output", type=bool, default=False)
    args = parser.parse_args()
    # Load datasets
    train_loader, val_loader, scaler = load_pair_datasets(args.pos_file, args.neg_file, test_size=args.test_size, quiet=args.quiet)
    # Initialize model
    model = HeuristicMLP()
    # Train the model
    best_model_state = train_model(
        model, train_loader, val_loader,
        num_epochs=args.num_epochs,
        learning_rate=args.learning_rate,
        threshold=args.threshold,
        quiet=args.quiet
    )
    best_model = {
        "state_dict": best_model_state,
        "scaler": scaler,
        "threshold": args.threshold,
    }
    joblib.dump(best_model, args.output)
    if not args.quiet:
        print(f"Best model saved to {args.output}")

if __name__ == "__main__":
    main()