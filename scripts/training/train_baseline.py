import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader


from scripts.model.dataset import ScGPTDrugDataset
from scripts.model.OncoMLP import OncoMLP

# --- 1. Configuration ---
FILE_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
TARGET_DRUG = "paclitaxel"
USE_REP = "X_pca"

BATCH_SIZE = 128
EPOCHS = 50

# --- 2. Data Loading (Using the Pre-computed Splits) ---
train_dataset = ScGPTDrugDataset(h5ad_path=FILE_PATH, target_drug=TARGET_DRUG, use_rep=USE_REP, split="train")
train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)

val_dataset = ScGPTDrugDataset(h5ad_path=FILE_PATH, target_drug=TARGET_DRUG, use_rep=USE_REP, split="val")
val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

# --- 3. Model Setup ---
sample_x, _ = train_dataset[0]
INPUT_DIM = sample_x.shape[0]

model = OncoMLP(input_dim=INPUT_DIM)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)

device = torch.device("mps" if torch.backends.mps.is_available() else "cuda" if torch.cuda.is_available() else "cpu")
print(f"Training on device: {device}")
model.to(device)

# --- 4. Training & Validation Loop ---
print(f"Starting Baseline (PCA) Training for {TARGET_DRUG}...")
for epoch in range(EPOCHS):
    # -- Training Phase --
    model.train()
    running_train_loss = 0.0

    for batch_X, batch_y in train_loader:
        batch_X, batch_y = batch_X.to(device), batch_y.to(device)

        optimizer.zero_grad()
        predictions = model(batch_X)
        loss = criterion(predictions, batch_y)
        loss.backward()
        optimizer.step()

        running_train_loss += loss.item() * batch_X.size(0)

    avg_train_loss = running_train_loss / len(train_dataset)

    # -- Validation Phase --
    model.eval()
    running_val_loss = 0.0

    with torch.no_grad():
        for batch_X, batch_y in val_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            predictions = model(batch_X)
            loss = criterion(predictions, batch_y)
            running_val_loss += loss.item() * batch_X.size(0)

    avg_val_loss = running_val_loss / len(val_dataset)

    if (epoch + 1) % 5 == 0 or epoch == 0:
        print(f"Epoch [{epoch+1:02d}/{EPOCHS}] | Train MSE: {avg_train_loss:.4f} | Val MSE: {avg_val_loss:.4f}")

print("scGPT Training complete!")