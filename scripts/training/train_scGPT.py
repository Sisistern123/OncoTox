from torch.utils.data import DataLoader

from scripts.model.OncoMLP import OncoMLP
from scripts.model.dataset import ScGPTDrugDataset
from scripts.training.training_utils import TrainConfig, train_model

FILE_PATH = "/Users/selin/Desktop/OncoTox/data/scRNAseq_SCP542/metadata/SCP542_CCLE_scGPT_human_embeddings_with_targets.h5ad"
TARGET_DRUG = "paclitaxel"
USE_REP = "X_scGPT"

BATCH_SIZE = 128


def main():
    train_dataset = ScGPTDrugDataset(
        h5ad_path=FILE_PATH, target_drug=TARGET_DRUG, use_rep=USE_REP, split="train"
    )
    val_dataset = ScGPTDrugDataset(
        h5ad_path=FILE_PATH, target_drug=TARGET_DRUG, use_rep=USE_REP, split="val"
    )

    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)

    sample_x, _ = train_dataset[0]
    input_dim = sample_x.shape[0]

    # scGPT embedding (512-dim) tolerates slightly more capacity than the PCA baseline,
    # but we keep depth modest to avoid the prior memorization regime.
    model = OncoMLP(
        input_dim=input_dim,
        hidden_dims=(128, 64),
        dropout_rate=0.5,
        input_dropout=0.1,
        norm="layer",
    )

    config = TrainConfig(
        epochs=50,
        lr=1e-3,
        weight_decay=1e-3,
        grad_clip=1.0,
        scheduler_patience=3,
        early_stop_patience=10,
        log_every=5,
        seed=42,
        loss="mse",
    )

    print(f"Starting scGPT Training for {TARGET_DRUG}...")
    train_model(model, train_loader, val_loader, config=config, tag="scGPT")


if __name__ == "__main__":
    main()
