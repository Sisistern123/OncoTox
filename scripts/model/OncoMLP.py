import torch
import torch.nn as nn
import torch.optim as optim

class OncoMLP(nn.Module):
    def __init__(self, input_dim, hidden_dim=64, dropout_rate=0.5):
        super(OncoMLP, self).__init__()

        self.network = nn.Sequential(
            # Layer 1
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            # Layer 2
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout_rate),

            # Output Layer (Single continuous value for viability)
            nn.Linear(hidden_dim // 2, 1)
        )

    def forward(self, x):
        return self.network(x)