import torch
import torch.nn as nn

import config


class BeamAutoencoder(nn.Module):
    """Small MLP autoencoder over the max-normalized trace features.

    Input width follows config.ae_input_dim() (included strips, plus the
    appended derivative when enabled). Trained on beam-dominated input
    only; reconstruction MSE is the anomaly score. Linear output (no
    squashing).
    """

    def __init__(self,
                 n_inputs=None,
                 hidden_dims=config.AE_HIDDEN_DIMS,
                 latent_dim=config.AE_LATENT_DIM):
        super().__init__()
        if n_inputs is None:
            n_inputs = config.ae_input_dim()
        enc = []
        prev = n_inputs
        for h in hidden_dims:
            enc.append(nn.Linear(prev, h))
            enc.append(nn.ReLU())
            prev = h
        enc.append(nn.Linear(prev, latent_dim))
        self.encoder = nn.Sequential(*enc)

        dec = []
        prev = latent_dim
        for h in reversed(hidden_dims):
            dec.append(nn.Linear(prev, h))
            dec.append(nn.ReLU())
            prev = h
        dec.append(nn.Linear(prev, n_inputs))
        self.decoder = nn.Sequential(*dec)

    def forward(self, x):
        return self.decoder(self.encoder(x))

    @torch.no_grad()
    def reconstruction_error(self, x):
        """Per-trace mean squared reconstruction error (the anomaly score)."""
        return torch.mean((self.forward(x) - x)**2, dim=1)
