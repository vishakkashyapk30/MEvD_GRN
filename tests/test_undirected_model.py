"""task='undirected' MeVD-GRN: permutation-invariant scores + a training step."""
import torch

from src.models.mevd_grn import MEvDGRN


def _toy(n=30, side=8, pair=6, seed=0):
    g = torch.Generator().manual_seed(seed)
    rna, atac = torch.randn(n, 3, generator=g), torch.randn(n, 3, generator=g)
    ei = torch.randint(0, n, (2, 80), generator=g)
    P = torch.randn(n, n, pair, generator=g)
    P = (P + P.transpose(0, 1)) / 2
    sig = torch.randn(n, 5, generator=g)
    m = MEvDGRN(hidden_dim=16, n_gnn_layers=2, task="undirected", side_in_dim=side, pair_dim=pair,
                use_motif=True)
    m.set_side_features(torch.randn(n, side, generator=g))
    m.set_pair_feats(P)
    return m, rna, atac, ei, sig


def test_scores_are_symmetric():
    m, rna, atac, ei, sig = _toy()
    m.eval()
    emb = m.encode(rna, atac, ei, ei, None, ei)
    u, v = torch.randint(0, 30, (200,)), torch.randint(0, 30, (200,))
    a = m.decode(emb, u, v, sig)
    b = m.decode(emb, v, u, sig)
    assert torch.allclose(a, b, atol=1e-5)


def test_forward_backward():
    m, rna, atac, ei, sig = _toy()
    m.train()
    emb = m.encode(rna, atac, ei, ei, None, ei)
    u, v = torch.randint(0, 30, (64,)), torch.randint(0, 30, (64,))
    loss = torch.nn.functional.binary_cross_entropy_with_logits(
        m.decode(emb, u, v, sig), torch.randint(0, 2, (64,)).float())
    loss.backward()
    assert all(torch.isfinite(p.grad).all() for p in m.parameters() if p.grad is not None)


def test_directed_default_unchanged():
    m = MEvDGRN(hidden_dim=16)
    assert m.task == "directed" and m.side_encoder is None and not hasattr(m, "node_fuse")
