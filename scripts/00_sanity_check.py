#!/usr/bin/env python
"""Sanity checks (plan Step 15). Runs on SYNTHETIC data so it needs no downloads
and no GPU — validates shapes, forward/backward, loss, metrics, checkpoint I/O,
the curriculum trainer, and (if present) real processed data.

Usage:  python scripts/00_sanity_check.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data import graph_builder as gb
from src.data.dataset import CellTypeData
from src.data.preprocessing import coexpression_signatures
from src.evaluation.metrics import compute_all_metrics
from src.models.mevd_grn import MEvDGRN
from src.training.curriculum import CurriculumStage
from src.training.losses import bce_weighted_loss, focal_loss
from src.training.trainer import MEvDTrainer

OK, FAIL = "\033[92mPASS\033[0m", "\033[91mFAIL\033[0m"


def check(name, cond):
    print(f"  [{OK if cond else FAIL}] {name}")
    assert cond, f"sanity check failed: {name}"


def make_synthetic(n_genes=300, n_tf=40, n_cells=200, sig_dim=16, hidden=32, seed=0):
    rng = np.random.default_rng(seed)
    # Structured expression: latent factors -> correlated genes (so co-expression is real).
    n_factors = 8
    loadings = rng.standard_normal((n_genes, n_factors))
    factors = rng.standard_normal((n_cells, n_factors))
    X = (factors @ loadings.T) + 0.3 * rng.standard_normal((n_cells, n_genes))  # cells x genes
    X = np.maximum(X, 0.0).astype(np.float32)                                    # non-negative "counts"

    signatures = coexpression_signatures(X, sig_dim)                             # (n_genes, d)
    mean = X.mean(0); var = X.var(0); detect = (X > 0).mean(0)
    rna = torch.tensor(np.stack([mean, var, detect], 1), dtype=torch.float32)
    atac = torch.rand(n_genes, 3, generator=torch.Generator().manual_seed(seed))
    openness = torch.tensor(rng.random((n_genes, 4)), dtype=torch.float32)
    sig = torch.tensor(signatures, dtype=torch.float32)

    tf_indices = list(range(n_tf))
    gene_index = {f"G{i}": i for i in range(n_genes)}
    loc = set()
    while len(loc) < 600:
        loc.add((int(rng.integers(n_tf)), int(rng.integers(n_genes))))
    loc = sorted(loc)
    pert = loc[:250]
    dual = pert[:80]
    ev = {
        "localization": torch.tensor(list(zip(*loc)), dtype=torch.long),
        "perturbation": torch.tensor(list(zip(*pert)), dtype=torch.long),
        "dual_evidence": torch.tensor(list(zip(*dual)), dtype=torch.long),
    }
    coexpr = gb.build_coexpression_graph(signatures, k=10)
    tf_cand = gb.build_prior_graph(gene_index, tf_indices, atac.numpy(), rna.numpy(),
                                   top_k=30, signatures=signatures,
                                   openness=openness.numpy(), evidence=ev,
                                   exclude_positives=True)
    neg_pool = gb.create_negative_pool(ev, tf_indices, n_genes, cap=20000, seed=seed)
    data = CellTypeData("SYN", rna, atac, tf_cand, coexpr, sig, openness, gene_index,
                        torch.tensor(tf_indices), ev, neg_pool)
    return data, hidden, X, tf_cand, tf_indices


def _encode(model, data):
    return model.encode(data.rna_features, data.atac_features,
                        data.coexpr_edges, data.tf_candidate_edges)


def _decode(model, emb, data, tf_idx, tg_idx):
    return model.decode(emb, tf_idx, tg_idx, data.rna_signature, data.openness)


def main():
    torch.manual_seed(0)
    print("\n=== 1. Model forward / shapes (plan 15.4) ===")
    data, hidden, X, tf_cand, tf_indices = make_synthetic()
    model = MEvDGRN(hidden_dim=hidden, n_gnn_layers=2, dropout=0.2)
    print(f"  params = {model.count_parameters():,}")
    emb = _encode(model, data)
    h_rna, h_atac = emb
    check("encode returns two channels", isinstance(emb, tuple) and len(emb) == 2)
    check("channel shape (n_genes, hidden)", tuple(h_rna.shape) == (data.n_genes, hidden))
    tf_idx = data.evidence["localization"][0][:16]
    tg_idx = data.evidence["localization"][1][:16]
    logits = _decode(model, emb, data, tf_idx, tg_idx)
    probs = torch.sigmoid(logits)
    check("decode logits shape", logits.shape == (16,))
    check("probs in [0,1]", bool((probs >= 0).all() and (probs <= 1).all()))

    print("\n=== 1b. Co-expression signatures ~ correlation (plan 15.x) ===")
    corr = np.corrcoef(X, rowvar=False)                       # (n_genes, n_genes)
    sig = data.rna_signature.numpy()
    approx = sig @ sig.T
    iu = np.triu_indices(sig.shape[0], k=1)
    r = float(np.corrcoef(approx[iu], corr[iu])[0, 1])
    print(f"  signature-dot vs true corr: pearson r = {r:.3f}")
    check("signature dot product tracks correlation (r>0.7)", r > 0.7)

    print("\n=== 1c. TF-candidate graph is TF-specific ===")
    # Each TF must NOT get an identical neighbour set (the old global-ranking bug).
    nbr = {}
    src = tf_cand[0].tolist(); dst = tf_cand[1].tolist()
    for s, d in zip(src, dst):
        nbr.setdefault(s, set()).add(d)
    tfs_with_edges = [t for t in tf_indices if t in nbr]
    distinct = {frozenset(nbr[t]) for t in tfs_with_edges}
    print(f"  {len(tfs_with_edges)} TFs, {len(distinct)} distinct neighbour sets")
    check("TF neighbourhoods differ across TFs", len(distinct) > 1)

    print("\n=== 1d. Accessibility gate in [0,1] ===")
    a = torch.sigmoid(model.decoder.acc_gate(h_atac[tg_idx]))
    check("gate in [0,1]", bool((a >= 0).all() and (a <= 1).all()))

    print("\n=== 2. Losses (plan 15.5) ===")
    labels = torch.cat([torch.ones(8), torch.zeros(8)])
    lb = bce_weighted_loss(logits, labels, pos_weight=5.0)
    lf = focal_loss(logits, labels)
    check("bce finite & positive", torch.isfinite(lb) and lb.item() > 0)
    check("focal finite & >=0", torch.isfinite(lf) and lf.item() >= 0)

    print("\n=== 3. One training step / gradients (plan 15.6) ===")
    model.zero_grad()
    lb.backward()
    grads = [p.grad for p in model.parameters() if p.requires_grad]
    check("some grads non-None", any(gr is not None for gr in grads))
    check("no NaN grads", all(torch.isfinite(gr).all() for gr in grads if gr is not None))

    print("\n=== 4. Metrics ===")
    y = np.array([1, 1, 0, 0, 1, 0])
    s = np.array([0.9, 0.8, 0.2, 0.1, 0.6, 0.3])
    m = compute_all_metrics(y, s, n_total_candidates=100)
    check("aupr in [0,1]", 0 <= m["aupr"] <= 1)
    check("epr computed", np.isfinite(m["epr"]))

    print("\n=== 5. Evidence nesting helper (plan 15.2) ===")
    rep = gb.verify_nesting(data.evidence, ["localization", "perturbation", "dual_evidence"])
    check("dual ⊆ pert", rep["dual_evidence⊆perturbation"] > 0.95)
    check("dual ⊆ loc", rep["dual_evidence⊆localization"] > 0.95)

    print("\n=== 6. No-leakage in splits (plan 15.3) ===")
    sp = gb.create_edge_splits(data.evidence["dual_evidence"], data.negative_pool,
                               0.7, 0.15, 5, 5, seed=42)
    tr = set(zip(sp["train"]["pos"][0].tolist(), sp["train"]["pos"][1].tolist()))
    te = set(zip(sp["test"]["pos"][0].tolist(), sp["test"]["pos"][1].tolist()))
    check("train/test positives disjoint", len(tr & te) == 0)
    pos_all = set(zip(*[data.evidence[t][i].tolist() for t in data.evidence for i in (0, 1)][:0])) \
        if False else gb.edge_set(data.evidence["localization"]) | \
        gb.edge_set(data.evidence["perturbation"]) | gb.edge_set(data.evidence["dual_evidence"])
    negset = gb.edge_set(data.negative_pool)
    check("negative pool has no positives", len(negset & pos_all) == 0)

    print("\n=== 7. Checkpoint round-trip (plan 15.7) ===")
    model.eval()
    pred_args = (data.rna_features, data.atac_features, data.coexpr_edges,
                 data.tf_candidate_edges, tf_idx, tg_idx, data.rna_signature, data.openness)
    with torch.no_grad():
        out1 = model.predict(*pred_args)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "m.pt"
        torch.save({"model": model.state_dict()}, path)
        model2 = MEvDGRN(hidden_dim=hidden, n_gnn_layers=2, dropout=0.2)
        model2.load_state_dict(torch.load(path)["model"])
        model2.eval()
        with torch.no_grad():
            out2 = model2.predict(*pred_args)
    check("reloaded model identical outputs", torch.allclose(out1, out2, atol=1e-6))

    print("\n=== 8. Mini curriculum run (2+1+1 epochs) ===")
    cfg = _mini_cfg()
    splits = {t: gb.create_edge_splits(data.evidence[t], data.negative_pool, 0.7, 0.15, 5, 5, 42)
              for t in data.evidence}
    trainer = MEvDTrainer(MEvDGRN(hidden_dim=hidden, n_gnn_layers=2, dropout=0.2), data, cfg, "cpu")
    stages = [
        CurriculumStage("Stage1_Localization", "localization", 2, 1e-3, 5, False),
        CurriculumStage("Stage2_Perturbation", "perturbation", 1, 3e-4, 5, False),
        CurriculumStage("Stage3_DualEvidence", "dual_evidence", 1, 1e-4, 10, True),
    ]
    res = trainer.run_full_curriculum(splits, stages)
    check("curriculum produced results for all stages", len(res) == 3)
    final = trainer.evaluate_split(splits["dual_evidence"]["test"])
    check("final eval AUPR finite", np.isfinite(final["aupr"]))

    print("\n=== 9. Real processed data (optional) ===")
    found = False
    for ct in ("K562", "ESC"):
        d = Path(f"data/processed/{ct}")
        if (d / "rna_features_aligned.npy").exists():
            found = True
            rna = np.load(d / "rna_features_aligned.npy")
            atac = np.load(d / "atac_features_aligned.npy")
            check(f"{ct} rna/atac same n_genes", rna.shape[0] == atac.shape[0])
            check(f"{ct} rna shape (N,3)", rna.shape[1] == 3)
            if (d / "rna_signature.npy").exists():
                sigf = np.load(d / "rna_signature.npy")
                check(f"{ct} signature rows == n_genes", sigf.shape[0] == rna.shape[0])
            if (d / "openness.npy").exists():
                opf = np.load(d / "openness.npy")
                check(f"{ct} openness len == n_genes", opf.shape[0] == rna.shape[0])
    if not found:
        print("  (no processed data yet — run scripts/02_preprocess.py)")

    print("\n\033[92mALL SANITY CHECKS PASSED\033[0m\n")


def _mini_cfg():
    return {
        "data": {"seed": 42, "tier_hierarchy": ["localization", "perturbation", "dual_evidence"]},
        "training": {"weight_decay": 1e-4, "batch_size": 256, "clip_grad_norm": 1.0,
                     "early_stopping_patience": 5, "eval_every": 1, "loss_type": "bce_weighted",
                     "focal_alpha": 0.75, "focal_gamma": 2.0},
        "curriculum": {"use_hard_negatives": True, "use_memory_replay": False, "replay_weight": 0.1},
        "checkpoint_dir": tempfile.mkdtemp(),
    }


if __name__ == "__main__":
    main()
