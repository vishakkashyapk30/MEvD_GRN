#!/usr/bin/env python
"""End-to-end preprocessing for one cell type (plan Step 4-5).

scRNA -> gene universe (RNA ∩ networks) -> scATAC gene-activity ->
evidence edges (+nesting check) -> prior graph -> negative pool -> edge splits.

Usage:
  python scripts/02_preprocess.py --config configs/k562.yaml
  python scripts/02_preprocess.py --config configs/esc.yaml

Graph-only rebuild (2026-10-08): if the config sets `paths.reuse_processed_from`
(an existing processed dir), every label-free artifact (features, signatures,
openness, gene index, co-expression graph, evidence, negative pool, FM
embeddings) is copied unchanged from that dir and ONLY the TF-candidate graph
is rebuilt, with the config's `data.prior_exclude_positives`. Splits are not
rewritten: they are regenerated in memory and checked edge-for-edge against
the existing `<splits_dir>/<cell_type>_<tier>_splits.pt`. Used by
configs/sweep/k562_fm_h384_l2_labelfree.yaml so the label-free run differs from
the paper's processed data in the graph alone (re-running the raw pipeline
would also re-draw the ARPACK start vector in `coexpression_signatures`).
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data import graph_builder as gb
from src.data import motif_scan
from src.data.preprocessing import _resolve_glob, preprocess_scatac, preprocess_scrna
from src.utils.io import ensure_dir, load_config, load_json, save_json

# Files copied verbatim by the graph-only rebuild (none depends on labels except
# evidence_*.pt / negative_pool.pt, which are the labels and are unchanged).
_REUSED_FILES = ["rna_features_aligned.npy", "atac_features_aligned.npy", "rna_signature.npy",
                 "openness.npy", "gene_index.json", "tf_indices.json", "coexpr_edges.pt",
                 "negative_pool.pt", "fm_gene_embeddings.npy"]


def rebuild_graphs_only(cfg: dict, cell_type: str, src_dir: Path, out_dir: Path,
                        splits_dir: Path) -> None:
    dcfg = cfg["data"]
    if src_dir.resolve() == out_dir.resolve():
        raise ValueError("reuse_processed_from must differ from processed_dir")
    mp = src_dir / "motif_edges.pt"
    if mp.exists() and torch.load(mp).shape[1]:
        raise NotImplementedError("graph-only rebuild cannot rebuild a non-empty motif graph "
                                  "(needs per-gene peak coordinates); run the full pipeline")
    exclude = bool(dcfg.get("prior_exclude_positives", False))
    print(f"\n########## Graph-only rebuild {cell_type}: {src_dir} -> {out_dir} "
          f"(prior_exclude_positives={exclude}) ##########", flush=True)

    gene_index = load_json(src_dir / "gene_index.json")
    tf_indices = load_json(src_dir / "tf_indices.json")
    rna = np.load(src_dir / "rna_features_aligned.npy")
    atac = np.load(src_dir / "atac_features_aligned.npy")
    sig = np.load(src_dir / "rna_signature.npy")
    openness = np.load(src_dir / "openness.npy")
    tiers = list(cfg["paths"]["networks"])
    evidence = {t: torch.load(src_dir / f"evidence_{t}.pt") for t in tiers
                if (src_dir / f"evidence_{t}.pt").exists()}
    top_k = int(dcfg.get("tf_candidate_topk", dcfg["top_k_prior_targets"]))

    # Same call as the full pipeline; with exclude=True it must reproduce the source graph.
    src_graph = torch.load(src_dir / "tf_candidate_edges.pt")
    if not exclude:
        check = gb.build_prior_graph(gene_index, tf_indices, atac, rna, top_k=top_k,
                                     signatures=sig, openness=openness, evidence=evidence,
                                     exclude_positives=True)
        same = check.shape == src_graph.shape and bool(torch.equal(check, src_graph))
        print(f"[check] rebuilding with exclude_positives=True reproduces the source "
              f"TF-candidate graph: {same}", flush=True)
        if not same:
            raise RuntimeError("source arrays do not reproduce the source TF-candidate graph")
    tf_cand = gb.build_prior_graph(gene_index, tf_indices, atac, rna, top_k=top_k,
                                   signatures=sig, openness=openness, evidence=evidence,
                                   exclude_positives=exclude)

    for name in _REUSED_FILES:
        if (src_dir / name).exists():
            shutil.copy2(src_dir / name, out_dir / name)
    for t in evidence:
        shutil.copy2(src_dir / f"evidence_{t}.pt", out_dir / f"evidence_{t}.pt")
    torch.save(tf_cand, out_dir / "tf_candidate_edges.pt")
    torch.save(torch.zeros((2, 0), dtype=torch.long), out_dir / "motif_edges.pt")

    # Splits: regenerate deterministically and require edge-identity with the
    # existing files, which the run will load. Nothing in splits_dir is written.
    neg_pool = torch.load(src_dir / "negative_pool.pt")
    active = {t: e for t, e in evidence.items() if e.shape[1] > 0}
    regen = gb.create_global_edge_splits(
        active, neg_pool, train_ratio=float(dcfg["train_ratio"]), val_ratio=float(dcfg["val_ratio"]),
        neg_train_ratio=int(dcfg["neg_train_ratio"]), neg_eval_ratio=5, seed=int(dcfg["seed"]))
    split_check = {}
    for t, sp in regen.items():
        old = torch.load(splits_dir / f"{cell_type}_{t}_splits.pt")
        for part in ("train", "val", "test"):
            for key in ("pos", "neg"):
                same = bool(torch.equal(sp[part][key], old[part][key]))
                split_check[f"{t}.{part}.{key}"] = same
                if not same:
                    raise RuntimeError(f"regenerated split {t}.{part}.{key} differs from "
                                       f"{splits_dir}/{cell_type}_{t}_splits.pt")
    print(f"[check] regenerated splits are edge-identical to {splits_dir}/{cell_type}_*: "
          f"{all(split_check.values())} ({len(split_check)} tensors)", flush=True)

    summary = load_json(src_dir / "summary.json")
    summary.update({
        "tf_candidate_edges": int(tf_cand.shape[1]), "motif_edges": 0,
        "prior_exclude_positives": exclude, "graph_only_rebuild_of": str(src_dir),
        "splits_dir": str(splits_dir), "splits_identical_to_existing": all(split_check.values()),
    })
    save_json(summary, out_dir / "summary.json")
    print(f"\n[done] {cell_type} graph-only rebuild -> {out_dir}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cell_type", default=None, help="override config cell_type")
    args = ap.parse_args()

    cfg = load_config(args.config)
    cell_type = args.cell_type or cfg["cell_type"]
    dcfg = cfg["data"]
    out_dir = ensure_dir(cfg["paths"]["processed_dir"])
    splits_dir = ensure_dir(cfg["paths"].get("splits_dir", "data/splits"))
    if cfg["paths"].get("reuse_processed_from"):
        rebuild_graphs_only(cfg, cell_type, Path(cfg["paths"]["reuse_processed_from"]),
                            out_dir, splits_dir)
        return
    network_paths = cfg["paths"]["networks"]

    # shared cfg for preprocessing functions
    pcfg = dict(dcfg)
    pcfg["genome"] = cfg.get("genome", "hg38")
    pcfg["atac"] = cfg.get("atac", {})

    print(f"\n########## Preprocessing {cell_type} ##########", flush=True)

    # 1. scRNA -------------------------------------------------------------
    rna_path = _resolve_glob(cfg["paths"]["single_cell"]["rna_glob"])
    print(f"[load] scRNA: {rna_path}", flush=True)
    rna_feats, rna_names, hvg, rna_sig = preprocess_scrna(rna_path, pcfg)
    rna_row = {g: i for i, g in enumerate(rna_names)}

    # 2. gene universe = RNA ∩ networks -----------------------------------
    universe, tf_names = gb.build_gene_universe(rna_names, network_paths)
    if not universe:
        raise RuntimeError("Empty gene universe — check gene-symbol casing / network paths.")
    gene_index = {g: i for i, g in enumerate(universe)}
    rna_aligned = np.stack([rna_feats[rna_row[g]] for g in universe]).astype(np.float32)
    sig_aligned = np.stack([rna_sig[rna_row[g]] for g in universe]).astype(np.float32)
    print(f"[universe] {len(universe)} genes (RNA ∩ networks); "
          f"{len(tf_names)} TF symbols in networks", flush=True)

    # 3. scATAC gene-activity + locus openness (already in universe order) -
    try:
        atac_path = _resolve_glob(cfg["paths"]["single_cell"]["atac_glob"])
        print(f"[load] scATAC: {atac_path}", flush=True)
    except FileNotFoundError:
        atac_path = None
        print("[load] no scATAC file found; ATAC features will be zero.", flush=True)
    atac_aligned, openness, gene_top_peaks = preprocess_scatac(atac_path, universe, gene_index, pcfg)

    # 4. evidence edges + nesting -----------------------------------------
    evidence = gb.build_evidence_edges(network_paths, gene_index)
    if dcfg.get("dual_evidence_mode") == "reconstruct":
        evidence["dual_evidence"] = gb.reconstruct_dual(evidence)
        print(f"[dual] reconstructed dual = loc ∩ pert: {evidence['dual_evidence'].shape[1]} edges", flush=True)
    nesting = gb.verify_nesting(evidence, dcfg["tier_hierarchy"])

    # 5. message-passing graphs -------------------------------------------
    #   (a) co-expression kNN (gene-gene); (b) TF-specific accessible+co-expressed
    #   candidates; (c) real TF motif hits (plan.md Section 2), if a genome
    #   FASTA + JASPAR motif file are configured -- empty otherwise.
    tf_indices = sorted(gene_index[g] for g in tf_names if g in gene_index)
    coexpr_edges = gb.build_coexpression_graph(sig_aligned, k=int(dcfg.get("coexpr_knn_k", 20)))
    tf_candidate_edges = gb.build_prior_graph(
        gene_index, tf_indices, atac_aligned, rna_aligned,
        top_k=int(dcfg.get("tf_candidate_topk", dcfg["top_k_prior_targets"])),
        signatures=sig_aligned, openness=openness, evidence=evidence,
        exclude_positives=bool(dcfg.get("prior_exclude_positives", False)))
    atac_cfg = cfg.get("atac", {}) or {}
    inv_gene_index = {i: g for g, i in gene_index.items()}
    tf_names_by_index = {i: inv_gene_index[i] for i in tf_indices}
    motif_edges = motif_scan.build_motif_graph(
        gene_index, tf_indices, tf_names_by_index, gene_top_peaks,
        fasta_path=atac_cfg.get("genome_fasta"), jaspar_path=atac_cfg.get("jaspar_pfm_path"),
        score_frac=float(atac_cfg.get("motif_score_frac", 0.75)),
        evidence=evidence, exclude_positives=bool(dcfg.get("prior_exclude_positives", False)),
        candidate_edges=tf_candidate_edges)

    # 6. negative pool: random TF×gene non-positive edges -----------------
    neg_pool = gb.create_negative_pool(evidence, tf_indices, len(universe),
                                       cap=int(dcfg.get("negative_pool_cap", 1_000_000)),
                                       seed=int(dcfg["seed"]))

    # 7. save processed ----------------------------------------------------
    np.save(out_dir / "rna_features_aligned.npy", rna_aligned)
    np.save(out_dir / "atac_features_aligned.npy", atac_aligned)
    np.save(out_dir / "rna_signature.npy", sig_aligned)
    np.save(out_dir / "openness.npy", openness)
    save_json(gene_index, out_dir / "gene_index.json")
    save_json(tf_indices, out_dir / "tf_indices.json")
    torch.save(tf_candidate_edges, out_dir / "tf_candidate_edges.pt")
    torch.save(coexpr_edges, out_dir / "coexpr_edges.pt")
    torch.save(motif_edges, out_dir / "motif_edges.pt")
    torch.save(neg_pool, out_dir / "negative_pool.pt")
    for tier, e in evidence.items():
        torch.save(e, out_dir / f"evidence_{tier}.pt")

    # 8. edge splits (GLOBAL, hierarchy-safe -- see graph_builder docstring) -
    # Evidence tiers are heavily nested (verified in step 4 above), so tiers
    # cannot be split independently without an edge held out in one tier
    # leaking into another tier's training positives. create_global_edge_splits
    # assigns every unique (TF, gene) pair to exactly one split ONCE, across
    # the union of all tiers, then derives each tier's train/val/test from
    # that single assignment.
    active_evidence = {t: e for t, e in evidence.items() if e.shape[1] > 0}
    splits_per_tier = gb.create_global_edge_splits(
        active_evidence, neg_pool,
        train_ratio=float(dcfg["train_ratio"]), val_ratio=float(dcfg["val_ratio"]),
        neg_train_ratio=int(dcfg["neg_train_ratio"]), neg_eval_ratio=5,
        seed=int(dcfg["seed"]))
    for tier in evidence:
        if tier not in splits_per_tier:
            print(f"[split] {tier}: 0 positives, skipping", flush=True)
            continue
        splits = splits_per_tier[tier]
        torch.save(splits, splits_dir / f"{cell_type}_{tier}_splits.pt")
        print(f"[split] {tier}: train_pos={splits['train']['pos'].shape[1]} "
              f"val_pos={splits['val']['pos'].shape[1]} test_pos={splits['test']['pos'].shape[1]}",
              flush=True)
    leak_report = gb.verify_no_cross_tier_leakage(splits_per_tier)

    # 9. summary -----------------------------------------------------------
    summary = {
        "cell_type": cell_type, "n_genes": len(universe), "n_tfs": len(tf_indices),
        "signature_dim": int(sig_aligned.shape[1]),
        "atac_nonzero_frac": float((atac_aligned[:, 0] > 0).mean()),
        # openness is now a (n_genes, 4) regulatory-potential descriptor (see
        # plan.md Section 2); "proximal" (col 1, max RP > 0.1) is the
        # meaningful accessibility stat -- a mere "any peak in the +/-100kb
        # window" gate (the old `openness_nonzero_frac`) passed ~90% of
        # genes and barely filtered anything, which is why it's kept here
        # only as a diagnostic contrast, not as the graph's accessibility gate.
        "openness_any_peak_in_window_frac": float((openness[:, 3] > 0).mean()),
        "openness_proximal_frac": float((openness[:, 1] > 0.1).mean()),
        "tf_candidate_edges": int(tf_candidate_edges.shape[1]),
        "motif_edges": int(motif_edges.shape[1]),
        # read by src/training/leak_guard.py; true = label-dependent graphs (leaky)
        "prior_exclude_positives": bool(dcfg.get("prior_exclude_positives", False)),
        "coexpr_edges": int(coexpr_edges.shape[1]),
        "negative_pool": int(neg_pool.shape[1]),
        "evidence_sizes": {t: int(e.shape[1]) for t, e in evidence.items()},
        "nesting": nesting,
        "cross_tier_leak_check": leak_report,
    }
    save_json(summary, out_dir / "summary.json")
    print(f"\n[done] {cell_type}: {summary}", flush=True)


if __name__ == "__main__":
    main()
