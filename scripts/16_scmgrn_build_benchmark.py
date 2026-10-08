#!/usr/bin/env python
"""Build scMultiomeGRN's benchmark for one cell type (ground truth, official
node features, official 10-fold splits) plus MeVD-GRN's processed inputs.

Steps (``--steps``, comma-separated, default all but ``fm``):
  gt        peaks >10% cells -> FASTA -> FIMO (HOCOMOCO v11) -> TFBS p<=1e-6 ->
            TF-promoter containment graph; MAESTRO Enhanced RP (ATAC node
            feature); node intersection -> x_graph.txt / x_node.txt /
            x_atac_rp_score_feature.txt (official file formats)
  official  GRNBoost2 scRNA node feature + tmp/scrna.csv (edge histograms),
            i.e. everything the official scMultiomeGRN code needs
  splits    the official 10-fold split files (seed 666)
  verify    re-run the OFFICIAL split_dataset (from the downloaded repo, in a
            subprocess with a stub pytorch_lightning) and check byte-identity
  processed MeVD-GRN inputs over the same node set (features from MeVD-GRN's
            own preprocessing; label-free message-passing graphs)
  fm        Geneformer embeddings for the node set (needs HF cache/internet)

Usage:
  python scripts/16_scmgrn_build_benchmark.py --config configs/scmgrn/lung.yaml \
      --root $SCMGRN_ROOT --cell_type Lymphoid_cells --fimo_bin fimo --n_jobs 10
"""
from __future__ import annotations

import argparse
import filecmp
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import scmgrn_groundtruth as gt                    # noqa: E402
from src.benchmarks import scmgrn_protocol as proto                    # noqa: E402
from src.data import graph_builder as gb                               # noqa: E402
from src.data.preprocessing import (load_mtx_dir, preprocess_scatac,   # noqa: E402
                                    preprocess_scrna)
from src.utils.io import load_config, save_json                        # noqa: E402


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def rp(root, rel):
    return str(Path(root) / rel) if rel is not None else None


def load_their_scrna(scrna_dir: Path, tf_nodes) -> pd.DataFrame:
    """official load_mtx(use_variable_feature=True, additional_genes=TFs):
    raw counts, genes x cells, restricted to var_features U TFs when a
    var_features.tsv exists."""
    X, genes = load_mtx_dir(str(scrna_dir))                    # cells x genes
    cells = pd.read_csv(scrna_dir / "barcodes.tsv", header=None)[0].astype(str).tolist()
    genes = [str(g) for g in genes]
    vf = scrna_dir / "var_features.tsv"
    keep = np.arange(len(genes))
    if vf.exists():
        want = set(pd.read_csv(vf, header=None)[0].astype(str)) | set(tf_nodes)
        keep = np.array([i for i, g in enumerate(genes) if g in want])
    Xk = sp.csc_matrix(X)[:, keep].T.toarray()                  # genes x cells
    df = pd.DataFrame(Xk, index=[genes[i] for i in keep], columns=cells)
    return df[~df.index.duplicated(keep="first")]


def step_gt(cfg, root, ct, gdir, args):
    ct_dir = Path(rp(root, cfg["sep_data"])) / ct
    work = gdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    log(f"loading ATAC {ct_dir/'atac'}")
    Xa, peaks = load_mtx_dir(str(ct_dir / "atac"))              # cells x peaks
    log(f"ATAC: {Xa.shape[0]} cells x {Xa.shape[1]} peaks")
    fp = gt.filter_peaks(Xa, peaks, float(cfg["peak_threshold"]))
    if args.chroms:
        fp = fp[fp["chrom"].isin(args.chroms.split(","))].reset_index(drop=True)
    fp.to_csv(work / "ATAC_peak.bed", sep="\t", header=False, index=False)
    log(f"peaks with rate > {cfg['peak_threshold']}: {len(fp)} "
        f"({(fp['end'] - fp['start']).sum() / 1e6:.2f} Mb)")
    fasta = work / "peak.fasta"
    if not fasta.exists() or args.overwrite:
        n = gt.write_peak_fasta(fp, rp(root, cfg["genome_fasta"]), str(fasta))
        log(f"wrote {n} sequences -> {fasta}")
    motif_dir = gt.unpack_motifs(rp(root, cfg["motif_zip"]), str(Path(root) / "raw" / "motifs"))
    only = args.motif_subset.split(",") if args.motif_subset else None
    fimo_dir = gt.run_fimo(str(fasta), motif_dir, str(work / "fimo-res" / str(cfg["fimo_thresh"])),
                           fimo_bin=args.fimo_bin, thresh=float(cfg["fimo_thresh"]),
                           n_jobs=args.n_jobs, only_motifs=only)
    tfbs = gt.parse_tfbs(fimo_dir, float(cfg["tfbs_thresh"]))
    tfbs.to_csv(work / "atac_tfbs_region_all.txt", sep="\t", header=False, index=False)
    names = gt.motif_tf_names(motif_dir)
    if only:
        names &= set(only)
    prom = gt.load_tf_promoters(rp(root, cfg["promoter_file"]), names)
    adj0 = gt.build_tf_graph(tfbs, prom)
    gt.write_their_adj(adj0, str(work / "adj.txt"))
    log(f"TFBS (p<={cfg['tfbs_thresh']}): {len(tfbs)}; promoter TFs: {prom['gene'].nunique()}; "
        f"graph: {adj0.shape[0]} nodes, {int(adj0.values.sum())} directed entries")

    # ATAC node feature (MAESTRO Enhanced RP, GRCh38 annotation as in the official code)
    from src.benchmarks.scmgrn_features import load_refgenes, maestro_enhanced_rp
    ref = load_refgenes(rp(root, cfg["maestro_refgenes"]))
    Xi8 = sp.csr_matrix(Xa, copy=True)
    Xi8.data = Xi8.data.astype(np.int8).astype(np.float64)       # official h5 writer casts to int8
    rpm = maestro_enhanced_rp(Xi8, peaks, ref, list(adj0.index), decay=float(cfg["maestro_decay"]))
    cells = pd.read_csv(ct_dir / "atac" / "barcodes.tsv", header=None)[0].astype(str).tolist()
    rpm.columns = cells

    # scRNA genes available (official extract_scrna_feature2 keeps TF nodes in data.index)
    _, rgenes = load_mtx_dir(str(ct_dir / "scrna"))
    nodes = gt.cross_nodes(adj0, rpm.index, set(map(str, rgenes)))
    x_adj = adj0.loc[nodes, nodes]
    gt.write_their_adj(x_adj, str(gdir / "x_graph.txt"))
    with open(gdir / "x_node.txt", "w") as f:
        f.write("".join(f"{n}\n" for n in nodes))
    rpm.loc[nodes].to_csv(gdir / "x_atac_rp_score_feature.txt", sep="\t")
    A = x_adj.values - np.diag(np.diag(x_adj.values))
    summ = {"cell_type": ct, "n_atac_cells": int(Xa.shape[0]), "n_peaks": int(Xa.shape[1]),
            "n_filtered_peaks": int(len(fp)), "n_tfbs": int(len(tfbs)),
            "graph_nodes_before_cross": int(adj0.shape[0]), "n_nodes": len(nodes),
            "adj_sum_directed": int(A.sum()), "n_undirected_edges": int(A.sum() // 2),
            "density_offdiag": float(A.sum() / max(len(nodes) * (len(nodes) - 1), 1)),
            "paper_initial_interactions": cfg["cell_types"].get(ct, {}).get("initial"),
            "restricted_chroms": args.chroms, "motif_subset": args.motif_subset}
    save_json(summ, gdir / "build_summary.json")
    log(f"cross graph: {summ['n_nodes']} nodes, {summ['n_undirected_edges']} undirected edges "
        f"(adj sum {summ['adj_sum_directed']}; paper 'Initial' = {summ['paper_initial_interactions']})")


def step_official(cfg, root, ct, gdir, args):
    from src.benchmarks.scmgrn_features import grnboost2_features
    nodes = [l.strip() for l in open(gdir / "x_node.txt") if l.strip()]
    ct_dir = Path(rp(root, cfg["sep_data"])) / ct
    data = load_their_scrna(ct_dir / "scrna", nodes)
    (gdir / "tmp").mkdir(exist_ok=True)
    data.to_csv(gdir / "tmp" / "scrna.csv")
    log(f"scRNA for GRNBoost2: {data.shape[0]} genes x {data.shape[1]} cells; {len(nodes)} TFs")
    net = grnboost2_features(data, nodes, seed=0, n_jobs=args.n_jobs)
    net.loc[nodes].to_csv(gdir / "x_scRNA_grnboost2_feature.txt", sep="\t")
    log(f"GRNBoost2 feature {net.shape} -> x_scRNA_grnboost2_feature.txt")


def step_splits(cfg, root, ct, gdir, sroot):
    s = cfg["split"]
    d = proto.make_splits(str(gdir / "x_graph.txt"), str(sroot), seed=int(s["seed"]),
                          rate=float(s["train_split_rate"]), n_splits=int(s["n_splits"]),
                          rest_all_train=bool(s["rest_all_train"]))
    log(f"splits -> {d}")


_VERIFY = r"""
import sys, types, random, numpy as np, torch
pl = types.ModuleType("pytorch_lightning")
def seed_everything(s):
    random.seed(s); np.random.seed(s); torch.manual_seed(s); torch.cuda.manual_seed_all(s); return s
pl.seed_everything = seed_everything
pl.LightningDataModule = object
sys.modules["pytorch_lightning"] = pl
from src.dataset import GraphDataModule
print(GraphDataModule.split_dataset(sys.argv[1], seed=int(sys.argv[3]), train_split_rate=float(sys.argv[4]),
      n_splits=int(sys.argv[5]), save_dir=sys.argv[2], rest_all_train=True))
"""


def step_verify(cfg, root, ct, gdir, sroot):
    s = cfg["split"]
    official = rp(root, cfg["official_repo"])
    ours = proto.split_dir(str(sroot), int(s["seed"]), int(s["n_splits"]), float(s["train_split_rate"]),
                           bool(s["rest_all_train"]))
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, "-c", _VERIFY, str(gdir / "x_graph.txt"), tmp,
                            str(s["seed"]), str(s["train_split_rate"]), str(s["n_splits"])],
                           cwd=official, env={**os.environ, "PYTHONPATH": os.pathsep.join(filter(None, [official, os.environ.get("PYTHONPATH")]))},
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"official split_dataset failed:\n{r.stderr[-2000:]}")
        theirs = r.stdout.strip().splitlines()[-1]
        a, b = sorted(os.listdir(ours)), sorted(os.listdir(theirs))
        same = a == b and all(filecmp.cmp(os.path.join(ours, f), os.path.join(theirs, f), shallow=False)
                              for f in a)
    save_json({"identical_to_official_split_dataset": bool(same), "n_files": len(a)},
              gdir / "split_verification.json")
    log(f"split verification vs official code: {'IDENTICAL' if same else 'DIFFERENT'} ({len(a)} files)")
    if not same:
        raise SystemExit("split files differ from the official implementation")


def step_processed(cfg, root, ct, gdir, pdir, mcfg_path):
    mcfg = load_config(mcfg_path)
    dcfg = mcfg["data"]
    nodes = [l.strip() for l in open(gdir / "x_node.txt") if l.strip()]
    gene_index = {g.upper(): i for i, g in enumerate(nodes)}
    ct_dir = Path(rp(root, cfg["sep_data"])) / ct
    pdir.mkdir(parents=True, exist_ok=True)
    gtf = gt.write_promoter_tss_gtf(rp(root, cfg["promoter_file"]), str(pdir / "tss_from_promoters.gtf"))
    pcfg = dict(dcfg)
    pcfg["genome"] = cfg.get("genome", "hg19")
    pcfg["atac"] = dict(mcfg.get("atac", {}) or {}, gtf_path=gtf, genome_fasta=None, jaspar_pfm_path=None)

    feats, names, _hvg, sig = preprocess_scrna(str(ct_dir / "scrna"), pcfg)
    row = {g: i for i, g in enumerate(names)}
    n = len(nodes)
    rna = np.zeros((n, feats.shape[1]), np.float32)
    signature = np.zeros((n, sig.shape[1]), np.float32)
    hit = 0
    for g, i in gene_index.items():
        if g in row:
            rna[i], signature[i] = feats[row[g]], sig[row[g]]
            hit += 1
    atac, openness, _ = preprocess_scatac(str(ct_dir / "atac"), list(gene_index), gene_index, pcfg)

    adj, _ = proto.load_adj(str(gdir / "x_graph.txt"))
    src_, dst_ = np.nonzero(adj)
    evidence = torch.tensor(np.stack([src_, dst_]), dtype=torch.long)
    off = ~np.eye(n, dtype=bool) & (adj == 0)
    ns, nd = np.nonzero(off)
    neg_pool = torch.tensor(np.stack([ns, nd]), dtype=torch.long)
    tf_idx = list(range(n))
    coexpr = gb.build_coexpression_graph(signature, k=int(dcfg["coexpr_knn_k"]))
    cand = gb.build_prior_graph(gene_index, tf_idx, atac, rna,
                                top_k=int(dcfg["tf_candidate_topk"]), signatures=signature,
                                openness=openness, evidence=None, exclude_positives=False)

    np.save(pdir / "rna_features_aligned.npy", rna)
    np.save(pdir / "atac_features_aligned.npy", atac)
    np.save(pdir / "rna_signature.npy", signature)
    np.save(pdir / "openness.npy", openness)
    save_json(gene_index, pdir / "gene_index.json")
    save_json(tf_idx, pdir / "tf_indices.json")
    torch.save(coexpr, pdir / "coexpr_edges.pt")
    torch.save(cand, pdir / "tf_candidate_edges.pt")
    torch.save(torch.zeros((2, 0), dtype=torch.long), pdir / "motif_edges.pt")
    torch.save(neg_pool, pdir / "negative_pool.pt")
    tier = mcfg["data"]["evidence_tiers"][0]
    torch.save(evidence, pdir / f"evidence_{tier}.pt")
    summ = {"cell_type": ct, "n_nodes": n, "rna_matched": hit, "tier": tier,
            "evidence_directed": int(evidence.shape[1]), "coexpr_edges": int(coexpr.shape[1]),
            "tf_candidate_edges": int(cand.shape[1]),
            "prior_exclude_positives": False,   # evidence=None above; read by leak_guard
            "atac_nonzero_frac": float((atac[:, 0] > 0).mean()) if n else 0.0,
            "note": "evidence/negative_pool are the FULL ground truth; per-fold training "
                    "uses only the fold's train split (scripts/17)."}
    save_json(summ, pdir / "summary.json")
    log(f"MeVD-GRN processed -> {pdir}: {summ}")


def step_fm(pdir):
    spec = importlib.util.spec_from_file_location("fm11", REPO / "scripts" / "11_extract_fm_embeddings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    gi = json.load(open(pdir / "gene_index.json"))
    names = [None] * len(gi)
    for g, i in gi.items():
        names[i] = g
    np.save(pdir / "fm_gene_embeddings.npy", mod.extract_embeddings(names))
    log(f"FM embeddings -> {pdir/'fm_gene_embeddings.npy'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--model_config", default=str(REPO / "configs/scmgrn/mevd_base.yaml"),
                    help="MeVD-GRN config whose `data` block drives the processed step")
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--steps", default="gt,official,splits,verify,processed")
    ap.add_argument("--fimo_bin", default="fimo")
    ap.add_argument("--n_jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--chroms", default=None, help="smoke test only: restrict peaks to these chromosomes")
    ap.add_argument("--motif_subset", default=None, help="smoke test only: comma list of TF names")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    root, ct = args.root, args.cell_type
    gdir = Path(rp(root, cfg["graph_dir"])) / ct
    sroot = Path(rp(root, cfg["splits_dir"])) / ct
    pdir = Path(rp(root, cfg["processed_dir"])) / ct
    gdir.mkdir(parents=True, exist_ok=True)
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    t0 = time.time()
    for s in steps:
        log(f"===== {ct}: step {s} =====")
        if s == "gt":
            step_gt(cfg, root, ct, gdir, args)
        elif s == "official":
            step_official(cfg, root, ct, gdir, args)
        elif s == "splits":
            step_splits(cfg, root, ct, gdir, sroot)
        elif s == "verify":
            step_verify(cfg, root, ct, gdir, sroot)
        elif s == "processed":
            step_processed(cfg, root, ct, gdir, pdir, args.model_config)
        elif s == "fm":
            step_fm(pdir)
        else:
            raise SystemExit(f"unknown step {s}")
    log(f"done {ct} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
