"""Label-free preprocessing of BEAR-GRN's INPUT.DATA into MeVD-GRN's processed
layout (the files src/data/dataset.py:load_celltype_data reads).

Leak-safety (docs/experiments/bear_grn_benchmark.md s6.5):
  * gene universe = the post-QC RNA gene list of INPUT.DATA (no network filter)
  * TF nodes = RNA genes that are a source in any label/GT file of the dataset
    (which genes are TFs -- not which edges exist)
  * co-expression kNN and TF-candidate graphs are built WITHOUT any edge label
    (`build_prior_graph(..., evidence=None, exclude_positives=False)`), so the
    second leak of scripts/02_preprocess.py (candidate graph excluding all
    positives, val/test included) cannot occur
  * the motif relation is empty
Label edges are saved separately (labels/<tier>.pt, gts/<name>.pt, in universe
indices); scripts/22 removes the held-out TFs from them per fold.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch

import contextlib

import scipy.sparse as sp

from src.data import graph_builder as gb
from src.data import preprocessing as _pp
from src.data.preprocessing import preprocess_scatac, preprocess_scrna
from src.utils.io import save_json


def resolve_input(root: Path, rel: str) -> str:
    """Path of an INPUT.DATA member: the extracted file if present, else
    'zip:<root>/zips/INPUT.DATA.zip.sparse::<rel>' (member kept compressed by
    scripts/20 --no_extract, so a 2.5 GB CSV never lands on disk)."""
    p = Path(root) / rel
    if p.exists():
        return str(p)
    sp_zip = Path(root) / "zips" / "INPUT.DATA.zip.sparse"
    if sp_zip.exists():
        import zipfile
        with zipfile.ZipFile(sp_zip) as z:
            if rel in z.namelist():
                return f"zip:{sp_zip}::{rel}"
    return str(p)


def read_csv_sparse(path: str, chunksize: int = 4000):
    """Features x cells CSV (BEAR INPUT.DATA layout: quoted header of barcodes,
    first column = feature name) -> (cells x features CSR float32, feature names),
    parsed in row chunks so a 2.5 GB mouse ATAC matrix never exists densely.
    `path` may be 'zip:<zipfile>::<member>' (streamed decompression)."""
    blocks, names = [], []
    src = path
    zf = None
    if str(path).startswith("zip:"):
        import zipfile
        zpath, member = str(path)[4:].split("::", 1)
        zf = zipfile.ZipFile(zpath)
        src = zf.open(member)
    for ch in pd.read_csv(src, index_col=0, chunksize=chunksize, engine="c"):
        names.extend(str(x) for x in ch.index)
        blocks.append(sp.csr_matrix(ch.to_numpy(dtype=np.float32)))
    X = sp.vstack(blocks).T.tocsr()                     # cells x features
    if zf is not None:
        zf.close()
    return X, names


@contextlib.contextmanager
def chunked_csv_loading():
    """Route src.data.preprocessing's CSV loading through read_csv_sparse (only
    for *.csv paths) without editing the shared module."""
    orig = _pp.load_features_by_cells

    def patched(path):
        if str(path).endswith(".csv"):
            return read_csv_sparse(str(path))
        return orig(path)

    _pp.load_features_by_cells = patched
    try:
        yield
    finally:
        _pp.load_features_by_cells = orig


def read_gt_pairs(path: str) -> pd.DataFrame:
    """BEAR GT file -> distinct upper-cased (Source, Target). 'NA' stays a gene
    name (matches R's read_tsv + paste behaviour, see bear_metrics.read_table)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    # SC-MO-GRN-DB RN_TSV files sometimes separate Target and Relationship with
    # spaces instead of a tab ("AASDHPPT  Transcriptional_Regulation"): keep the
    # first whitespace token (gene symbols contain no spaces).
    tgt = df["Target"].str.strip().str.split(r"\s+", n=1, regex=True).str[0]
    df = pd.DataFrame({"Source": df["Source"].str.upper().str.strip(),
                       "Target": tgt.str.upper()}).dropna()
    return df.drop_duplicates().reset_index(drop=True)


def pairs_to_index(df: pd.DataFrame, gene_index: Dict[str, int]) -> torch.Tensor:
    s = df["Source"].map(gene_index)
    t = df["Target"].map(gene_index)
    ok = s.notna() & t.notna() & (s != t)
    e = np.stack([s[ok].astype(np.int64).to_numpy(), t[ok].astype(np.int64).to_numpy()])
    e = np.unique(e, axis=1) if e.shape[1] else e
    return torch.from_numpy(e.astype(np.int64))


def prepare_dataset(ds_cfg: dict, root: Path, out_dir: Path, pcfg: dict, gtf_path: str | None,
                    genome: str, compendium: str | None = None) -> dict:
    """Build out_dir/{rna_features_aligned,atac_features_aligned,rna_signature,
    openness}.npy, gene_index.json, tf_indices.json, {coexpr,tf_candidate,motif}_edges.pt,
    negative_pool.pt (empty placeholder), labels/*.pt, gts/*.pt, summary.json."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rna_path = resolve_input(root, ds_cfg["rna"])
    atac_path = resolve_input(root, ds_cfg["atac"]) if ds_cfg.get("atac") else None

    # 1. RNA (BEAR's QC'd cells; our normalisation: CP10k + log1p)
    with chunked_csv_loading():
        rna_feats, rna_names, _hvg, rna_sig = preprocess_scrna(rna_path, pcfg)
    universe = list(rna_names)                        # post-QC RNA genes, in file order
    gene_index = {g: i for i, g in enumerate(universe)}

    # 2. TF nodes: sources of any label / GT file (gene-level fact only).
    # Regime L2 (s6.2): the ONLY training labels are the non-specific compendium
    # (RN005/RN006); the cell-type networks are not loaded as labels at all.
    if compendium:
        label_dfs = {"compendium": read_gt_pairs(str(root / compendium))}
    else:
        label_dfs = {t: read_gt_pairs(str(root / p)) for t, p in ds_cfg.get("labels", {}).items()}
    gt_dfs = {}
    for name, p in ds_cfg.get("gts", {}).items():
        fp = root / p
        if fp.exists():
            gt_dfs[name] = read_gt_pairs(str(fp))
    tf_names = set()
    for d in list(label_dfs.values()) + list(gt_dfs.values()):
        tf_names |= set(d["Source"])
    tf_indices = sorted(gene_index[g] for g in tf_names if g in gene_index)

    # 3. ATAC regulatory potential, label-free
    gcfg = dict(pcfg)
    gcfg["genome"] = genome
    gcfg["atac"] = dict(pcfg.get("atac", {}), gtf_path=gtf_path)
    with chunked_csv_loading():
        atac_feats, openness, top_peaks = preprocess_scatac(atac_path, universe, gene_index, gcfg)

    # 4. message-passing graphs: NO labels anywhere
    coexpr = gb.build_coexpression_graph(rna_sig, k=int(pcfg.get("coexpr_knn_k", 20)))
    tf_cand = gb.build_prior_graph(gene_index, tf_indices, atac_feats, rna_feats,
                                   top_k=int(pcfg.get("tf_candidate_topk", 500)),
                                   signatures=rna_sig, openness=openness,
                                   evidence=None, exclude_positives=False)

    np.save(out_dir / "rna_features_aligned.npy", rna_feats.astype(np.float32))
    np.save(out_dir / "atac_features_aligned.npy", atac_feats.astype(np.float32))
    np.save(out_dir / "rna_signature.npy", rna_sig.astype(np.float32))
    np.save(out_dir / "openness.npy", openness.astype(np.float32))
    save_json(gene_index, out_dir / "gene_index.json")
    save_json(tf_indices, out_dir / "tf_indices.json")
    torch.save(coexpr, out_dir / "coexpr_edges.pt")
    torch.save(tf_cand, out_dir / "tf_candidate_edges.pt")
    torch.save(torch.zeros((2, 0), dtype=torch.long), out_dir / "motif_edges.pt")
    torch.save(torch.zeros((2, 0), dtype=torch.long), out_dir / "negative_pool.pt")
    # per-gene top-K RP-weighted peaks (chrom, start, end) for the motif pair features (s12.6 M3)
    save_json({str(g): [list(map(str, p)) for p in pk] for g, pk in top_peaks.items()},
              out_dir / "gene_top_peaks.json")

    (out_dir / "labels").mkdir(exist_ok=True)
    (out_dir / "gts").mkdir(exist_ok=True)
    label_sizes, gt_sizes = {}, {}
    for t, d in label_dfs.items():
        e = pairs_to_index(d, gene_index)
        torch.save(e, out_dir / "labels" / f"{t}.pt")
        label_sizes[t] = {"edges_in_universe": int(e.shape[1]), "edges_file": int(len(d)),
                          "tfs_in_universe": int(len(set(e[0].tolist())))}
    for name, d in gt_dfs.items():
        e = pairs_to_index(d, gene_index)
        torch.save(e, out_dir / "gts" / f"{name}.pt")
        gt_sizes[name] = {"edges_in_universe": int(e.shape[1]), "edges_file": int(len(d)),
                          "tfs_file": int(d["Source"].nunique()),
                          "tfs_in_universe": int(len(set(e[0].tolist()))),
                          "targets_file": int(d["Target"].nunique())}
    summary = {
        "n_genes": len(universe), "n_tf_nodes": len(tf_indices),
        "n_tfs_named_in_labels": len(tf_names),
        "atac_nonzero_frac": float((atac_feats[:, 0] > 0).mean()),
        "openness_proximal_frac": float((openness[:, 1] > 0.1).mean()),
        "coexpr_edges": int(coexpr.shape[1]), "tf_candidate_edges": int(tf_cand.shape[1]),
        "labels": label_sizes, "gts": gt_sizes,
        "label_free_graphs": True,     # build_prior_graph called with evidence=None
        # zero-leak guard fields (read by the shared training-entry guard)
        "prior_exclude_positives": False,
        "motif_exclude_positives": False,
        "prior_exclude_positives": False,   # read by src/training/leak_guard.py
        "regime": "L2" if compendium else "L1", "compendium": compendium,
    }
    save_json(summary, out_dir / "summary.json")
    return summary
