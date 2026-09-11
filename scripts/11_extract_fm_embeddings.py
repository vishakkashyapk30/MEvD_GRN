#!/usr/bin/env python
"""Extract pretrained Geneformer gene embeddings for a cell type's gene
universe (plan.md Section 5, Workstream 4).

Pure inference, no fine-tuning: Geneformer's input token-embedding layer
already assigns every gene a learned vector, from pretraining on tens of
millions of real single cells. We reuse that vector directly as an
additional per-gene node feature. Unlike our hand-built, per-cell-type
co-expression signature, this embedding is CELL-TYPE-INVARIANT -- it is the
same vector for a given gene regardless of which cell type we're looking at
-- which is exactly the property that should help cross-cell-type transfer
(see plan.md Section 5's literature note on arXiv 2605.08128).

Usage:
  python scripts/11_extract_fm_embeddings.py --cell_type K562
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.preprocessing import canon
from src.utils.io import load_json

GENEFORMER_REPO = "ctheodoris/Geneformer"
GENEFORMER_SUBFOLDER = "Geneformer-V2-104M"
DICT_VERSION = "gc104M"


def _load_geneformer_dicts():
    """(token_dict: Ensembl ID -> token id, name_dict: CANON gene symbol ->
    Ensembl ID). name_dict is re-keyed through `canon()` so it matches the
    UPPERCASE symbols used everywhere else in this codebase."""
    from huggingface_hub import hf_hub_download

    tok_path = hf_hub_download(GENEFORMER_REPO, f"geneformer/token_dictionary_{DICT_VERSION}.pkl")
    name_path = hf_hub_download(GENEFORMER_REPO, f"geneformer/gene_name_id_dict_{DICT_VERSION}.pkl")
    with open(tok_path, "rb") as f:
        token_dict = pickle.load(f)
    with open(name_path, "rb") as f:
        name_dict_raw = pickle.load(f)
    name_dict = {canon(k): v for k, v in name_dict_raw.items()}
    return token_dict, name_dict


def _load_geneformer_embedding_matrix() -> np.ndarray:
    """(vocab_size, hidden_size) float32 input-embedding weight matrix --
    the static, per-gene-token pretrained embedding, extracted without
    running a forward pass (pure weight lookup, no expression data needed)."""
    from huggingface_hub import snapshot_download
    from transformers import AutoModel

    local_dir = snapshot_download(GENEFORMER_REPO,
                                  allow_patterns=[f"{GENEFORMER_SUBFOLDER}/*"])
    model = AutoModel.from_pretrained(f"{local_dir}/{GENEFORMER_SUBFOLDER}")
    return model.get_input_embeddings().weight.detach().numpy().astype(np.float32)


def extract_embeddings(gene_names: list) -> np.ndarray:
    """gene_names: UPPERCASE symbols in universe order (see gene_index.json).
    Returns (n_genes, hidden_size); genes Geneformer doesn't recognize get an
    all-zero row (same graceful-degradation convention used elsewhere in
    this codebase for missing ATAC/TSS data)."""
    token_dict, name_dict = _load_geneformer_dicts()
    emb_matrix = _load_geneformer_embedding_matrix()
    hidden = emb_matrix.shape[1]
    out = np.zeros((len(gene_names), hidden), dtype=np.float32)
    n_hit = 0
    for i, g in enumerate(gene_names):
        ens = name_dict.get(g)
        tok = token_dict.get(ens) if ens is not None else None
        if tok is not None:
            out[i] = emb_matrix[tok]
            n_hit += 1
    print(f"[geneformer] matched {n_hit}/{len(gene_names)} genes "
          f"({100 * n_hit / max(len(gene_names), 1):.1f}%) to a pretrained embedding "
          f"(hidden_size={hidden})", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--processed_dir", default=None,
                    help="defaults to data/processed/<cell_type>")
    args = ap.parse_args()

    d = Path(args.processed_dir or f"data/processed/{args.cell_type}")
    gene_index = load_json(d / "gene_index.json")
    gene_names = [None] * len(gene_index)
    for g, i in gene_index.items():
        gene_names[i] = g

    fm_emb = extract_embeddings(gene_names)
    out_path = d / "fm_gene_embeddings.npy"
    np.save(out_path, fm_emb)
    print(f"[done] saved {fm_emb.shape} -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
