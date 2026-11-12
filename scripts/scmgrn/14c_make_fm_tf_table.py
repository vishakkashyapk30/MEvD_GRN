#!/usr/bin/env python
"""Precompute Geneformer input-token embeddings for every TF that can be a node of the
scMultiomeGRN graphs (the HOCOMOCO v11 TF names), so Ada (whose outbound internet is ~10 KB/s)
never has to download the 400 MB model. Same lookup as scripts/11 (a static weight lookup, no
forward pass). Output: npz with `names` (UPPERCASE symbols) and `emb` (n, 768) float32; scripts/16
`--steps fm` uses it when <root>/raw/fm_tf_table.npz exists (identical result to scripts/11).

Usage: python scripts/scmgrn/14c_make_fm_tf_table.py --motif_zip <.../HOCOMOCOv11.zip> --out fm_tf_table.npz
"""
import argparse
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_groundtruth as gt   # noqa: E402

spec = importlib.util.spec_from_file_location("fm11", REPO / "scripts" / "sc_mo_grn_db" / "11_extract_fm_embeddings.py")
fm11 = importlib.util.module_from_spec(spec); spec.loader.exec_module(fm11)

ap = argparse.ArgumentParser()
ap.add_argument("--motif_zip", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
with tempfile.TemporaryDirectory() as tmp:
    names = sorted(gt.motif_tf_names(gt.unpack_motifs(a.motif_zip, tmp)))
emb = fm11.extract_embeddings(names)
np.savez(a.out, names=np.array(names), emb=emb)
print(f"{len(names)} TF names -> {a.out} ({emb.nbytes / 1e6:.1f} MB)")
