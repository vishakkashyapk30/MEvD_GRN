#!/usr/bin/env python
"""Make an ablation copy of the official scMultiomeGRN code (the original copy stays untouched).

  python scripts/23_patch_official_ablation.py --src <root>/official/ScmultiomeGRN-main \
      --dst <root>/official_ablation/ScmultiomeGRN-main

Every patch is switched on by an extra keyword argument that flows through Trainer.add_argparse_args
(parser.set_defaults(**kw)) into GraphDataModule(**config) and model_cls(**config). With none of them given,
the patched code computes exactly what the original does (same layers, same order, same seeds).

  mp_graph = "train" (default) | "none" | "external"
      graph used for message passing (adj_norm, and the edges whose histogram features are computed).
      "none": no TF-TF links at all, self-loops only. "external": a label-free graph read from mp_adj_file
      (TF x TF table like x_graph.txt, e.g. STRING), reindexed to x_node.txt, made symmetric 0/1, no diagonal.
      Training labels and pos_weight are unchanged (they come from the fold's training file as before).
  use_edge_attr = True (default) | False
      False removes the pair-histogram edge features completely (edge_channels = 0, the conv on them is skipped,
      and they are not computed).
  modalities = ["scrna", "atac"] (default) | ["scrna"] | ["atac"]
      which node-feature modalities enter GraFrank. A single modality keeps GraFrank's message passing (the
      official scMultiomeGRNVariant instead drops GraFrank completely, which removes more than one thing).
  gf_feature_file = None (default) | TSV path (TF x d, rows in x_node.txt order, first column = name)
      Geneformer gene embeddings as a third modality next to scRNA and ATAC.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def sub(text: str, old: str, new: str, path: str, every: bool = False) -> str:
    n = text.count(old)
    assert (n >= 1) if every else (n == 1), f"{path}: expected {'a' if every else 'exactly one'} match for:\n{old}"
    return text.replace(old, new)


def patch(dst: Path):
    # ------------------------------------------------------------------ dataset.py
    p = dst / "src/dataset.py"
    t, tail = p.read_text().split("class DeepTFniDataModule", 1)     # patch GraphDataModule only
    tail = "class DeepTFniDataModule" + tail
    t = sub(t, "        self.mask_type = mask_type\n",
            "        self.mask_type = mask_type\n"
            "        # --- ablation switches (scripts/23_patch_official_ablation.py) ---\n"
            "        self.mp_graph = kwargs.get('mp_graph', 'train')\n"
            "        self.mp_adj_file = kwargs.get('mp_adj_file', None)\n"
            "        self.use_edge_attr = bool(kwargs.get('use_edge_attr', True))\n"
            "        self.gf_feature_file = kwargs.get('gf_feature_file', None)\n", str(p))
    t = sub(t, "        adj_norm = preprocess_graph(adj)\n",
            "        adj_norm = preprocess_graph(self.message_passing_adj(adj))\n", str(p))
    t = sub(t, "        edge_attr = self.compute_edge_attr(self.scrna_file, adj_norm)\n"
               "        edge_attr = torch.tensor(edge_attr).float()\n",
            "        if self.use_edge_attr:\n"
            "            edge_attr = self.compute_edge_attr(self.scrna_file, adj_norm)\n"
            "            edge_attr = torch.tensor(edge_attr).float()\n"
            "        else:\n"
            "            edge_attr = torch.zeros((adj_norm._nnz(), 16, 16)).float()   # never read (edge_channels = 0)\n"
            "        extra = {}\n"
            "        if self.gf_feature_file:\n"
            "            gf = pd.read_csv(self.gf_feature_file, sep='\\t', index_col=0)\n"
            "            assert list(gf.index) == list(self.node_names), 'gf_feature_file rows must follow x_node.txt'\n"
            "            extra['gf_feature'] = torch.tensor(gf.values).float()\n", str(p))
    t = sub(t, "                         edge_index=edge_index, label=label, edge_attr=edge_attr)\n",
            "                         edge_index=edge_index, label=label, edge_attr=edge_attr, **extra)\n", str(p))
    t = sub(t, "        self.atac_feat_dim = self.all_data.data['atac_feature'].shape[-1]\n",
            "        self.atac_feat_dim = self.all_data.data['atac_feature'].shape[-1]\n"
            "        self.gf_feat_dim = self.all_data.data['gf_feature'].shape[-1] if 'gf_feature' in self.all_data.data else 0\n",
            str(p))
    t = sub(t, "    def compute_edge_attr(self, scrna_file, adj_norm, bins=16):\n",
            "    def message_passing_adj(self, adj):\n"
            "        if self.mp_graph == 'train':\n"
            "            return adj\n"
            "        if self.mp_graph == 'none':\n"
            "            return np.zeros_like(adj)\n"
            "        assert self.mp_graph == 'external' and self.mp_adj_file, 'mp_graph=external needs mp_adj_file'\n"
            "        ext = pd.read_csv(self.mp_adj_file, sep='\\t', index_col=0)\n"
            "        ext = ext.reindex(index=self.node_names, columns=self.node_names).fillna(0).values\n"
            "        ext = ((ext + ext.T) > 0).astype(int)\n"
            "        np.fill_diagonal(ext, 0)\n"
            "        return ext\n\n"
            "    def compute_edge_attr(self, scrna_file, adj_norm, bins=16):\n", str(p))
    p.write_text(t + tail)

    # ------------------------------------------------------------------ trainer.py
    p = dst / "src/trainer.py"
    t = p.read_text()
    t = sub(t, "        config['scrna_feat_dim'] = data.scrna_feat_dim\n",
            "        config['scrna_feat_dim'] = data.scrna_feat_dim\n"
            "        config['gf_feat_dim'] = getattr(data, 'gf_feat_dim', 0)\n", str(p), every=True)
    p.write_text(t)

    # ------------------------------------------------------------------ models/model.py
    p = dst / "src/models/model.py"
    t = p.read_text()
    t = sub(t, "        self.model = GraFRankModelAE(atac_feat_dim, scrna_feat_dim)\n",
            "        self.model = GraFRankModelAE(atac_feat_dim, scrna_feat_dim,\n"
            "                                     modalities=kwargs.get('modalities', ('scrna', 'atac')),\n"
            "                                     gf_feat_dim=kwargs.get('gf_feat_dim', 0),\n"
            "                                     use_edge_attr=bool(kwargs.get('use_edge_attr', True)))\n", str(p))
    t = sub(t, "        return self.model(atac_feature, scrna_feature, adj_norm, edge_attr)\n",
            "        return self.model(atac_feature, scrna_feature, adj_norm, edge_attr, gf_feature=kwargs.get('gf_feature'))\n",
            str(p))
    p.write_text(t)

    # ------------------------------------------------------------------ models/GraFRank.py
    p = dst / "src/models/GraFRank.py"
    t = p.read_text()
    t = sub(t, "    def __init__(self, atac_feat_dim, scrna_feat_dim, fusion_dim=512, edge_channels=36, num_layers=2, hidden_dim=256, dropout=0):\n"
               "        super(GraFRankModelAE, self).__init__()\n"
               "        self.fusion = GraFrank(atac_feat_dim+scrna_feat_dim, fusion_dim, edge_channels, num_layers, [scrna_feat_dim, atac_feat_dim])\n",
            "    def __init__(self, atac_feat_dim, scrna_feat_dim, fusion_dim=512, edge_channels=36, num_layers=2, hidden_dim=256, dropout=0,\n"
            "                 modalities=('scrna', 'atac'), gf_feat_dim=0, use_edge_attr=True):\n"
            "        super(GraFRankModelAE, self).__init__()\n"
            "        dims = {'scrna': scrna_feat_dim, 'atac': atac_feat_dim, 'gf': gf_feat_dim}\n"
            "        self.modalities = list(modalities) + (['gf'] if gf_feat_dim else [])\n"
            "        if not use_edge_attr:\n"
            "            edge_channels = 0\n"
            "        self.fusion = GraFrank(sum(dims[m] for m in self.modalities), fusion_dim, edge_channels, num_layers,\n"
            "                               [dims[m] for m in self.modalities])\n"
            "        if not use_edge_attr:\n"
            "            for convs in self.fusion.modality_convs:\n"
            "                for conv in convs:\n"
            "                    conv.no_edge = True\n", str(p))
    t = sub(t, "    def forward(self, atac_feature, scrna_feature, adj_norm, edge_attr):\n"
               "        output = self.encode(torch.cat((scrna_feature, atac_feature), dim=1), adj_norm, edge_attr)\n",
            "    def forward(self, atac_feature, scrna_feature, adj_norm, edge_attr, gf_feature=None):\n"
            "        feats = {'scrna': scrna_feature, 'atac': atac_feature, 'gf': gf_feature}\n"
            "        output = self.encode(torch.cat([feats[m] for m in self.modalities], dim=1), adj_norm, edge_attr)\n",
            str(p))
    p.write_text(t)

    # ------------------------------------------------------------------ models/layers.py
    p = dst / "src/models/layers.py"
    t = p.read_text()
    t = sub(t, "        edge_attr = torch.unsqueeze(edge_attr, 1)\n"
               "        a = self.edge_conv(edge_attr)\n"
               "        a = a.reshape(a.shape[0], -1)\n",
            "        if getattr(self, 'no_edge', False):\n"
            "            a = x_j.new_zeros((x_j.shape[0], 0))\n"
            "        else:\n"
            "            edge_attr = torch.unsqueeze(edge_attr, 1)\n"
            "            a = self.edge_conv(edge_attr)\n"
            "            a = a.reshape(a.shape[0], -1)\n", str(p))
    p.write_text(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--dst", required=True)
    a = ap.parse_args()
    src, dst = Path(a.src).expanduser(), Path(a.dst).expanduser()
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "lightning_logs", "*.pyc"))
    patch(dst)
    (dst / "PATCHED_BY_MEVD").write_text("scripts/23_patch_official_ablation.py\n")
    print("patched copy ->", dst)


if __name__ == "__main__":
    main()
