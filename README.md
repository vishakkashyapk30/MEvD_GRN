# MeVD-GRN

MeVD-GRN infers gene regulatory networks from paired single-cell RNA and ATAC data. It combines per-gene RNA and
ATAC summaries, a pretrained Geneformer gene embedding and label-free gene graphs in a two-tower graph network,
and scores regulator pairs with a symmetric or role-aware decoder. Every evaluation is zero-leak: held-out pairs never
enter training, the message-passing graph or negative sampling.

## Benchmarks

| Benchmark | Task | Ground truth | Doc |
|---|---|---|---|
| scMultiomeGRN (Xu et al., NAR 2025): 9 fetal-lung + 3 blood cell types | undirected TF-TF link prediction | motif in an open promoter | `docs/benchmarks/scmultiomegrn.md` |
| BEAR-GRN (Nat Commun 2026) | TF -> target | ChIP-seq, knockouts | `docs/benchmarks/bear_grn.md` |
| PBMC10k Multiome | TF -> target | Cistrome blood ChIP (LINGER style) | `docs/benchmarks/pbmc10k_linger.md` |
| SC-MO-GRN-DB (original pipeline) | TF -> target | tiered evidence | `docs/sc_mo_grn_db_pipeline.md` |

On scMultiomeGRN the repo also provides two evaluations beyond the official folds: a degree-matched test (negatives
matched to the positives' training degrees) and a held-out-TF split (whole TFs hidden; pairs classed C2/C3 as in
Park & Marcotte 2012), both in `src/benchmarks/scmgrn_protocol.py`.

## Branches

| Branch | Contents |
|---|---|
| `main` | The model (`src/`), data building and evaluation for every benchmark, MeVD-GRN training and evaluation, docs, paper |
| `scmgrn-reproduction` | `main` + end-to-end reproduction of scMultiomeGRN with its official code (`scripts/scmgrn/18_scmgrn_official_repro.py`, worker kind `official`) |
| `ablations` | `scmgrn-reproduction` + MeVD-GRN ablation configs and tables, scMultiomeGRN ablations (patched official code, STRING prior), run plans |
| `baselines` | `scmgrn-reproduction` + statistical and simple baselines, link-prediction heuristics, the look-alike probe, DeepTFni and GENELink |
| `archive/*` (branches and tags) | The pre-restructure branches, unchanged |

## Layout

```
src/                 model, data processing, training, evaluation, benchmark protocols (scmgrn_*, bear_*, pbmc_*)
configs/             sc_mo_grn_db/ (original pipeline + sweeps), scmgrn/, bear/, pbmc/
scripts/             sc_mo_grn_db/ (00-13), scmgrn/ (14-32), bear/ (20-24), pbmc/ (25-29), common/
slurm/               Ada batch jobs per benchmark; slurm/scmgrn/kinds/ = unit kinds of the shared GPU work queue
tests/               unit tests (pytest)
docs/                benchmark records, compiled results, references, figures (index: docs/README.md)
paper/               LaTeX manuscript (September 2026 draft; its numbers predate the leak fixes)
archive/             superseded docs, meeting decks, weekly updates, old SC-MO-GRN-DB results
```

## Running

Environment: `pip install -r requirements.txt && pip install -e .` (Python 3.10+, PyTorch, PyTorch Geometric).
Everything heavy runs as Slurm batch jobs on Ada (`docs/reference/ada.md`). For the scMultiomeGRN benchmark:

```bash
sbatch slurm/scmgrn/scmgrn_setup.sh          # environments, official code, Geneformer cache
bash slurm/scmgrn/scmgrn_submit_all.sh       # raw data -> benchmark build -> GPU work queue
sbatch slurm/scmgrn/scmgrn_cpu_job.sh scripts/scmgrn/30_scmgrn_make_tfholdout_splits.py --root '$SHARE'
sbatch --export=ALL,PLAN=slurm/scmgrn/plans/<plan>.plan slurm/scmgrn/scmgrn_work_plan.sh   # 4 GPU workers
python scripts/scmgrn/25_scmgrn_baseline_table.py --root ~/scmgrn --ref mevd_fm_h384_obs_dm --out tables
```
