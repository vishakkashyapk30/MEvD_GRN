# Weekly Update: 5 Sep 2026

**GENE REGULATORY NETWORK INFERENCE FROM INTEGRATIVE MULTI-OMICS DATA**  
By Vishak Kashyap K, UG4 CND  
Advisor: Dr Vinod PK


## Slide 1. Title

Ablation story of MEvD-GRN: the numbers look competitive, the design does not hold up


## Slide 2. The surface story (what a results table wants you to believe)

On K562, the shipped model looks strong in the usual comparison slide:

| Method | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---:|---:|---:|
| MEvD-GRN (ours) | 0.420 | **0.575** | **0.558** |
| scMultiomeGRN | **0.831** | 0.423 | 0.384 |
| GRNBoost2 | 0.535 | 0.206 | 0.173 |
| RegDiffusion | 0.546 | 0.185 | 0.165 |
| GMF-GAE | 0.526 | 0.218 | 0.160 |

Reading only this table: we beat RNA-only baselines on perturbation and dual, and we beat scMultiomeGRN on dual with ~13× fewer parameters.

That is the marketing row. The ablations below are the autopsy.


## Slide 3. What the model claims to be

On paper, MEvD-GRN sells four special ideas:

1. Multi-omic: RNA + ATAC, with ATAC as an accessibility gate  
2. Two structural priors: co-expression graph + TF-candidate graph  
3. Curriculum: localization → perturbation, dual held out as the hard exam  
4. Compact GNN link predictor that generalizes to dual-evidence zero-shot  

If those ideas were working, ablations that remove them should hurt in a clean, intended way. Mostly, they do not.


## Slide 4. Ablation 1: remove the GNN

**What we did.** Turn off message passing (`no_gnn`). Features go straight to the decoder.

| | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---:|---:|---:|
| full model | 0.410 | 0.569 | 0.546 |
| no GNN | 0.498 | 0.341 | **0.241** |

**Insight.** Dual collapses. Message passing is doing real work. This is the one ablation that clearly supports part of the design.

**Critical note.** This only proves "some graph smoothing helps." It does **not** prove that our specific two-prior regulatory story is correct. Later ablations show the GNN is mostly riding co-expression.


## Slide 5. Ablation 2: remove ATAC

**What we did.** Zero ATAC features and openness (`rna_only`).

| | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---:|---:|---:|
| full model | 0.410 | 0.569 | 0.546 |
| RNA only | 0.414 | 0.540 | **0.543** |

**Insight.** Dual barely moves (0.546 → 0.543). The multi-omic story is not earning its place. Chromatin is in the diagram; it is not in the performance gap.

**Critical note.** If ATAC can be deleted with almost no cost, then:
- the role-aware "accessibility gate" is mostly theater under current features
- any prior that depends on openness (TF-candidate graph) is built on a weak signal


## Slide 6. Ablation 3: change how RNA and ATAC are fused

**What we did.** Replace role-aware decoding with legacy gated fusion or plain concat.

| | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---:|---:|---:|
| role-aware (default) | 0.410 | 0.569 | 0.546 |
| gated fusion | 0.419 | 0.588 | 0.546 |
| concat fusion | 0.421 | **0.595** | **0.555** |

**Insight.** Fancy fusion does not win. Concat is as good or slightly better. The biologically motivated decoder is not the reason dual works.

**Critical note.** We are paying complexity for an integration story that ablations say is second-order at best.


## Slide 7. Ablation 4: train on one tier only

**What we did.** `loc_only`, `pert_only`, `dual_only`.

| Ablation | Loc AUPR | Pert AUPR | Dual AUPR | Reading |
|---|---:|---:|---:|---|
| loc_only | **0.940** | 0.187 | 0.663 | learns binding landscape; fails functional edges |
| pert_only | 0.407 | 0.496 | 0.485 | never gets the easy warm start |
| dual_only | 0.779 | 0.197 | 0.479 | trains on the gold set; still not magical |

**Insight.** No single tier is enough for the full story. Localization alone transfers somewhat to dual (0.663) but dies on perturbation. Perturbation alone never becomes great. Dual alone is small and does not invent the rest of the network.

**Critical note.** This supports using multiple tiers, but it does **not** prove that our sequential curriculum is the right way to combine them.


## Slide 8. Ablation 5: curriculum vs train together

**What we did.** Compare shipped sequential training to joint training and to a heavier consolidation recipe.

| Setup | Loc AUPR | Pert AUPR | Dual AUPR | Honest label |
|---|---:|---:|---:|---|
| shipped loc → pert | 0.410 | 0.569 | 0.546 | default / paper recipe |
| all_at_once (loc+pert together) | **0.937** | 0.348 | **0.788** | no sequential curriculum |
| with_replay (3-stage + consolidation) | 0.765 | 0.411 | **0.940** | trains on dual too |

**Insight.**
- Sequential curriculum **destroys localization** (AUROC falls below chance after Stage 2).
- Training loc and pert **together** keeps localization and improves zero-shot dual (0.788).
- The flashy 0.940 dual needs Stage 3 on dual itself, so it is not the same claim as zero-shot.

**Critical note.** The method's namesake idea (evidence curriculum) is one of its weakest parts. The default recipe looks like a design failure next to `all_at_once`.


## Slide 9. Ablation 6 (this week): turn off each prior graph

**What we did.** Keep the same training recipe, but let the GNN see only co-expression, only TF-candidate, or both.

| Run | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---:|---:|---:|
| both graphs | 0.404 | **0.507** | **0.451** |
| coexpr only | 0.418 | 0.441 | 0.383 |
| TF-candidate only | 0.418 | 0.374 | **0.258** |

**Insight.**
- Co-expression alone already carries most of the dual score (0.383).
- TF-candidate alone falls to 0.258: essentially `no_gnn` territory.
- "Two priors" is mostly **one prior** with a weak second relation attached.

**Critical note.** This is the cleanest proof that the architecture is overclaiming multi-relational biology. The GNN is a co-expression smoother, not a balanced coexpr + TF-candidate reasoner.


## Slide 10. The contradiction in one sentence

**Numbers say:** MEvD-GRN beats strong baselines on perturbation and dual-evidence.  

**Ablations say:** it does that while (a) barely using ATAC, (b) barely using the TF-candidate prior, (c) using a curriculum that forgets localization, and (d) relying on co-expression neighborhood smoothing more than on the "role-aware multi-omic" story.


## Slide 11. Why ATAC features make this failure predictable

How we build gene activity:

```text
gene_activity = peak_matrix  x  incidence_matrix
```

Peak linked to gene if midpoint is within ±100 kb of the TSS (binary). Then we keep only mean, variance, and openness.

That design is blunt:
- near and far peaks count the same
- no promoter / enhancer distinction
- most cell-level ATAC detail is discarded

So the `rna_only` flatline and the dead TF-candidate graph are not surprises. They are what you get when accessibility enters as a coarse TSS-window sum.


## Slide 12. Full ablation map (quick reference)

| Ablation | Dual AUPR | What it exposed |
|---|---:|---|
| no_gnn | 0.241 | need message passing |
| tf_cand_only | 0.258 | TF-candidate prior is weak alone |
| coexpr_only | 0.383 | co-expression does most GNN work |
| both graphs | 0.451 | second prior adds little |
| shipped curriculum | 0.546 | competitive headline, broken loc |
| rna_only | 0.543 | ATAC not used |
| gated / concat fusion | 0.546 / 0.555 | role-aware decoder not the win |
| loc_only | 0.663 | binding signal transfers some, not enough |
| all_at_once | 0.788 | curriculum order is the wrong default |
| with_replay (trains dual) | 0.940 | different claim; not zero-shot |

Pattern: every "special" module except generic graph smoothing fails its stress test.


## Slide 13. Bottom line

We should not tell the advisor: "ablations confirm the architecture."

We should tell the advisor:

> Head-to-head numbers make MEvD-GRN look better than RNA-only baselines and competitive on dual. Controlled ablations show that advantage does not come from the parts we advertise. ATAC, the TF-candidate prior, role-aware fusion, and sequential curriculum are not carrying the model. Co-expression GNN smoothing is. Until those pieces survive ablation, the design is not validated.

That is the actual picture.


## Slide 14. What to fix (design, not spin)

1. Rebuild ATAC gene activity or drop ATAC claims until it moves metrics.  
2. Redesign or remove the TF-candidate prior.  
3. Replace sequential curriculum as the default (`all_at_once` or real consolidation).  
4. Keep dual zero-shot only when the supporting modules survive ablations.  
5. Stop treating "beats GRNBoost2 on dual" as proof that the architecture is right.


## Slide 15. Thank you

Questions and feedback welcome.
