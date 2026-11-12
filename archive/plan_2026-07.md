# MEvD-GRN Project Plan — Getting to a Publishable Result

Written 2026-09-11. This file explains, in plain language, where the project
stands, what was broken and is now fixed, and what we are going to build next
to make this work strong enough for a top venue (ICLR/ICML, or a strong
biology journal as a fallback). Keep this file updated as we go — it is the
single place to look for "what is the plan right now."

**Latest status entry: Section 11 (2026-10-08).** It covers the two leak
fixes and their reruns, the PBMC10k and BEAR-GRN benchmarks, and the Ada
lockout. It supersedes Section 3b's "Resolved (2026-09-18)" call and
Section 10 item 1.

---

## 1. Where we are right now

### 1.1 The bugs we already found and fixed

Before this plan started, the model's results looked bad, and the worry was
data leakage. We checked carefully and **did not find leakage** in the parts
people usually worry about (train/val/test splits, the graphs used for
message passing, tier isolation). Instead we found two real, more subtle bugs:

1. **A metric bug.** "Early Precision Ratio" (EPR) was being computed by
   comparing a small test set to the wrong baseline number (the size of the
   whole genome-wide TF-gene space, instead of the actual test set). This
   made EPR numbers look huge and meaningless (as high as 449). Fixed — EPR
   numbers are now small and sensible (roughly 1-5x better than random, which
   is normal).

2. **A training bug that looked like "the model forgot everything."** Our
   training curriculum teaches the model on "localization" evidence first,
   then "perturbation" evidence second. After the second stage, the model's
   score on localization data dropped *below random* (AUROC 0.33, where 0.5 =
   random guessing). That is not normal forgetting — normal forgetting drifts
   *toward* random, not *past* it. We found the real cause: during the second
   training stage, the code was picking "hard negative" examples (edges we
   tell the model are wrong) from the **validation and test set** of the
   first stage, not just its training set. In other words, the model was
   being actively taught that some of the exact edges it would later be
   graded on (in the localization test set) were wrong. Fixed by restricting
   "hard negatives" to only the training portion of the lower tier.

**Effect of both fixes**, confirmed on our laptop GPU and again on the Ada
cluster: the model's zero-shot score on the hardest, most trustworthy tier
("dual-evidence", never trained on directly) went from AUPR 0.558 → **0.870**,
and AUROC from 0.838 → **0.968**. This is a big, real improvement — the model
is much better than we thought.

### 1.2 What we now know does and doesn't help (from ablations, i.e.
turning one piece off at a time and re-measuring)

- **The graph neural network (GNN) message-passing step matters a lot.**
  Turning it off collapses the "dual-evidence" score badly. Keep it.
- **The chromatin accessibility (ATAC) data is currently close to useless**,
  and we now know exactly why (see Section 2) — it's a broken pipeline, not
  a fact about biology.
- **Having two different "prior graphs" (co-expression, and a
  TF-to-candidate-target graph) barely helps over having just the
  co-expression one.** The TF-candidate graph is weak because it is built
  directly from the same broken ATAC signal.
- **Correction (2026-09-11), important**: we initially believed "train on
  both tiers at once" (no forgetting, since there's no tier-switch moment)
  was a strictly better replacement for the two-stage schedule, because it
  beat the old schedule on every measure. That comparison was made using the
  *buggy* code. After fixing the leakage bug in Section 1.1, we re-ran both
  properly, full budget, side by side:

  | | tier-1 score | tier-2 score | hardest tier, never trained on |
  |---|---|---|---|
  | train one tier, then the next | 0.57 | **0.63** | **0.87** |
  | train on both at once | **0.94** | 0.35 | 0.79 |

  So training on both at once genuinely removes the "forgets tier 1" problem
  and wins big on tier 1 — but it does *worse* on tier 2 and on the score
  that matters most for the paper (how well it does on the hardest,
  never-trained-on tier). Our best guess why: tier 1 has about 4.6x more
  positive examples than tier 2, so mixing them into one training pool lets
  tier 1 dominate and drowns out tier 2's signal — whereas training on tier
  2 by itself, after tier 1, gives tier 2 the model's full, undivided
  attention for a while, and that attention is what carries over to the
  hardest tier. **We are keeping the one-after-another schedule as the
  default because of this** — the small amount of "forgetting" it causes is
  a reasonable trade-off, not a bug, now that the real bug (leakage) is
  fixed. "Train on both at once" is still kept around as a comparison, since
  it's still a genuinely different, forgetting-free way to train and might
  win under a different setup later (e.g. once Section 5/6's richer data
  changes how much signal each tier carries).

### 1.3 Cluster access

We now have access to Ada, IIIT-H's compute cluster: nodes with **4 RTX
2080 Ti GPUs (11 GB each) and 40 CPU cores**, much more than the single
8 GB laptop GPU we started on. Code and data now live at
`/share1/vishakkashyap.k/mevd_grn` (permanent storage) and get copied to
local `/scratch` on whichever node we're using (fast, temporary storage) for
actual training runs. This removes compute as a bottleneck — see Section 3
for how we should actually use that extra room.

---

## 2. Why ATAC (chromatin accessibility) isn't helping — and how we fix it

We traced the accessibility data all the way from the raw files to the final
model score, and found a chain of small problems that each throw away
information, stacking up into "basically no signal left":

1. A gene "sees" a peak (an open region of DNA) if that peak is anywhere
   within 100,000 base pairs of the gene's start — with **no regard for how
   close or far** it is within that window. A peak right next to the gene
   counts exactly the same as one 99,999 base pairs away. Real biology says
   nearby peaks matter much more.
2. All the peaks near a gene get squashed down into just two numbers (an
   average and a variance across cells). We measured these two numbers on
   real data: they are 96% correlated with each other. That means they carry
   almost exactly the same information — it's really only *one* number, not
   two.
3. There's a separate "openness" score meant to say how accessible a gene's
   region is. It turns out to be a coarse, rounded version of the same one
   number above (64% correlated with it), and it only takes 77 different
   values across all 22,943 genes — most genes get lumped into the same
   handful of buckets.
4. That openness score gets used as a yes/no switch (open or not) when we
   decide which genes a TF can even possibly reach. But 90% of genes pass
   this "open" test, so the switch barely filters anything.
5. The tiny neural network that reads in this data (`ATACEncoder`) is just
   one layer, so it can't fix any of the above — it just passes the weak
   signal through.
6. At the very end, the same weak number gets used **twice** — once as a
   multiplier and once as an add-on term — which doesn't add new information,
   just repeats the same weak one.

**The fix (done, 2026-09-11)**: rebuild the accessibility pipeline properly.
- Weight each peak by how close it is to the gene (closer = more weight),
  instead of a blunt yes/no window. This is the standard approach used in the
  field (called "regulatory potential", using an exponential distance decay).
- Keep a small set of numbers per gene describing the *shape* of its nearby
  peaks (how strong, how close, how many), not just one collapsed average.
  `openness` is now 4 numbers per gene instead of 1.
- Give the little neural network reading this data a real hidden layer so it
  can actually learn something from the richer input, instead of one bare
  linear pass-through.
- (Real TF motif scanning — does this specific TF's known DNA pattern appear
  in this specific peak — is still open, see the note at the end of this
  section.)

**Result, measured on real K562 data, before vs. after:**

| | before | after |
|---|---|---|
| "openness" score, unique values across 22,943 genes | 77 | 20,694 / 6,875 / 20,205 (its 3 richest columns) |
| correlation between openness and the gene-activity average | 0.64 (redundant) | 0.53 (still related, but genuinely less redundant, and now backed by 3 more columns that aren't just a re-derivation) |

**Then we re-ran the actual "does ATAC help" test** (train the full model vs.
training with ATAC zeroed out, same data, same budget). Before the fix,
these two were nearly identical — that's what told us ATAC was contributing
almost nothing. After the fix:

| | ATAC zeroed out | ATAC on (full model) |
|---|---|---|
| tier 1 score | 0.559 | 0.573 |
| tier 2 score | 0.592 | 0.641 |
| hardest tier, never trained on (zero-shot) | 0.856 | **0.881** |

Every single score is now better with ATAC turned on — a small but real and
*consistent* gain (not huge, but real), where before there was essentially
none. This is a clean result for the paper either way it had gone, and it
happened to go the good way: we found a real bug, fixed it, and it paid off.

**Still open**: the real TF motif scanning step (does this specific TF's
known DNA pattern appear in this specific peak) is not done yet — the fix
above only rebuilt the *accessibility* signal (RP scores replacing
openness), not the *motif* signal. That's the next piece: it should give the
model a genuinely different, second graph (Section 4) instead of one built
from the same accessibility numbers as everything else.

---

## 3. Model size — how big should MEvD-GRN actually be?

**Short answer: bigger than it is now, but not by guessing — by testing, and
only after the input data itself gets richer.**

Right now the model is tiny: about **306,000 parameters**. That is not a
typo — for comparison, a small BERT model has over 100 million. Our model is
small on purpose right now, because the data going *into* it is also very
thin: each gene's RNA information is just 3 numbers (average expression,
variance, how often it's detected), and its ATAC information is just 2
numbers (soon to be more, per Section 2). A bigger neural network reading 2-3
numbers per gene doesn't learn more — it just has more ways to overfit
"noise" in those same few numbers. That's exactly what we saw: turning ATAC
off barely changes anything, because there was almost nothing there to lose.

This changes once we do two things later in this plan:
- **Section 2's fix** turns 2 ATAC numbers per gene into roughly 7 richer
  ones (a proper regulatory-potential summary, plus real motif-based edges).
- **Section 5's foundation-model embeddings** (Workstream 4) would add a
  pretrained embedding per gene — typically 256 to 512 numbers per gene,
  learned from tens of millions of real cells by someone else's huge model,
  which we get to reuse for free (no need to train anything that big
  ourselves).

Once the input per gene grows from "a handful of numbers" to "hundreds of
numbers," a 128-dimensional hidden layer becomes a genuine bottleneck, not a
sensible small model. So the plan is:

- **Do not blindly scale up the model just because Ada has more GPU memory
  than the laptop.** The graph itself (about 23,000 genes) is not "big data,"
  and some of our tiers have very few positive examples to learn from (the
  hardest tier, "dual-evidence," has only ~8,000 positive edges total) — a
  much bigger model risks overfitting those few examples, not learning more
  from them.
- **Do scale the model size together with the richer features**, once
  Section 2 and Section 5's foundation-model embeddings land. A reasonable
  new target is roughly **1-3 million parameters** (up from 306,000) — still
  small by modern standards, but sized to match the richer input, not picked
  out of thin air.
- **Use Ada's 4 GPUs to actually test this instead of guessing.** Run a small,
  proper sweep — a few hidden-layer sizes (128 / 256 / 384) crossed with a
  few depths (2 / 3 GNN layers) — in parallel, one setting per GPU, and pick
  the winner by validation score. This turns "we picked hidden_dim=256"
  from a guess into a measured, defensible choice — exactly the kind of
  rigor a top-tier reviewer wants to see.
- **Also use Ada's extra headroom for things that were previously
  memory-constrained**, like training the scMultiomeGRN comparison model with
  its full original settings (it ran out of memory on the 8 GB laptop and
  needed a smaller setting there) and running several experiments side by
  side instead of one at a time.

**Update (2026-09-12): the sweep ran, here's what it found.** Five configs
(hidden_dim x n_gnn_layers) were trained on Ada, one per GPU, on K562 with
the FM-embedding input (the richer input this section predicted we'd need
before scaling up made a difference):

| config  | params    | hidden | layers | dual AUPR | dual AUROC | dual EP |
|---------|-----------|--------|--------|-----------|------------|---------|
| h128l3  | 562,567   | 128    | 3      | 0.9515    | 0.9889     | 0.8725  |
| h256l2  | 1,516,295 | 256    | 2      | 0.9572    | 0.9903     | 0.8809  |
| h256l3  | 2,042,631 | 256    | 3      | 0.9595    | 0.9908     | 0.8830  |
| h384l2  | 3,257,479 | 384    | 2      | 0.9628    | 0.9916     | 0.8888  |
| h384l3  | 4,440,199 | 384    | 3      | 0.9638    | 0.9917     | 0.8911  |

Bigger is monotonically better across all five, confirming this section's
prediction that a richer input (RP-weighted ATAC + FM embeddings) can
actually make use of more capacity, unlike the original 306K-parameter model.
But the gains flatten out fast: going from h384l2 to h384l3 costs 36% more
parameters for +0.001 dual AUPR — essentially noise-level. **Chosen size for
the paper: hidden_dim=384, n_gnn_layers=2 (3.26M parameters)** — it ties the
largest config tested at a fraction of the parameters, which is also the
cleaner "we tested, then picked the smallest model that's within noise of
the best" story for reviewers, rather than just reporting the single biggest
number.

### 3b. Closing the last gap vs. scMultiomeGRN (2026-09-18)

The scMultiomeGRN baseline finally finished a full 3-tier run (see
`results.md` Section 7): it beats our small, no-FM base model on
localization and perturbation, though not on the harder zero-shot
dual-evidence tier. Re-scoring against h384/l2 + Geneformer instead of the
base model closes most of that gap immediately — MEvD-GRN already wins
perturbation (0.834 vs. 0.655 AUPR) and dual-evidence (0.963 vs. 0.847 AUPR)
outright, and only trails on localization (0.751 vs. 0.906 AUPR).

That remaining localization gap has a specific, already-understood cause,
not a sign the architecture doesn't work: Section 3's earlier
sequential-vs-all_at_once comparison already showed `all_at_once` alone (no
FM, small model) pushes localization to 0.937 AUPR — above scMultiomeGRN's
0.906 — because it removes the Stage-2 fine-tuning step that damages what
Stage 1 learned. The catch, at that size, was perturbation dropping to
0.351. The open question is whether that tradeoff still holds once the
model has FM embeddings and 10x the capacity to work with: a bigger,
richer model may have enough room to hold onto localization-relevant
structure while still fitting perturbation well, where the tiny 307K-param
model didn't.

**Resolved (2026-09-18): the tradeoff does not hold at this size.**
`all_at_once` + FM + h384/l2 beats scMultiomeGRN on AUPR across all three
tiers — localization 0.967 vs. 0.906, perturbation 0.690 vs. 0.655, and
zero-shot dual-evidence 0.950 vs. 0.847 (AUROC: wins 2/3, within 0.011 on
perturbation). The perturbation-dilution mechanism from Section 3 was real
at 307K params; a 3.26M-param model with FM embeddings has enough capacity
to fit both tiers without one crowding out the other, so joint training no
longer costs anything on perturbation while still avoiding the sequential
curriculum's localization damage.

**`all_at_once` + FM + h384/l2 is now the recommended configuration for
the paper's headline numbers**, replacing sequential as the default
protocol at this model size (Section 3's "sequential is the default" call
was correct for the small model it was measured on, not a general result).
Practical next step: update `configs/k562.yaml` / a new default config to
this protocol + size once the rest of the pipeline (motif graph, other
cell types) is re-validated against it, so the paper reports one
consistent configuration rather than mixing sequential and all_at_once
numbers across sections.

> **Correction (2026-10-08): the "Resolved" call above was a misread
> comparison.** It compared `all_at_once` + FM + h384/l2 with scMultiomeGRN,
> not with the sequential curriculum at the same size. Against sequential +
> FM + h384/l2 (localization 0.751, perturbation 0.834, dual 0.963 AUPR),
> `all_at_once` wins only localization and loses perturbation (0.696) and
> dual evidence (0.952), so the dilution tradeoff did hold at this size. All
> of these numbers were also trained with the negative-sampling leak. After
> the fix, `all_at_once` loses to scMultiomeGRN on perturbation and the
> curriculum leads on perturbation and dual evidence
> (`docs/experiments/leakfix_rerun.md`; Section 11.1). The "recommended
> configuration" paragraph just above is superseded until the headline
> decision is made (Section 11.7).

A secondary test, `with_replay` at the same size, does NOT close the gap
(localization 0.890, still below both `all_at_once` and scMultiomeGRN) and
its dual-evidence number trains directly on dual_evidence (not zero-shot),
so it isn't a candidate for the headline config. `fm_only` at this size
ties `with_fm` almost exactly (0.963 vs 0.963 dual AUPR) — a secondary
finding that hand-crafted RNA features add ~nothing once Geneformer + a
big-enough model are both present, not the primary question this section
was testing.

---

## 4. The two "prior graphs" — are they really adding anything?

Right now, the model can message-pass over up to two hand-built graphs: genes
that are co-expressed with each other, and TFs linked to genes they might
regulate (picked by accessibility and co-expression together). Testing shows
the co-expression graph does almost all of the work, and the TF-candidate
graph adds very little — which makes sense, since it is built from the same
broken accessibility signal described in Section 2.

Two things will fix this:
1. **Once Section 2's real, motif-based TF→gene graph is built** (based on
   actual DNA binding evidence, not accessibility guesswork — this part is
   still open, see Section 2's closing note), we'll have a graph that might
   carry real, different information from the co-expression graph — worth
   testing properly instead of assuming it's weak forever.
2. **Done (2026-09-11): the model now learns how much to trust each graph**,
   instead of us adding their signals together with equal, fixed weight. We
   added a small, optional switch (`combine_mode: "gated"`) that gives the
   model one learned number per graph, per layer, and lets it decide the mix
   itself — starting from "trust both equally" and adjusting from there
   during training. With only the two current graphs this doesn't move the
   accuracy numbers much yet (as expected — there isn't a genuinely
   different second signal to lean on until the motif graph above exists),
   but the *mechanism* is built, tested, and already saves the learned
   numbers to a results file — so as soon as the motif graph lands, we get
   the "the model learned to rely mostly on co-expression, and only a little
   on the motif graph, and here's the number that proves it" result for
   free, instead of having to build this afterward.

---

## 5. Using pretrained foundation models

The user specifically asked: can we use foundation model embeddings instead
of (or alongside) our hand-built features? Yes, and recent papers (from
2025-2026) specifically show this helps a model generalize to cell types it
has never seen — which lines up with the multi-cell-type work in Section 6.

Plan: use **Geneformer**, a model trained on tens of millions of real
single cells, to get a ready-made "embedding" (a list of a few hundred
numbers) for each gene in our data, with no training of our own required —
just look the gene up and use the numbers it comes with. We add this as a
new, optional input channel alongside our existing RNA/ATAC features, and
test whether it helps:
- on our normal K562 results (within one cell type), and, more importantly,
- on transferring to a cell type the model has never trained on (Section 6) —
  since a foundation model's embeddings don't change from one cell type to
  another the way our hand-built, per-cell-type statistics do, this is where
  we most expect it to help.

**Result (done, 2026-09-11) — this is the single biggest improvement found
in this whole project so far.** We downloaded Geneformer's pretrained gene
embedding table (no training needed, just a lookup: 768 numbers per gene,
matched 17,668 of our 22,943 genes, ~77% coverage), added it as a new input
channel, and re-ran the "does it help" test on real K562 data, full training
budget:

| | hardest tier, never trained on (zero-shot) |
|---|---|
| our hand-built features only (Section 2's fixed-ATAC baseline) | 0.881 |
| Geneformer embedding ALONE (hand-built RNA features zeroed out) | **0.951** |
| Geneformer + hand-built features together | 0.949 |

Two things stand out. First, the size of the jump — 0.881 to 0.951 is a much
bigger gain than anything else we've tried today, including the ATAC fix.
Second, and more surprising: the embedding *by itself* does slightly BETTER
than combining it with our hand-built features. In other words, once the
pretrained embedding is available, our own carefully-built RNA
mean/variance/co-expression numbers stop adding anything — they're now the
weak link, the same way the old ATAC pipeline used to be. This is worth
digging into later (Section 2's fix pattern — richer per-peak features
instead of collapsed averages — may apply to the RNA side too), but for now
the practical takeaway is: **the pretrained embedding is doing most of the
real work, and should be a first-class default input, not an optional
add-on.**

We ran this on K562 only so far (within one cell type). The more important
test — since a pretrained embedding doesn't change from one cell type to
another, unlike our hand-built per-cell-type statistics — is whether it
also makes transfer to a NEW, never-trained-on cell type better. That's
Section 6, next.

We have not yet tried a second foundation model (e.g. scGPT, which is more
citation-heavy in this exact space but harder to set up) — given how well
Geneformer already works, this is a lower priority than finishing the
cell-type-transfer test above.

---

## 6. Testing on more than one cell type

Right now the whole project only trains and tests on one cell type, K562.
The dataset we're using (SC-MO-GRN-DB) actually has several other cell types
with both RNA and ATAC data available, which — as far as we can tell — no
one else has used for this kind of study yet. Showing the model works on
more than one cell type, or better yet that it *transfers* from one cell type
to another it never saw, is one of the strongest things we can show a
reviewer.

Two cell types are close to ready:
- **MCF7** (a breast cancer cell line): has the richest reference network
  after K562 (262 TFs, 1.8 million known edges) and matching RNA+ATAC data is
  already downloaded — though the downloaded file is currently corrupted and
  needs to be repaired or re-downloaded first.
- **Macrophage**: smaller reference network, but its RNA+ATAC data is
  already unpacked and ready to use right now, no fixing needed.

**Macrophage result (done, 2026-09-11) — a genuine surprise that complicates
Section 5's story.** We preprocessed Macrophage (16,202 genes, 22 TFs,
112,383 localization edges) and ran the zero-shot transfer test for real —
train on K562, test on Macrophage with NO retraining at all, using the
transfer script that existed but had never actually been run before today.
For context, we also trained a model from scratch directly on Macrophage's
own data, as a reference point for "how good can a model get on this cell
type at all":

| | AUPR | AUROC |
|---|---|---|
| trained from scratch ON Macrophage (reference ceiling) | 0.841 | 0.915 |
| K562 → Macrophage, zero-shot, WITHOUT the FM embedding | 0.607 | 0.692 |
| K562 → Macrophage, zero-shot, WITH the FM embedding | **0.447** | **0.648** |

This is the opposite of what we expected going in. The pretrained embedding
that gave the single biggest improvement of the whole project on
*within*-K562 generalization (Section 5) actually makes transfer to a
genuinely new cell type WORSE, not better — and by a large margin (0.607 →
0.447 AUPR). We double-checked this isn't a loading bug (the checkpoint
loader uses strict mode, which would hard-error on any architecture
mismatch, and it didn't).

Our best working explanation: the small neural network that reads the
pretrained embedding (`FMEncoder`) gets trained ALONGSIDE the rest of the
model on K562 only, so it — and everything downstream of it — calibrates
itself to K562-specific patterns in a way that doesn't carry over. The
pretrained embedding itself is the same for a gene regardless of cell type,
but *what our model learns to do with it* is not. Our older, simpler,
hand-built features may be "worse but more generic," while the FM-boosted
model is "much better but more of a K562 specialist." This is a real,
interesting finding in its own right, not a failure — it means the paper's
story should be: **foundation-model embeddings are a big win in-domain, but
naive fine-tuning on top of them does not automatically transfer, and that
gap is itself worth studying** (e.g. freezing more of the network, or
training on more than one cell type at once, might fix it — see the
still-open MCF7 work below).

**MCF7 result (done, 2026-09-11) — confirms the Macrophage finding was not a
fluke.** The downloaded zip turned out to be genuinely truncated (an
incomplete download, not real corruption) — deleting it and re-downloading
fresh fixed it. This also used up almost all remaining disk space (the
machine hit 100% full mid-unzip); freed ~5GB by deleting already-extracted
zips, an old backup, and unused alternate-format files that our pipeline
never reads. MCF7 preprocessed cleanly (16,730 genes, 250 TFs, 1.2M
localization edges — the richest network after K562), and we ran the exact
same three-way comparison as Macrophage:

| | AUPR | AUROC |
|---|---|---|
| trained from scratch ON MCF7 (reference ceiling) | 0.933 | 0.921 |
| K562 → MCF7, zero-shot, WITHOUT the FM embedding | 0.758 | 0.720 |
| K562 → MCF7, zero-shot, WITH the FM embedding | **0.597** | **0.553** |

Same pattern, even more pronounced: the FM-augmented model transfers worse
(AUROC actually falls close to random, 0.553), while the plain model
degrades much more gracefully from its in-domain ceiling. Seeing the exact
same direction of effect on TWO independent, unrelated cell types (a
myeloid/immune cell and a breast cancer cell line) makes us confident this
is a real, repeatable pattern, not a one-off. (One implementation detail
along the way: `evaluate_split` used to score an entire eval set in one
un-batched pass, which worked fine at K562/Macrophage's scale but ran out of
GPU memory on MCF7's much larger localization set — fixed by batching the
decode step, which is also just a good general robustness fix.)

**Joint training result (done, 2026-09-11) — a big win, with an honest
nuance.** We built the small joint-training addition described above (new
`scripts/12_joint_train.py`: one shared model, one gradient step per cell
type per epoch, each cell type keeping its own graph/features) and trained
on K562 + MCF7 together, then tested zero-shot on Macrophage — a THIRD cell
type, never seen during this training run at all. We ran this both with and
without the FM embedding, to properly separate "does joint training help"
from "does joint training specifically fix the FM problem":

| | AUPR (zero-shot on Macrophage) | AUROC |
|---|---|---|
| K562 alone, no FM → Macrophage | 0.607 | 0.692 |
| K562 alone + FM → Macrophage (the problem from before) | 0.447 | 0.648 |
| K562 + MCF7 joint, no FM → Macrophage | **0.743** | **0.835** |
| K562 + MCF7 joint + FM → Macrophage | 0.703 | 0.782 |

Two honest findings here, not one: (1) **joint training on more than one
cell type is a big, unambiguous win for transfer, with or without the FM
embedding** (no-FM: 0.607 → 0.743; with-FM: 0.447 → 0.703) — probably the
single most reliable lever we've found for the transfer problem specifically.
(2) Joint training substantially shrinks the FM transfer penalty (the gap
between "with FM" and "without FM" shrinks from -0.160 AUPR at single-cell-type
to only -0.040 AUPR once jointly trained) but does **not fully reverse
it** — in this run, no-FM still edges out FM on the specific zero-shot
transfer number. So we should NOT claim "joint training fixes FM and makes
it best-of-both" — the more defensible, and still very publishable, claim
is "joint multi-cell-type training is the strongest lever for
cross-cell-type transfer we've found, and it makes the FM-embedding penalty
on transfer much smaller, even if not fully gone." In-domain performance on
both training cell types also improved under joint training regardless of
FM (e.g. no-FM joint K562 test AUPR 0.940 vs. 0.881 single-cell-type).

**This is still the strongest, most complete result in the project** — a
model trained jointly across more than one cell type (something nobody else
appears to have done with this dataset) generalizes zero-shot to a third,
completely unseen cell type close to the ceiling of what a model trained
directly on that cell type can achieve (0.743 vs. 0.841 from-scratch, no
FM). This should be the headline result of the paper, with the FM-embedding
interaction reported as an honest, interesting secondary finding rather
than folded into the headline number.

**Robustness check (done, 2026-09-12) — held out all three cell types in
turn, not just Macrophage.** Good news and a nuance:
- **The headline finding holds for every holdout choice**: joint training
  beats single-cell-type transfer every time we have a baseline to compare
  against (Macrophage held out: 0.743 vs 0.607; MCF7 held out: 0.837 vs
  0.758) — this is now a well-replicated result, not a one-off.
- **The "FM slightly hurts" finding does NOT hold for every holdout choice.**
  Holding out K562 instead (train on Macrophage+MCF7), FM actually helps
  (0.744 vs 0.704) — the opposite direction from the Macrophage/MCF7-held-out
  cases. All six with/without-FM numbers across the three holdout choices sit
  within a tight 0.70-0.84 band, so this is a small, second-order,
  holdout-dependent wobble on top of the much bigger, consistent
  joint-training win — worth reporting honestly as "mixed," not smoothing it
  into either "FM helps" or "FM hurts" as a fixed rule under joint training.
  See `results.md` Section 6 for the full 6-row table.

**Still open**: investigating WHY the FM-transfer interaction flips sign by
holdout choice (e.g. does freezing more of the FM-reading network, or a
third training cell type, stabilize it either way).

---

## 7. Order of work

We do the cheap, safe, already-proven fix first, then the harder, riskier
pieces, and pilot the two open-ended ideas (Sections 5 and 6) in parallel
before committing fully to either:

1. ~~Fix the training schedule~~ — **done, and corrected.** We fixed the
   real leakage bug, then properly re-tested "train on both tiers at once"
   against the original one-after-another schedule. Once the leakage bug
   was gone, one-after-another turned out to still be the better choice for
   the score that matters most (see Section 1.2's correction). We're keeping
   it as the default. "Train on both at once" stays available as a
   comparison option, not the default.
2. ~~Fix the ATAC/accessibility pipeline~~ — **mostly done.** The RP
   (distance-weighted) rewrite is done and *measured* to help: with ATAC
   turned off the model now scores clearly worse on every tier, where before
   the fix there was almost no difference at all (Section 2's before/after
   table). Real TF motif scanning (the other half of this workstream) is
   still open.
3. ~~Let the model learn how much to trust each prior graph~~ — **done**
   (the `combine_mode: "gated"` switch, Section 4). Doesn't move the numbers
   yet since there's only one genuinely-different graph so far; will matter
   once the motif graph above exists.
4. **Foundation model embeddings** (Section 5) and **multi-cell-type
   testing** (Section 6) — next up. Pilot both for about a week each, keep
   whichever (or both) shows a real result, drop or shrink whichever doesn't.
5. Update the model size (Section 3) once Sections 2 and 5 give it richer
   input to work with, and prove the new size choice with a small sweep on
   Ada rather than guessing.

Along the way, we rerun the full set of ablation tests after every change
(using the test harness that already exists, `scripts/06_ablation.py`), so
we always know exactly what helped and what didn't, and we keep writing
results into `results.md` and the paper draft (`paper/main.tex`) honestly —
including when something *doesn't* work, since that is still a useful,
publishable finding, and it's the style this project has already been
written in.

---

## 8. What "done" looks like

By the end of this plan, we should be able to say, with real numbers behind
each claim:
- The model's headline result (zero-shot score on the hardest tier) is
  correct and reproducible (already true today: AUPR 0.87).
- We know exactly why chromatin accessibility did or didn't help, backed by
  a clear before/after comparison, not a guess.
- The model learns, and can show, how much it trusts each of its input
  graphs, instead of us picking one by hand.
- We tested whether a large pretrained model's gene embeddings help, and
  have a clear yes/no answer with numbers.
- We tested the model on at least one more cell type, ideally showing it
  transfers to a cell type it never trained on.
- The model's size was chosen by testing a few options, not guessed.
- Every one of these claims has an ablation test behind it in
  `results/ablations/`, and is written up plainly in `results.md` and the
  paper draft.

---

## 9. Independent critical review, and fixing what it found (2026-09-18)

Two fresh-context subagents (mine and a concurrent teammate session working
on the same repo) independently reviewed the codebase for publication
readiness, playing skeptical reviewer against both an ICLR/ICML-caliber
methods bar and a PLOS One/Bioinformatics-caliber journal bar. Both
converged on nearly the same findings (see `docs/critical_review_independent.md`
for the full writeup), which is itself a useful signal that
the findings are real rather than an artifact of one review's framing.
Headline verdict from both: solid, honestly-reported empirical work, not
publication-ready as-is, better suited to a bioinformatics venue than a
top-tier ML conference given the architecture recombines known components
rather than introducing a new one.

**Fixed as a direct result of this review:**
- **A new, previously-undiscovered bug**: MEvD-GRN's own training scripts
  (`03_train.py`, `06_ablation.py`) never seeded model weight
  initialization anywhere -- only the scMultiomeGRN DDP script did. Every
  "same seed 42" rerun was silently using different random initial
  weights. Fixed by adding `torch.manual_seed`/`np.random.seed`/
  `random.seed` before model construction in both scripts, plus a
  `--seed` CLI override for multi-seed sweeps. Verified locally (same seed
  -> identical init, different seed -> different init).
- **The `all_at_once` split-safety issue** (Section 3b above): fixed by a
  teammate session to build the merged train/val/test split from the union
  of each tier's already-persisted, leak-checked splits, with an explicit
  runtime assertion, instead of re-splitting raw evidence independently
  (which was leak-free only by seed/nesting/sort-order coincidence).
  Verified the merged train-positive count is unchanged (964,360), so this
  does not invalidate the numbers already in Section 7 -- it makes the
  guarantee behind them real instead of coincidental.
- **`results/` was entirely gitignored** -- nothing under it had ever been
  tracked in git, so every number in this file, `results.md`, and the paper
  citing a `results/ablations`/`results/baselines` JSON was unreproducible
  from a fresh clone regardless of what was on any one machine's disk.
  Fixed by carving out `.gitignore` exceptions for the small JSON result
  files (checkpoints/logs stay ignored). While syncing, discovered several
  results already written up in this file had never actually been
  transferred from Ada to any local disk at all -- not just untracked, but
  missing outright (the full model-size sweep, most of the multi-cell-type
  transfer runs). All synced and committed now.
- **`paper/main.tex`'s entire body was pre-bugfix** and directly
  contradicted its own Abstract (dual AUPR 0.558 in one table, 0.881 in the
  Abstract, for the same claim). Full rewrite pass: current numbers
  throughout, EPR description fixed to match the actual (already-corrected)
  code, baseline-fairness caveats added explicitly, Limitations expanded
  from 4 to 8 items to cover everything both reviews raised. Compiles
  cleanly.

**Running now, not yet landed:** a 5-seed rerun of the recommended headline
configuration (`all_at_once` + FM + h384/l2), to finally answer whether the
close comparisons throughout this document (e.g. the 0.011 AUROC gap vs.
scMultiomeGRN on perturbation) are real or within run-to-run noise; plus
the `motif_graph_only`/`gated_relations` ablations rerun with the real
28,204-edge motif graph for the first time, and 2 additional scMultiomeGRN
seeds for the baseline side of the same question. Results will replace the
single-run point estimates in `results.md` Section 7 and the paper once
they land.

> **Update (2026-09-18, later the same day; recorded here 2026-09-29):**
> the 5-seed rerun (seeds 42-46) landed. It shows std ≤ 0.0017 on every
> tier and metric. It confirmed that MEvD-GRN beats scMultiomeGRN on 5 of 6
> metrics by margins far above the measured noise. It also showed that
> perturbation AUROC is a small, real loss (0.9096±0.0004 vs. 0.919), not a
> near-tie. Full numbers are in `results.md` Section 7, and the paper was
> updated with them. `motif_graph_only` with the real motif graph also
> landed (`results.md` Section 2). `gated_relations`, however, was
> accidentally rerun without `use_motif: true`, so the 3-relation question
> is still open (Section 10). The 2 extra scMultiomeGRN seeds were still
> running at the last update, so the baseline side of the comparison is
> still single-seed.
>
> **Correction (2026-10-08):** the "5 of 6 metrics" result above compared a
> leaky MEvD-GRN with a non-leaky scMultiomeGRN. With the negative-sampling
> fix, `all_at_once` wins 4 of 6 and loses perturbation on both AUPR and
> AUROC (Section 11.1).
>
> A second, separately spawned review wrote its own report at
> `.claude/worktrees/agent-ac4cbf373d32c1f24/critical_review.md`. That
> worktree was never merged into main. It mostly repeats
> `docs/critical_review_independent.md`. Its one extra point is the
> dual-evidence-as-selection-metric concern in the next paragraph.

**Confirmed true, not yet fixed (a judgment call, not just an edit):**
dual-evidence AUPR has been the deciding metric for roughly six major
project decisions (curriculum protocol chosen and reversed, ATAC fix
validated, model size chosen, FM adopted, joint-training recipe chosen) --
not literal leakage, but the multiple-comparisons pattern where the
supposedly-untouched generalization metric was also the tuning signal.
Flagged explicitly as Limitation (2) in the paper. The real fix (a truly
untouched final-check split, never used for any prior decision) is a
bigger methodological call than a same-session edit, and is left as an
open recommendation rather than something unilaterally implemented.

---

## 10. Open items, deferred work, and working notes (2026-09-18 snapshot, merged 2026-09-29, updated 2026-10-08)

This section takes in the unique content of the former local-only
`status.md` TODO tracker (last updated 2026-09-18, removed in the
2026-09-29 docs cleanup), so open items now live in this one tracked file.
Where that tracker had since been overtaken by later work, this section
says so. The 2026-10-08 update annotates items 1-5 and adds items 6-12;
Section 11 has the background.

**Still open, in priority order:**

1. **Make the recommended configuration the shipped default.** The
   recommended configuration is `all_at_once` + FM + h384/l2 (Section 3b).
   `configs/default.yaml` still ships `protocol: sequential`,
   `hidden_dim: 128`, `use_fm: false` (checked 2026-09-29). To do:
   (a) set `curriculum.protocol: all_at_once`, `model.use_fm: true`,
   `model.hidden_dim: 384` (with `n_gnn_layers: 2`) for K562;
   (b) rerun `results.md` Section 1/2's "current architecture" ablation
   table at this config, because those sections still show the old
   307K-param, sequential, no-FM model; (c) re-validate Macrophage/MCF7
   transfer and joint training at the new config, since all of them were
   tuned and tested at the old size.
   **Update (2026-10-08): (a) is superseded.** The `all_at_once`
   recommendation rested on a misread comparison, and the leak fix widened
   the curriculum's lead (Section 11.1). Do not switch the default to
   `all_at_once` before the headline decision (Section 11.7). If the
   curriculum is chosen, `configs/default.yaml` already ships
   `protocol: sequential`, but a shipped config that reproduces the headline
   still needs `use_fm: true`, `hidden_dim: 384`, `n_gnn_layers: 2`
   (`docs/mevd_vs_scmultiomegrn.md` §8, D8). (b) and (c) still stand, and
   both must now be run with both leak fixes on.
2. **Motif graph + gated relation combiner.** The motif scan finished on
   2026-09-12 (noticed 2026-09-18). It took 35 minutes once a ~20x-oversized
   scan was fixed, and produced 28,204 PWM-backed motif edges, compared with
   871,450 co-expression and 112,500 TF-candidate edges. As of 2026-09-18
   the populated `motif_edges.pt` was on Ada
   (`/share1/.../data/processed/K562/`) and had not yet been synced back to
   the local repo. `motif_graph_only` has since run (`results.md`
   Section 2). The `gated_relations` rerun was launched without
   `use_motif: true`, so whether the combiner does anything useful with 3
   real relations is still open. The rerun is job
   2701095, queued as of 2026-09-18. When it lands, make the
   relation-weight interpretability figure. `use_motif` is still `false` in
   every shipped config.
   **Update (2026-10-08):** job 2701095 ran on the pre-upgrade Ada and has
   not been checked since. The motif scan is restricted to TF-candidate
   pairs, so the motif relation is a *subset* of the TF-candidate relation,
   not a "genuinely different" third relation (`docs/mevd_vs_scmultiomegrn.md`
   §8, D13).
3. **scMultiomeGRN multi-seed.** 2 more seeds were running as of
   2026-09-18. Until they land, the baseline side of the `results.md`
   Section 7 comparison is still n=1.
   **Update (2026-10-08):** still n=1 in `docs/experiments/leakfix_rerun.md`.
4. **Open question:** why does the Geneformer transfer effect flip sign
   depending on which cell type is held out under joint training
   (Section 6)?
   **Update (2026-10-08):** every transfer and joint-training number is
   pre-fix, n=1 and at h128. Rerun them with both leak fixes before asking
   why the sign flips.
5. **Methodological call still to be made:** add a truly untouched
   final-check split, because dual-evidence AUPR has also been used as the
   selection metric (Section 9).
   **Update (2026-10-08):** dual *val* positives are also part of the
   merged val set used for early stopping (`docs/mevd_vs_scmultiomegrn.md`
   §6, R4), which strengthens the case for this split.
6. **(New 2026-10-08) Finish the label-free graph rerun** (Section 11.2).
   It is running on the laptop. When it lands, its numbers replace the
   `leakfix_rerun.md` numbers everywhere (`results.md` banner, the paper).
7. **(New 2026-10-08) Headline-model decision** (curriculum vs
   `all_at_once`, Section 11.7), then the paper rewrite in
   `docs/paper_revision_plan.md`.
8. **(New 2026-10-08) PBMC10k: make the LINGER comparison like-for-like**
   (Section 11.3). Collect the LINGER re-run (job 2399) and re-score every
   method on LINGER's own target set (`linger_tg`); finish GRNBoost2 (job
   475); re-run SCENIC+ and scTFBridge if feasible; find out why the
   target-disjoint regime fails (AUROC 0.455).
9. **(New 2026-10-08) BEAR-GRN: finish the runs** (Section 11.4): the 5-seed
   K562/Macrophage runs and scoring on Ada, the planned indegree-residual
   control (not built yet), and phase 2 (iPS + mouse), keeping BEAR to at
   most 6 submitted job records.
10. **(New 2026-10-08) Regain Ada access** (Section 11.5) and check every
    job that was queued or running at the lockout.
11. **(New 2026-10-08) Hub and identity controls on the K562 benchmark**
    (`docs/mevd_vs_scmultiomegrn.md` §6, R5-R7): degree-only and gene-ID
    baselines, a TF-disjoint split, a learned gene-ID embedding as an FM
    control, and `rna_only` at h384 + FM with several seeds. PBMC10k and
    BEAR-GRN already have trivial baselines; K562 does not.
12. **(New 2026-10-08) Text fixes outside the paper.** `docs/citations.md`
    still describes the TF-candidate graph as leak-free (§8, D4), motivates
    the gated combiner with the Relational Graph Transformer instead of HAN
    (§17, D12), calls the motif graph a genuinely different relation (§19,
    D13) and describes replay and hard negatives as on in the headline (§10
    and §12, D14). The RESCAL and Duren DOI fixes (D10, D11) were made on
    2026-10-08.

**Resolved since that tracker's snapshot** (details in Section 9 and
`docs/critical_review_independent.md`): the `paper/main.tex` body was
resynced to current numbers, including its stale editorial note; the
5-seed variance for the headline configuration landed; the `all_at_once`
split-safety coincidence was fixed with an explicit leak assertion, and
the unused `data` parameter was dropped from `build_all_at_once_stage`
and both of its call sites. Other small fixes made the same day for the
2026-09-18 review:
- The stale `src/training/curriculum.py::build_all_at_once_stage`
  docstring, which still claimed "all_at_once beats sequential on every
  metric", now describes the size-dependent finding.
- The `docs/figures/architecture_diagram_v3_1` subtitle said 306,691 base
  params. It now says 315,015, which was verified by instantiating
  `MEvDGRN` and matches `results.md`/`paper/main.tex`. The +Geneformer
  figure of 430,471 was already correct. The PNG was re-rendered.
- The single-seed caveat was added at the top of `results.md`.

**Caveats to keep stating explicitly** (both reviews): MEvD-GRN trains
45 epochs in total (30 + 15) versus scMultiomeGRN's 2000 per tier, so the
comparison is not compute-matched. GRNBoost2, RegDiffusion and GMF-GAE
are unsupervised, RNA-only methods scored on a supervised multi-omic
benchmark, which makes them a floor rather than peers. The role-aware
decoder and the gated relation combiner show no measurable ablation gap,
so frame them as interpretability choices, not accuracy contributions.

**Lower priority / explicitly deferred:** a second foundation model
(scGPT), more cell types beyond Macrophage/MCF7 (GM12878, HepG2), and the
K562 scCRISPR angle.

**Standing practice:** sync result artifacts (JSONs, logs) from Ada to the
local repo *before* writing their numbers into `results.md`, this file, or
the paper, not after. The 2026-09-18 review found cited numbers that had
no matching file in the repo.

**Workstream numbering** (used in code comments and older docs):
Workstream 1 = curriculum protocol (Sections 1.2, 3b); Workstream 2 = ATAC
pipeline (2a = the RP-weighting rewrite, and the motif-scanning half,
Sections 2 and 7); Workstream 3 = gated relation combiner (Section 4);
Workstream 4 = Geneformer embeddings (Section 5); Workstream 5 =
multi-cell-type transfer and joint training (Section 6).

---

## 11. Status update (2026-10-08): two leaks, two external benchmarks, Ada locked

Everything in this section happened on branch `leak-fix-and-benchmarks`.
The weekly write-up for 27 Sept - 3 Oct
(`docs/weekly_updates/27-sept-to-3rd-oct-2026.md`) covers the same ground
in slide form. Numbers below are copied from the experiment docs named in
each subsection.

### 11.1 The negative-sampling leak: fixed, rerun, and it changes the headline

**The bug.** The trainer drew random negatives from the *whole* 1M-pair
negative pool every epoch, and that pool also supplies every tier's val/test
negatives. On K562 the localization stage asks for about 4M negatives per
epoch from a 1M pool, so it trained on every val/test negative, as a label-0
example, every epoch. The baselines never had this leak.

**The fix** (2026-09-30, commit c26c03d): `MEvDTrainer.restrict_negative_pool`
removes all 300,000 val/test negatives from the pool before training
(1,000,000 → 700,000), and `train_stage` refuses to run without it.
`exclude_eval_negatives: false` reproduces the old behaviour.

**The rerun** (Ada job 467, finished 2026-10-02; full detail in
`docs/experiments/leakfix_rerun.md`): K562, Geneformer + h384/l2, 5 seeds
each of `all_at_once` and `full_curriculum` (the sequential curriculum),
plus 2 legacy-pool sanity runs. With the fix off, the rerun reproduces the
paper's seed-42 numbers exactly, to 4 decimals.

| Run (mean of 5 seeds) | loc AUPR | pert AUPR | dual AUPR | dual AUROC |
|---|---|---|---|---|
| `all_at_once`, pre-fix (paper headline) | 0.9676 | 0.6955 | 0.9517 | 0.9880 |
| `all_at_once`, fixed | 0.9544 | 0.5722 | 0.9017 | 0.9781 |
| `full_curriculum`, fixed | 0.7450 | **0.8165** | **0.9562** | **0.9902** |
| scMultiomeGRN (n = 1, never leaky) | 0.9059 | 0.6550 | 0.8466 | 0.9705 |

What it means:
1. The leak inflated the paper's headline: for `all_at_once`, perturbation
   AUPR drops by 0.123 and dual-evidence AUPR by 0.050.
2. After the fix, `all_at_once` loses to scMultiomeGRN on perturbation
   (AUPR 0.572 vs 0.655, AUROC 0.880 vs 0.919). It wins 4 of 6 metrics, not
   5 of 6.
3. After the fix, the curriculum beats `all_at_once` on perturbation and
   dual evidence and beats scMultiomeGRN on both by wide margins. Its open
   weakness is forgetting localization (0.745 vs 0.954).
4. **The curriculum was ahead before the fix too.** At FM + h384/l2 the
   pre-fix sequential run gives perturbation 0.834 and dual 0.963, against
   `all_at_once`'s 0.696 and 0.952. Section 3b's "Resolved (2026-09-18)"
   call and the matching claims in `results.md` Sections 3 and 7 misread
   that comparison (corrected in place, dated 2026-10-08). The leak hit
   `all_at_once` much harder (seed 42: perturbation −0.120 vs −0.017 for
   the curriculum). A plausible mechanism: `all_at_once` turns hard
   negatives off, so all its negatives come from the leaky random pool.

### 11.2 The second leak (label-dependent TF-candidate graph): quantified, rerun running locally

From `docs/experiments/labelfree_graph_rerun.md` §1:
- `scripts/02_preprocess.py` builds the TF-candidate graph (top-500
  proximally accessible, most co-expressed genes per TF) with every known
  positive of every tier and split excluded (`prior_exclude_positives:
  true`).
- On the paper's graph (112,500 edges), **0** val or test positives of any
  tier are graph edges, against **2.8-3.0%** of the val/test negatives. So
  "this pair is an edge of the input graph" implies "label 0" with
  certainty.
- In the label-free graph, positives and negatives are edges at nearly the
  same rate (about 2.3% vs 2.1%).

The fix keeps every other input identical: a new config,
`configs/sweep/k562_fm_h384_l2_labelfree.yaml`
(`prior_exclude_positives: false`), and a graph-only rebuild that copies
every other processed artifact and checks that the splits are
edge-identical (sanity checks in that doc's §3). The rerun (`full_curriculum`
and `all_at_once`, both leak fixes on) is running one run at a time on the
laptop's RTX 4060, because Ada is locked (11.5). **Its numbers are
pending.** When they land, they, not the 11.1 numbers, are the K562 numbers
for the paper.

The PBMC10k and BEAR-GRN pipelines (11.3, 11.4) built their graphs
label-free and kept the negative-pool fix on from the start, so neither leak
affects them.

### 11.3 PBMC10k Multiome vs LINGER: the main external benchmark (interim, provisional)

Chosen after a survey of about 130 papers
(`docs/reference/multiome_grn_benchmark_consensus.md`): 10x
`pbmc_granulocyte_sorted_10k` Multiome, scored with LINGER's own Cistrome
ChIP-seq evaluation code (20 ChIP datasets, 10 TFs in 4 cell types; the
headline is the mean over the 19 datasets LINGER reports). LINGER, KEGNI,
scTFBridge and regX report on it. Full detail:
`docs/experiments/pbmc10k_linger_benchmark.md`.

Design, pre-registered on 2026-10-01 before any result: CollecTRI labels
with all 10 evaluation TFs removed as regulators, TF-disjoint splits,
degree-matched negatives, label-free graphs, model selection on validation
only, 5 seeds.

Interim result (job 474, all 220 units, 2026-10-02 04:00 IST; §10-11 of that
doc), candidate space `expressed`:

| Method | AUROC (19 datasets) | AUPR ratio |
|---|---|---|
| MeVD-GRN `fm_h384` (CollecTRI, TF-disjoint, 5 seeds) | **0.7343 ± 0.0046** | 2.135 |
| LINGER (published, Table S7) | 0.7143 | **2.2526** |
| MeVD-GRN `fm_h384_rnaonly` | 0.7197 ± 0.0301 | 2.178 |
| MeVD-GRN `base` (no FM) | 0.6701 ± 0.0083 | 1.722 |
| MeVD-GRN `base_rnaonly` | 0.5438 ± 0.0111 | 1.337 |
| MeVD-GRN `fm_h384_uniformneg` | 0.5973 ± 0.0136 | 1.698 |
| Trivial baselines (degree, gene-ID, Pearson) | 0.50-0.59 | |

- MeVD-GRN is above LINGER on AUROC and below it on AUPR ratio.
- The trivial baselines stay at 0.50-0.59, so the signal is not just hub
  structure.
- ATAC adds +0.126 AUROC without Geneformer and only +0.015 with it (within
  the RNA-only seed spread).
- The target-disjoint regime fails (AUROC 0.455 / 0.458): the model does not
  generalise to target genes it never saw as targets.
- **Provisional.** LINGER's published numbers were computed on LINGER's own
  candidate genes; ours are on the `expressed` gene universe. The
  like-for-like comparison needs the LINGER re-run on our cells (job 819 was
  OOM-killed at 30 GB and resubmitted as 2399 with 120 GB) and the
  `linger_tg` re-score. GRNBoost2 (job 475) had scored 1 of 19 datasets.
  SCENIC+ and scTFBridge are quoted from their papers, not re-run. None of
  these jobs can be checked while Ada is locked.

### 11.4 BEAR-GRN: the secondary benchmark, and what its metrics reward

BEAR-GRN (Karamveer *et al.*, *Nat Commun*, 18 Sep 2026,
doi 10.1038/s41467-026-77838-w) is from the SC-MO-GRN-DB lab and benchmarks
9 multiome methods on SC-MO-GRN-DB data against ChIP, knockout, union and
intersection ground truths. Full detail: `docs/experiments/bear_grn_benchmark.md`.

- **Built:** BEAR's exact inputs and ground truths, and a Python port of
  its R scoring. Re-scoring the released GRNs reproduces the paper's AUPRC
  exactly in all 10 checked cases; AUROC agrees within the random
  down-sampling noise.
- **Design:** TF-disjoint 5-fold cross-fitting (the design the BEAR authors
  name in their peer-review file), degree-matched negatives, label-free
  graphs, and trivial baselines scored the same way.
- **Artifacts on K562** (§10.1): a constant score for every measured pair
  already gets ChIP AUPRC 0.431, equal to the best published method (LINGER
  0.430). A TF-disjoint target in-degree ranker beats every BEAR method on
  ChIP (AUROC 0.652, AUPRC 0.484). Against the knockout and intersection
  ground truths nothing is above random AUPRC.
- **First MeVD-GRN K562 result** (§10.2; seed 42 only, a direction, not a
  result): ChIP AUROC / AUPRC 0.602 / 0.468, above every BEAR method but
  below the in-degree ranker; Union 0.634 / 0.347, only marginally above it
  (0.626 / 0.344). So the margin over BEAR's methods is mostly the hub
  prior that any TF-disjoint supervised model gets.
- **Pending:** the 5-seed Ada runs and scoring (queued behind the LINGER
  re-run as of 2026-10-02), the indegree-residual control (planned, not
  built), and phase 2 (iPS + 4 mouse embryo sets + naive mESC, cancelled on
  2026-10-01 and to be resubmitted).
- Note for the headline decision: BEAR's L1 protocol was pre-registered as
  `all_at_once` (§6.2 of that doc). Changing it after seeing results would
  break the pre-registration.

### 11.5 Ada lockout (2026-10-08)

- Both gateways (`ada-gw1`, `ada-gw2`) accept our SSH key and then reply
  "Permission denied". The cause is unconfirmed (`docs/reference/ada.md`,
  top section).
- The 1 Oct setup ran conda installs, 20+ GB downloads and polling scripts
  on ada-gw1. That breaks the login-node limit the HPC admins announced on
  1 Oct (1 CPU core and 2 GB of memory per user on the gateways).
- Consequence: no queued or running job (LINGER 2399, GRNBoost2 475, the
  BEAR-GRN chain) can be checked, and nothing on /share1 can be read. The
  label-free rerun moved to the laptop.
- To do: ask the HPC admins to restore access. Once it is back, run
  installs, downloads and polling loops only through `sbatch` or
  `sinteractive`, never on the gateways.

### 11.6 Novelty audit and the paper

- `docs/mevd_vs_scmultiomegrn.md` (2026-10-01) compares MeVD-GRN with
  scMultiomeGRN component by component. The claimable novelty is the
  evidence-tier protocol, the Geneformer transfer finding and multi-cell-type
  joint training; none of it is a new architecture block. It also lists 18
  places where the code and the docs or paper disagree.
- Every number in `paper/main.tex` is pre-fix. The full list of changes it
  needs, with locations, is in `docs/paper_revision_plan.md` (2026-10-08).
  `main.tex` itself has not been edited.
- Fixed on 2026-10-08: `docs/citations.md` now credits the asymmetric
  bilinear score to RESCAL rather than DistMult, and the Duren *et al.* DOI
  is corrected there and in `paper/references.bib`.

### 11.7 Decisions pending (the user's)

1. **Headline model.** `docs/experiments/leakfix_rerun.md` recommends the
   curriculum (option (a)); the alternative is to keep `all_at_once`. Best
   made once the label-free numbers (11.2) are in.
2. **Story.** Lead with PBMC10k vs LINGER, with SC-MO-GRN-DB (K562) and
   BEAR-GRN as supporting results, and report BEAR's hub and coverage
   artifacts openly.
3. **When to rewrite the paper:** after the label-free rerun and the
   `linger_tg` re-score, following `docs/paper_revision_plan.md`.
