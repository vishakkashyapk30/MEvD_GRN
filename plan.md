# MEvD-GRN Project Plan — Getting to a Publishable Result

Written 2026-09-11. This file explains, in plain language, where the project
stands, what was broken and is now fixed, and what we are going to build next
to make this work strong enough for a top venue (ICLR/ICML, or a strong
biology journal as a fallback). Keep this file updated as we go — it is the
single place to look for "what is the plan right now."

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

**Natural next steps** (not yet done): hold out MCF7 or K562 instead of
Macrophage, to make sure this pattern isn't specific to which cell type gets
held out; and investigate WHY FM still slightly underperforms even jointly
trained (e.g. does freezing more of the FM-reading network, or adding a
third training cell type, close the remaining gap).

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
