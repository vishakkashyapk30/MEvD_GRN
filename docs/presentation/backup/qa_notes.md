# Backup notes for questions (not part of the slides)

**Why do the K562 numbers differ from the earlier weekly updates?**
The earlier updates (up to 12 Sept) used an evaluation setup in which some test data had
reached training in two places: (1) the pool of negative examples contained the validation and
test negatives; (2) the input network was built by removing known edges, including test edges,
so a pair's absence from it hinted at its label. Both were fixed on 30 Sept and 8 Oct. The slides
show only results from the fixed setup. Effect on the old single-stage headline model, seed 42,
same data: perturbation AUPR 0.70 to 0.59 and zero-shot test 0.95 to 0.90. The step-by-step model
barely moved (0.83 to 0.81 and 0.96 to 0.96). The old claim "wins 5 of 6 metrics vs
scMultiomeGRN" becomes "step-by-step wins 2 of 3 tests".
The comparison chart is in `backup/fig2_leak_effect.png`; details in
`docs/experiments/leakfix_rerun.md` and `labelfree_graph_rerun.md`.

**Is the scMultiomeGRN comparison fair?**
Same data, same splits, same training tiers, and it was never affected by the issues above. But it is
our re-implementation (one run), not their released code, and not their own dataset.

**Why not run scMultiomeGRN's own dataset (fetal lung)?**
Its answer key comes from DNA-motif hits in the same open-DNA data the model reads (circular) and the
data are probably not truly paired. The pipeline is built (`scripts/14-19`) if we decide to run it.

**Is 0.734 vs LINGER's 0.714 a real win?**
Only as a number. A ranking using open-DNA data alone scores 0.814, and a different TF's ranking
scores the same as the right one, so it does not show TF-specific regulation. We have not re-run LINGER.

**Can we say we beat LINGER on BEAR-GRN?**
Only under the leave-one-TF-out design the BEAR-GRN authors propose. Our model uses ChIP answers, which
BEAR-GRN's stated scope excludes. A no-answer-key model is the proposed next step.

**How many runs behind each number?**
K562: 4 runs (seeds 42-45). PBMC: 5 runs per setting, 220 in total. BEAR-GRN: 1 run (K562 only); the
5-run grid on 9 datasets is not done.
