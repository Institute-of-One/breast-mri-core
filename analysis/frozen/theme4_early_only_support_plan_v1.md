# Is the analysis support late-dependent? The rule, fixed before the result (2026-11-11)

Referee 1 point 3 and Referee 2 point 2 of BPEX-107031 ask whether the paired
comparison is biased because both arms are measured inside an analysis support
that the archive built with the full dynamic series, late image included.

This document fixes the rule **before any number from part B or part C
exists**, because every constant in it could otherwise be tuned until the
answer looked favourable. Where the published specification supplies a
constant, the specification's value is used and cited. Where it does not,
**no value is chosen**: the result is reported across a declared grid and the
claim is only what holds over the whole grid.

Frozen by `scripts/theme4_early_only.py --freeze`, which pins the sha256 of
this document and of the analysis script into
`docs/theme4_early_only_lock_v1.json`. The script refuses to run if either has
changed.

---

## 0. What was already known when this was written

**Part A had already been measured.** In an exploratory run on 2026-11-11,
before this document existed, the containment fraction of §2 was computed on
232 examinations and returned a mean of 0.99925, a minimum of 0.983, two
examinations below 0.99, and 0.058 % of 5 650 468 support voxels failing. That
is recorded here as a prior observation, not as a prespecified result, and the
committed script recomputes it on all 234 examinations of the published cohort.

The freeze therefore protects parts B and C, which had produced no comparison
and no overlap statistic when it was taken. It protects part A only in the weak
sense of fixing how the number is computed; **part A has no constant to tune**,
because the enhancement threshold is read per examination from the archive's
own stored label (§2), and the only other choice — what to do with a voxel that
has no precontrast signal — is settled against the result by counting such
voxels as failures.

**The design changed, and an earlier run is why.** The first version of this
plan made part C the primary experiment: redraw the VOI from orthogonal
projections with no reader input, and run the comparison inside it. A smoke run
on 2026-11-11 put the redrawn support on the **contralateral breast** for
participant ACRIN-6698-103939, with a Dice of 0.000 against the archive's
support. That is a property of the rule, not a bug, and it is reported in
part C. But an experiment whose support is sometimes the wrong organ cannot
answer the referee, so part C became the secondary, negative result and the
measurements that do answer the referee became parts A and B. The superseded
design is in §7.

## 1. What the supplied support actually bundles

`theme3_ftv.load_case` takes the eligible set as the **inverse** of the
archive's `VOLSER Analysis Mask` — the zero-valued voxels. The archive's mask
therefore already encodes four things at once, and the manuscript says so in
§2.2:

1. the manually drawn rectangular VOI;
2. the percent-enhancement selection, nominally PE > 70 %;
3. the precontrast background exclusion, nominally 60 % of the within-VOI 95th
   percentile;
4. neighbour-based connectivity, plus any OMIT regions.

`volume()` then applies only the SER threshold inside that set. The referee's
objection is that items 1 to 4 might depend on the late image, in which case
measuring both arms inside that set would favour the arm that has the late
image. §2.2 currently says the support was not reconstructed; **this experiment
reconstructs it**, and §2.2 will be rewritten to say what was found rather than
why it was not tried.

## 2. Part A: containment

**The question.** Does the supplied support contain a single voxel that the
late image could have put there?

**The measurement.** For every examination, for every voxel of the supplied
support, percent enhancement is computed from the precontrast and early images
alone,

> `PE = 100 (early - pre) / pre` where `pre > 0`, undefined elsewhere,

and tested against the threshold at which **that examination** was analysed.

**The threshold is read, not assumed.** The archive stores the threshold in the
label of its own FTVpe measurement, for example `FTVpe(PEearly>70%)`. It is not
always 70 %: in the published cohort six participants carry 40, 50 or 60 %. An
earlier run that assumed 70 % reported a containment of 0.178 on one of them,
which was an error of mine and not a property of the archive. The script parses
the percentage out of the stored label, and asserts that the stored SER label
names the same threshold.

**A voxel with no precontrast signal counts as failing.** It has no percent
enhancement, and the reading that cannot flatter the result is the one that
calls it a failure rather than a pass.

**Two readings are reported.** The specification's reading, `PE > threshold`,
and the stricter one the archive's stored PE maps imply, `floor(PE) >
threshold`, because those maps are floored to whole percent
(`data/derived/theme3_pilot_v1/pe_reference_audit.json`). The specification's
reading is primary.

**Independently of the voxels**, the script checks that the supplied support's
voxel count equals the archive's stored FTVpe count. If it does, the archive
itself is asserting that its support *is* the early percent-enhancement
selection; the voxelwise test then confirms that assertion from the raw images
rather than taking the label's word for it.

**What part A can and cannot show.** If the supplied support is contained in
the early-only selection, then no voxel in it was chosen by late behaviour, and
items 2, 3 and 4 of §1 cannot be carrying late information into the comparison.
It says nothing about item 1, the VOI, which is reader input by design — the
specification has a person draw it on an early postcontrast image or its
subtraction, and part B tests what happens inside it.

## 3. Part B: the paired comparison inside a late-independent support

Inputs: the `pre` and `early` arrays, the voxel geometry, and the supplied
support's bounding box. **The `late` array is not read.** The script asserts
this by holding the arrays behind a wrapper that raises on any access to `late`
while a support is being built; the late image is unlocked afterwards, once,
for the SER measurement, which is the step that is *supposed* to use it.

1. **Box.** The axis-aligned bounding box of the supplied support. This is the
   one thing taken from the archive, and it is the one thing the specification
   says a person supplies. By part A the supplied support is itself an
   early-only selection, so its bounding box carries no information the late
   image could have put there.
2. **Margin.** The box is dilated by `m` voxels in every direction and clipped
   to the array. **`m` is not chosen.** It runs over the declared grid
   `{0, 2, 4, 8}` and every result below is reported at every `m`. A statement
   that holds at some margins and not others is reported as that, not as a
   result.
3. **Enhancement selection.** `PE > threshold`, with the threshold read per
   examination exactly as in part A. *Published threshold [9], per-examination
   value from the archive.*
4. **Background exclusion.** Inside the box, keep voxels with
   `pre > 0.60 * p95(pre within the box)`. *Published rule [9], applied within
   the dilated box rather than the archive's VOI.*

The result is the **late-independent support**. Inside it, with the published
`volume()` function unmodified and the archive's own stored SER threshold and
cap:

* `V_sub` = the support's volume, since the substituted ratio is 1 wherever
  early exceeds precontrast, and 1 > 0.9;
* `V_real` = the voxels of the support whose SER exceeds the stored threshold.

Endpoints are the published four: T0, T1, absolute change, relative reduction.
The comparison is the published estimator, `theme4_pcr_auc.compare`, with the
published seed 20260918 and 10 000 replicates, and the published sign
convention, including its rule that the relative reduction is **undefined, not
zero**, when the T0 volume is zero.

The connectivity operation of §1 item 4 is **not** reconstructed. The archive's
support is not one connected object — in a dense breast it is scattered
enhancing parenchyma — so keeping only the largest component would build
something the archive never computed. This is stated as a limit, not worked
around.

## 4. Part C: the fully automatic redraw

The same four steps, except that the box comes from the images instead of from
the archive: for each of the three array axes, the maximum-intensity projection
of `early - pre`; on each projection the pixels above that projection's 95th
percentile; the largest 4-connected component; its bounding box; the three
boxes intersected. *The 95th percentile is the percentile the published
background rule already uses [9]; it is not a value chosen here.*

This uses no reader input of any kind, and **it is expected to fail**, because
the largest bright component is the contralateral breast whenever that side
enhances more. It is reported with its Dice distribution and the count of
examinations whose redrawn support is disjoint from the archive's. It is kept
because the honest answer to the referee includes what a fully unsupervised
reconstruction does.

## 5. Cohort

**All 117 participants and all 234 examinations.** Every extracted pair file is
on disk, so there is no subset to select and therefore no selection to defend.
The referee asked for a subset; the whole cohort is strictly better and avoids
a choice.

## 6. What would be unfavourable

Declared now, so that the result cannot be read as favourable by construction.
Each is judged and reported in these words, as pass or fail, before any
narrative is written about the result.

* **F1.** The supplied support is **not** an early-only selection: at any
  examination, fewer than **95 %** of its voxels satisfy `PE > threshold` from
  the precontrast and early images alone. The referee's objection would then be
  sustained on its own terms, and part B's box would be contaminated.
* **F2.** At any margin in the grid, the AUC difference between arms in the
  rebuilt support has a 95 % interval that **excludes zero**. The manuscript's
  conclusion would then depend on the supplied support, which is the referee's
  objection sustained by a different route.
* **F3.** At any margin, the rebuilt support is empty for **more than 10 %** of
  examinations at either visit. The rule would then not work well enough for
  any conclusion to rest on it, favourable or not.
* **F4.** The rebuilt and supplied supports overlap by a median Dice below
  **0.30** at every margin. The rebuilt support would then not be the same
  anatomy, and the comparison would be between two different measurements
  rather than two supports for one.

F4 is not a criterion on part C. Part C is expected to overlap badly; that is
its result.

## 7. The superseded design

The first version of this plan, written and never frozen, made part C primary
and added a reader-seeded variant B in which the component containing the
projection of the archive support's centroid was kept instead of the largest
one. Both are dropped: variant A because of the contralateral-breast failure
of §0, variant B because seeding from the archive's centroid takes more from
the archive than part B's bounding box does while answering less. The declared
margin grid `{0, 2, 4, 8}` and the falsification criteria on emptiness and Dice
are carried over unchanged from that version; the containment criterion F1 is
new, and F2 is the old F1 restated for the rebuilt support.

## 8. Reported quantities

Part A: the distribution of the containment fraction over examinations under
both readings, its minimum, the number of examinations below 0.99, the share of
all support voxels that fail, the number of support voxels with no precontrast
signal, the thresholds actually used and how many examinations each covers, and
the count of examinations whose support size equals the stored FTVpe count.

Part B, for each margin: the number of examinations with a non-empty support;
median and quartiles of the support volume beside the supplied support's;
median Dice and Jaccard against the supplied support; the fraction of the
supplied support the rebuilt one contains; the four AUC differences with
intervals; and the retained fraction.

Part C: the Dice distribution, the number of examinations disjoint from the
supplied support, and the same four AUC differences.

Plus the four falsification criteria.

## 9. Limits, to be stated in the same sentence as the result

* It is not a reader study. No radiologist drew anything, and the rule is an
  algorithm, not a person following the specification's intent. The manuscript
  must call it an algorithmic reconstruction of the documented rule.
* **For non-Philips participants the arrays are already the archive's published
  DCE crop** (§2.2), so parts B and C test the support *within* a crop the
  archive produced, not the crop itself.
* OMIT regions exist only as the archive's mask and cannot be recovered.
* The specification's neighbour-based connectivity is not reconstructed (§3).
* Agreement between the rebuilt and supplied supports does not prove the
  archive's VOIs were drawn without the late image; it shows that a
  late-independent support reaches the same conclusion.

## 10. Amendments

Each amendment records what changed, why, what was known at the time, and the
sha256 of the plan and script it supersedes. The lock it replaces is kept.

### Amendment 1 (2026-11-11): the two examinations with an empty supplied support

**What stopped.** The frozen run halted at participant 49 of 117 on
`assert n, "the supplied support is empty"` in `containment()`. It did not
produce a wrong number; it refused to produce one.

**What was found.** Two examinations of the published cohort have an **empty
supplied support**: `ACRIN-6698-417474` T1 and `ACRIN-6698-490805` T1 both
carry a stored `FTVpe(PEearly>70%)` of 0 voxels and 0 cc, so the archive's
analysis mask excludes every voxel. Both are follow-up visits at which the
archive found no enhancing tissue. They are the same two examinations already
reported as the cohort's zero-reference examinations in the Referee 1 point 4
analysis (`docs/theme4_revision_results.md` §R1.4), and they were identified
here by a stored volume of zero, **with no reference to the outcome label, to
any AUC, or to any quantity this experiment estimates.**

**The rule, now fixed.**

1. *Part A.* The containment fraction is 0/0 and is reported as **undefined**.
   Those examinations are excluded from the distribution, which is therefore
   over **232** of the 234 examinations, and they are listed by name in the
   output. Reporting 1.0 would add two free successes; reporting 0.0 would add
   two free failures. Neither is a measurement.
2. *Parts B and C.* No box exists, so no support can be rebuilt. The
   examination is carried with a volume of **zero in both arms**, which is
   exactly the value the published analysis gives it
   (`real_T1 = sub_T1 = 0` in `theme4_pcr_v1/paired_scores.csv`). The
   participant therefore stays in the paired comparison and the comparison
   remains on all 117 participants. The alternative — dropping the
   participant — would change the cohort of the comparison relative to the
   published one, and the choice between them is settled by fidelity to the
   published analysis, not by either answer.
3. Dice, Jaccard and the retained fraction are undefined against an empty
   support and are excluded from those distributions rather than set to zero.
4. The count of examinations carried this way is reported at every margin, so
   that the emptiness criterion F3 cannot be satisfied by hiding them among the
   rebuilt supports: F3 counts **rebuilt** supports that came out empty, and
   these two never had a box to rebuild.

**Superseded hashes.** Plan
`d3308df81438cf6e9158af5107dffacedf3542e8045e5709d84c30100a2dad24`, script
`0d103a128c15c48888726a8b55203fece7c60cd76ff5a554a392d0ddc5f915ff`. The lock
they were frozen in is kept as `docs/theme4_early_only_lock_v1.superseded-1.json`.

**What had been seen when this was written.** The containment fractions of the
first 48 participants, printed by the halted run: all between 0.992 and 1.000.
No part B comparison, no Dice distribution and no part C result existed. The
falsification criteria of §6 are unchanged.

---

[9] The ACRIN 6698 / I-SPY 2 analysis specification cited in the manuscript's
reference list as the source of the PE threshold, the background rule, the
connectivity operation and the VOI description.
