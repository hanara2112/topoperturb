# Data: what it is, and what it is not

The analysis uses the **prepared** CausalBench tables published on Hugging Face
(`huggingbahl21/causalbench`, one gzipped TSV per cell line), pinned to a single
revision and verified by SHA-256 in [`topoperturb/data.py`](../topoperturb/data.py).

## What one table contains

- One **cell per row**; one **gene per column**, plus a final `intervention` column.
- `intervention` is the knocked-down gene (an Ensembl ID) or `non-targeting` for a
  control. A small number of rows are labelled `excluded` and are dropped.
- Values are already **library-normalised and `log1p`-transformed** counts
  (finite, non-negative). We add no normalisation, log, imputation or filtering.

`data.audit()` scans the whole table once and refuses it unless every value is
finite and non-negative, every cell has a label, and controls are present.

## Upstream selection (important for interpretation)

The published tables are **already efficacy-filtered**: CausalBench selected strong
perturbations and kept perturbed cells partly by their target-gene response. A
detected difference therefore describes *selected responders*, not every cell
assigned a guide. The control-vs-control calibration (step 4) does **not**
reproduce this selection mechanism, so it bounds only the ordinary false-positive
rate, not any bias introduced by responder selection.

## What the table does **not** contain

- **No guide ID, no batch, no cell barcode.** The only grouping is `intervention`.
- `source_rows` saved with the clouds is the zero-based row index in the TSV, kept
  for traceability — it is *not* a biological identifier.

Consequence: we cannot stratify permutations by guide or batch, and cannot do
per-guide validation, from these inputs alone. Doing so would need richer upstream
data. This is a limitation we state, not one we have worked around.

## From cells to clouds (step 1)

`experiments/01_prepare_data.py` keeps the controls plus the perturbations we study
(the eight K562 pilot knockdowns; for RPE1 the four shared knockdowns plus one
strong positive control) and projects every cell onto **30 principal components**.

PCA is fit **label-blind** with `fit="balanced"`: an equal number of cells per
group define the axes, so one large group cannot dominate the representation. We
fix this one representation deliberately rather than searching over PCA/graph
settings for favourable p-values; `data.embed` exposes `n_pcs`, `fit` and
`drop_genes` for anyone who wants to probe sensitivity. The control covariance
spectrum saved alongside the clouds is what the synthetic benchmark (step 2a) uses
to make its base Gaussian realistic.
