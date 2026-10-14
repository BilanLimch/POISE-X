# POISE-X

Plug-and-play Optimizer for Imbalance-aware Subtyping of Endometrial cancer.
A leakage-controlled benchmarking pipeline for multi-omics feature selection
and reconstruction of the TCGA molecular subtypes of endometrial carcinoma
(TCGA-UCEC), in which the metaheuristic search engine (X) is an
interchangeable component. Every data-dependent step (prefiltering,
standardization, wrapper search, classifier choice) is fitted inside the outer
training fold; the choice among optimizers is made on the outer test folds and
is evaluated separately by leave-one-fold-out selection (`ucec_app.selection_rule`). Nineteen optimizers from the MEALPY library are
provided as instantiations (`POISE-GWO`, `POISE-DE`, ...).

## Requirements

Python 3.11 and the packages in `requirements.txt`. The reported results were
produced with NumPy 1.26.0, SciPy 1.16.3, pandas 2.x, scikit-learn 1.9.0,
XGBoost 3.2.0 and MEALPY 3.0.3.

    python -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt

## Data

Download the following files from the GDC TCGA-UCEC hub of UCSC Xena
(https://xenabrowser.net/datapages/) into one directory, referred to below as
`UCEC_RAW_DIR`:

| File | Content |
|---|---|
| `TCGA-UCEC.star_fpkm.tsv.gz` | gene expression, log2(FPKM+1) |
| `TCGA-UCEC.mirna.tsv.gz` | miRNA expression, log2(RPM+1) |
| `TCGA-UCEC.methylation450.tsv.gz` | Illumina 450K beta values |
| `TCGA-UCEC.gene-level_ascat3.tsv.gz` | gene-level copy number (ASCAT3) |
| `TCGA-UCEC.somaticmutation_wxs.tsv.gz` | whole-exome somatic mutations |
| `TCGA-UCEC.clinical.tsv.gz` | clinical table |

Molecular-subtype labels come from cBioPortal (study
`ucec_tcga_pan_can_atlas_2018`, column *Subtype*). A copy is provided in
`data/TCGA-UCEC.molecular_subtype.tsv`. The strict-intersection cohort table
used by the pipeline is provided in `data/final_cohort.tsv` (369 patients); it
can be rebuilt from the raw files with

    cp data/TCGA-UCEC.molecular_subtype.tsv $UCEC_RAW_DIR/
    UCEC_RAW_DIR=$UCEC_RAW_DIR python scripts/build_cohort.py

Either way, `final_cohort.tsv` must be present in `UCEC_RAW_DIR`.

The revision analyses additionally use the MLH1 probes of the 450K array;
their list (from the public GDC hg38 probe map) is provided in
`data/mlh1_probes.tsv`, and `scripts/extract_mlh1_probes.py` extracts their
beta values from the raw methylation file.

## Quick check without data

    python -m ucec_app.baseline_smoke_test   # every optimizer respects the evaluation budget
    python -m ucec_app.smoke_test            # end-to-end on synthetic data (a few minutes)

## Running the pipeline

All commands are run from this directory. Paths are controlled by three
environment variables: `UCEC_RAW_DIR` (raw files and `final_cohort.tsv`),
`UCEC_CACHE_DIR` (default `./data`) and `UCEC_RESULTS_DIR` (default
`./results`).

1. Build the label-free feature cache once (reads the raw files, keeps one
   primary-tumor sample per patient, drops features with missing or constant
   values; ~30 min, mostly for the methylation file):

       export UCEC_RAW_DIR=/path/to/raw
       python -m ucec_app.data --build

2. Run one unit of the nested cross-validation: one optimizer, one outer fold,
   one seed. Run index `i` in 0..19 uses outer fold `i // 4` and seed
   `2026 + 98*i`. Each unit takes between 0.1 and 2.5 h depending on the
   optimizer and writes `results/<ALG>/fold<k>_seed<s>.json`.

       python -m ucec_app.nested_cv --algorithm GWO --run 0

   `run_all.sh` loops over all optimizers and runs (380 units) and skips units
   whose result file exists, so it can be restarted:

       JOBS=8 OMP_NUM_THREADS=2 ./run_all.sh

3. The no-selection reference (each classifier on the full candidate pool of
   every outer fold):

       python -m ucec_app.full_pool_baseline

4. Aggregate: per-run table, per-optimizer summary, Friedman test and
   Holm-corrected paired Wilcoxon tests against the best optimizer, feature
   selection frequencies, classifier frequencies, Jaccard stability and
   convergence curves, all written to `results/summary/`:

       python -m ucec_app.aggregate

5. Consensus-panel experiment (fixed panels of the most frequently selected
   features of each fold's own searches, compared with filter and random
   panels of the same size):

       python -m ucec_app.consensus_validation --alg GWO --alg ALL

## Revision analyses

`run_revision.sh` runs everything below after the main experiment (steps 1-4
above) has finished; the individual steps are:

1. Sensitivity searches of POISE-GWO and POISE-ACO (160 units): parsimony
   weight alpha = 0.8, 0.95, 0.99 and an omics-only candidate pool without the
   five clinical variables (Section 4.8, Tables S4 and S5). Same folds and
   seeds as the main experiment; with `--alpha 0.9` a unit reproduces the
   corresponding unit of the main experiment exactly.

       python -m ucec_app.search_variant --algorithm GWO --run 0 --alpha 0.95
       python -m ucec_app.search_variant --algorithm ACO --run 0 --omics-only

2. Reference analyses without search: single omics blocks, leave-one-block-out
   and two-block pools, the surrogate of the clinical molecular classification
   (POLE exonuclease-domain hotspot, mismatch-repair deficiency, TP53
   mutation), omics-only and clinical-only LR, and the noise of the inner-fold
   macro-F1 estimate (Section 4.8, Tables S2-S4, Section 3.5):

       python scripts/extract_mlh1_probes.py
       python -m ucec_app.revision_baselines

3. Leave-one-fold-out evaluation of the configuration-selection rule
   (Section 4.7, Table S1):

       python -m ucec_app.selection_rule

4. Chance reference of the run-to-run Jaccard index (Section 4.5):

       python -m ucec_app.jaccard_chance

5. Summary of the sensitivity searches and Figure S1:

       python -m ucec_app.revision_aggregate

Results are written to `results/revision/` and `results/summary/`. Logistic
regression fits can differ in the third decimal with the number of
linear-algebra threads; the reference models were fitted with the library
defaults.

## Layout

    ucec_app/data.py                  raw-file loading, primary-sample selection, feature cache
    ucec_app/prefilter.py             in-fold candidate pool (ANOVA F top-500 per block, mutation >= 3%)
    ucec_app/feature_selection.py     encoding, decoding, classifier pool, macro-F1 fitness (inner 3-fold CV)
    ucec_app/baselines.py             the 19 MEALPY optimizers under one evaluation budget
    ucec_app/nested_cv.py             one (optimizer, fold, seed) unit of the outer 5-fold protocol
    ucec_app/full_pool_baseline.py    no-selection reference
    ucec_app/aggregate.py             summary tables and statistical tests
    ucec_app/consensus_validation.py  consensus-panel experiment
    ucec_app/search_variant.py        revision: searches with another alpha or without clinical variables
    ucec_app/revision_baselines.py    revision: data-layer, surrogate, omics-only and inner-noise analyses
    ucec_app/selection_rule.py        revision: leave-one-fold-out evaluation of the selection rule
    ucec_app/jaccard_chance.py        chance reference of the panel Jaccard index
    ucec_app/revision_aggregate.py    revision: summary of the sensitivity searches, Figure S1
    ucec_app/*_smoke_test.py          quick self-checks
    scripts/build_cohort.py           strict-intersection cohort table
    scripts/extract_mlh1_probes.py    MLH1 probe beta values for the surrogate
    data/                             cohort table, subtype labels, MLH1 probe list
    run_all.sh                        runs all 380 units
    run_revision.sh                   runs all revision analyses

## Adding an optimizer

Register a MEALPY class in `ALGORITHMS` in `ucec_app/baselines.py` (module
path and class name); any optimizer that works on a box-constrained
continuous problem can be plugged in without changing the rest of the
pipeline.
