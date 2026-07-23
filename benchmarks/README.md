# Public External-Label Benchmark

This benchmark evaluates the cosine score/representative configuration of
FASC on public, non-SPMS datasets. External class labels are withheld from
preprocessing, parameter construction, seeding, capacity selection, and
clustering. They are used only to calculate evaluation metrics.

The experiment measures agreement with the supplied reference
classifications. It is not a distribution-free proof that arbitrary labels
can be recovered from an unlabeled feature matrix.

## Datasets

The reported suite contains:

| Dataset | Source | Reported sample size |
|---|---|---:|
| Digits | scikit-learn | 1,797 |
| Optdigits | OpenML 28 | 2,000 |
| Pendigits | OpenML 32 | 2,000 |
| Segment | OpenML 36 | 2,000 |
| Vehicle | OpenML 54 | 846 |
| ISOLET | OpenML 300 | 2,000 |
| MFeat-Factors | OpenML 12 | 2,000 |
| Letter | OpenML 6 | 2,000 |
| MNIST test | official test split | 2,000 |

Larger datasets are deterministically subsampled without replacement before
labels are inspected. The MNIST loader uses rows 60,001 through 70,000 of a
local 70,000-row export, corresponding to the official test split, and then
applies the same maximum sample count.

The `pilot` suite contains Iris, Wine, Breast Cancer, and Digits. It is for
local protocol checks and is not included in the manuscript results.

## Common Representation

Nonnegative image and trajectory arrays use feature-wise max-absolute
scaling. Other continuous arrays are standardized feature-wise. Rows are
then L2 normalized for the reference methods; FASC applies the equivalent
normalization inside its cosine score.

For two unit vectors, squared Euclidean distance is twice one minus cosine
similarity. This gives the methods a common input representation without
making their objectives, representative updates, or cluster-selection rules
identical.

## FASC Settings

- Base capacity:
  `Kmax = min(100, max(12, ceil(sqrt(N))))`
- Initial seed budget: `min(8, Kmax)`
- Minimum returned support: `max(2, ceil(0.001 * N))`
- Assignment threshold: quantile 0.01 of each observation's unlabelled
  10th-nearest-neighbour cosine similarity
- Merge-threshold quantiles: 0.50, 0.75, and 0.90
- Primary merge quantile: 0.90
- DASS relative-support setting: `alpha = 0.05`,
  `lambda = alpha / N`, and `phi(n) = n`
- Iteration safeguard: 500 transitions

SF is the basic FASC assignment rule. DASS is an optional support preference.
The primary capacity multiplier is 1.0. Multipliers 0.75 and 1.5 are reported
in a separate public-data sensitivity file.

No FASC configuration is selected by maximizing an external-label metric.

## Reference Methods

All adaptive reference settings are label-blind.

- HDBSCAN: minimum-cluster fractions 0.005, 0.01, and 0.02; primary 0.01.
  `min_samples = max(5, ceil(0.005 * N))` is held independent of the
  minimum-cluster-size grid.
- OPTICS: the same minimum-cluster fractions and fixed `xi = 0.05`, with the
  same independent `min_samples` rule; primary fraction 0.01.
- BIRCH: `n_clusters=None`; radius thresholds are unlabelled
  10th-neighbour distance quantiles 0.25, 0.50, and 0.75; primary 0.50.
- Deterministic batch DP-means: the squared-distance penalty is an
  unlabelled pairwise-distance quantile 0.10, 0.25, or 0.50; primary 0.25.
  The implementation uses canonical initialization, farthest-point
  promotion, and tie handling, and never receives the reference class count.
- Oracle-K K-means: receives the reference class count and is repeated over
  ten fixed seeds. It is a separate feature-separability reference, not a
  like-for-like adaptive baseline.

The sensitivity ranges are descriptive. They show how each method moves
between coverage, fragmentation, and agreement operating points; they are
not used to choose a best score after labels are observed.

## Metrics

Primary reported quantities are:

- coverage and outlier fraction;
- all-sample ARI and AMI, with each rejected observation represented as a
  distinct singleton;
- pairwise precision, recall, and F1, where a positive predicted pair
  requires two assigned observations in the same nonzero cluster; and
- returned cluster count and minimum returned support.

The CSV also records assigned-only ARI/AMI, homogeneity, completeness,
V-measure, macro and minimum class coverage, runtime, termination state,
cycle length, and the full parameter object.

Majority-vote mapped accuracy is excluded because a fragmented partition can
map many small clusters to one class and obtain a favorable score.

## Reproduce the Fixed Grid

Install NumPy, SciPy, and scikit-learn 1.3 or later. From the repository root:

```bash
python benchmarks/run_external_benchmark.py \
  --suite public \
  --mnist-dir "PATH/TO/MNIST/rawData" \
  --overwrite
```

The script defaults match the reported 2,000-observation fixed grid. Without
`--mnist-dir`, the eight datasets available from scikit-learn/OpenML are run.

Run the FASC capacity sensitivity separately:

```bash
python benchmarks/run_external_benchmark.py \
  --suite public \
  --mnist-dir "PATH/TO/MNIST/rawData" \
  --fasc-only \
  --threshold-quantiles 0.01 \
  --merge-threshold-quantiles 0.90 \
  --dass-alphas 0.05 \
  --fasc-capacity-multipliers 0.75,1.0,1.5 \
  --output benchmarks/results/fasc_capacity_sensitivity.csv \
  --overwrite
```

Run the local benchmark tests with:

```bash
python -m unittest discover -s benchmarks -p "test_*.py" -v
```

The script writes one CSV row per dataset, method, and parameter
configuration. Complete versioned result files are in `results/`.
