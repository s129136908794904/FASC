# Versioned Benchmark Results

## Files

`external_label_benchmark.csv` contains 252 rows:

- 27 FASC-SF rows: 9 datasets x 3 merge quantiles
- 27 FASC-DASS rows: 9 datasets x 3 merge quantiles
- 27 HDBSCAN rows: 9 datasets x 3 minimum-cluster fractions
- 27 OPTICS rows: 9 datasets x 3 minimum-cluster fractions
- 27 BIRCH rows: 9 datasets x 3 radius quantiles
- 27 batch DP-means rows: 9 datasets x 3 penalty quantiles
- 90 oracle-K K-means rows: 9 datasets x 10 seeds

`fasc_capacity_sensitivity.csv` contains 54 rows:

- 2 FASC assignment variants
- 3 capacity multipliers: 0.75, 1.0, and 1.5
- 9 datasets

The complete parameter object is stored in the `parameters` column of every
row.

## Fixed Primary Settings

Macro means across the nine datasets are:

| Method | Coverage | ARI | AMI | Pairwise F1 | Mean Kout |
|---|---:|---:|---:|---:|---:|
| FASC-SF | 0.990 | 0.333 | 0.587 | 0.370 | 39.2 |
| FASC-DASS | 0.990 | 0.334 | 0.589 | 0.372 | 38.7 |
| HDBSCAN | 0.600 | 0.289 | 0.311 | 0.370 | 7.1 |
| OPTICS-xi | 0.449 | 0.090 | 0.084 | 0.159 | 4.4 |
| BIRCH, no global K | 1.000 | 0.262 | 0.537 | 0.293 | 119.6 |
| Batch DP-means | 1.000 | 0.206 | 0.266 | 0.319 | 7.6 |
| Oracle-K K-means | 1.000 | 0.447 | 0.570 | 0.514 | 12.6 |

Oracle-K K-means receives the reference class count and is not an adaptive
baseline. Labels are otherwise used only to calculate the recorded metrics.

The neighboring parameter settings and public-data capacity analysis are
reported in the Supplementary Information. The protocol and interpretation
of rejected observations are documented in `../README.md`.
