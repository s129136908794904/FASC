
# FASC: Flexible Adaptive Stable Clustering Algorithm

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Patent](https://img.shields.io/badge/Patent-Granted-brightgreen)](https://cpquery.cnipa.gov.cn/)

**FASC** (Flexible Adaptive Stable Clustering) is a novel, high-performance clustering algorithm designed specifically for high-dimensional mass spectrometric big data.

Developed at Southern University of Science and Technology (SUSTech), FASC addresses the challenges of clustering single-particle mass spectrometry (SPMS) data by combining density-based initialization with an adaptive resonance theory-inspired learning process. It effectively handles million-scale datasets with high stability and speed.

## Key Features

* **Adaptive Clustering:** Automatically determines the number of clusters based on data density and similarity thresholds, eliminating the need for pre-defined cluster counts.
* **Dual-Cosine Metric:** Incorporates a specialized "Dual-Cosine" similarity algorithm (Positive/Negative ion modes) optimized for chemical fingerprinting of atmospheric particles.
* **High Performance:** Fully vectorized and parallelized implementation using MATLAB's `parfor` and `parfeval` to accelerate processing of large datasets.
* **User-Friendly GUI:** Includes a complete App Designer interface (`FASC_GUI`) for easy data loading, parameter tuning, and real-time monitoring.
* **Advanced Visualization:** Built-in tools for generating publication-quality similarity heatmaps, cluster distribution histograms, and convergence plots.

## System Requirements

* **MATLAB:** R2021b or later (Recommended for best App Designer compatibility).
* **Required Toolboxes:**
    * *Statistics and Machine Learning Toolbox* (for `pdist2`, distance metrics).
    * *Parallel Computing Toolbox* (for parallel acceleration).

## Installation

1.  Clone this repository to your local machine:
    ```bash
    git clone [https://github.com/YourUsername/FASC.git](https://github.com/YourUsername/FASC.git)
    ```
2.  Add the `FASC` folder (and subfolders) to your MATLAB path.

## Usage

### 1. Graphical User Interface (Recommended)
The GUI is the easiest way to explore your data and tune parameters.

```matlab
% In MATLAB Command Window:
FASC_GUI
````

  * **Load Data:** Supports `.mat` files (containing a numeric matrix) and `.csv` files.
  * **Parameters:** Adjust `Sim (intra)` and `Sim (inter)` directly in the sidebar.
  * **Algorithm:** Select `dual-cosine` for standard SPMS data or other metrics (e.g., `euclidean`, `cosine`) for general data.
  * **Run:** Click **RUN ANALYSIS**. The interface runs asynchronously, allowing you to monitor progress logs.

-----

### 2\. Command Line Interface (Scripting)

For batch processing or integration into existing pipelines, use the core function directly.

#### Input Data Format

  * **Rows:** data vector, such as Individual particles (spectra).
  * **Columns:** features, such as Mass-to-charge ratio (m/z) .

#### Example Script


```matlab
% 1. Load your data (n_particles x n_features)
load('datamatrix.mat', 'dataMatrix'); 

% 2. Define Indices for Dual-Cosine only (if applicable)
% Example: Features 1-300 are Positive spectrum, 301-600 are Negative
idx_pos = 1:300; 
idx_neg = 301:600;

% 3. Run FASC
[centers, counts, labels, info] = FASC(...
    dataMatrix, ...       % Input data (Single or Double precision)
    idx_pos, idx_neg, ... % Column indices for dual-ion mode
    0.7, ...              % Inter-cluster similarity threshold (Merger)
    0.7, ...              % Intra-cluster similarity threshold (Vigilance)
    8, ...                % Initial cluster seed limit
    50, ...               % Maximum cluster number
    200, ...              % Maximum iterations
    'DASS', ...           % Strategy: 'DASS' (Density-Adaptive) or 'SF' (Sim-First)
    2, ...                % Minimum cluster volume to be valid
    'dual-cosine');       % Similarity algorithm
```


#### Output Arguments

  * `centers`: The resulting cluster centroid matrix.
  * `counts`: The number of particles in each cluster.
  * `labels`: Vector containing the cluster assignment ID for each particle (0 = outlier).
  * `info`: Structure containing iteration history and convergence metrics.

## License & Patent Notice

### Copyright © 2025 Shao Shi, Southern University of Science and Technology (SUSTech).

#### 1\. Open Source License (GPLv3)

This software is provided under the **GNU General Public License v3.0 (GPLv3)**.

  * **You are free to:** Run, study, share, and modify the software.
  * **Conditions:** If you distribute modified versions or software based on FASC, you must open-source them under the same GPLv3 license.
  * See the `LICENSE` file for the full text.

#### 2\. Patent Notice

> **The core algorithm implemented in this software is protected by China Invention Patent.**
>
>   * **Patent No.:** ZL 2025 1 0839139.3
>   * **Status:** Granted
>   * **Commercial Use:** Any commercial use of this algorithm (including but not limited to integration into commercial mass spectrometry software, hardware, or paid analytical services) is **strictly prohibited** without a separate commercial license. Please contact the author for licensing inquiries.

## Citation

If you use FASC or its GUI in your research, please cite our paper:

> **[Paper Title Placeholder]**
> Shao Shi, [Other Authors].
> *[Journal Name]*, 202X.
> DOI: [DOI Placeholder]

## Contact

**Shao Shi (石邵)**
School of Environmental Science and Engineering
Southern University of Science and Technology (SUSTech)
Shenzhen, China
Email: [Your Email Address]
GitHub: [Your GitHub Profile URL]

