# FASC: Flexible Adaptive Stable Clustering Algorithm

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Patent](https://img.shields.io/badge/Patent-Granted-brightgreen)](https://cpquery.cnipa.gov.cn/)

**FASC** (Flexible Adaptive Stable Clustering) is a novel, high-performance clustering algorithm designed specifically for high-dimensional mass spectrometric big data.

Developed at **Southern University of Science and Technology (SUSTech)**, FASC addresses the challenges of clustering single-particle mass spectrometry (SPMS) data by combining density-based initialization with an adaptive resonance theory-inspired learning process. It effectively handles million-scale datasets with high stability and speed.

## Key Features

* **Adaptive Clustering:** Automatically determines the number of clusters based on data density and similarity thresholds, eliminating the need for pre-defined cluster counts.
* **Dual-Cosine Metric:** Incorporates a specialized "Dual-Cosine" similarity algorithm (Positive/Negative ion modes) optimized for chemical fingerprinting of atmospheric particles.
* **High Performance:** Fully vectorized and parallelized implementation using MATLAB's `parfor` and `parfeval` to accelerate processing of large datasets.
* **User-Friendly GUI:** Includes a complete App Designer interface for easy data loading, parameter tuning, and real-time monitoring.
* **Advanced Visualization:** Built-in tools for generating publication-quality similarity heatmaps, cluster distribution histograms, and convergence plots.

## Repository Structure

To facilitate reproducibility, the project is organized as follows:

```text
.
├── src/                # Core algorithm source code (FASC.m, FASC_GUI.m)
├── utils/              # Helper functions for visualization and analysis
├── data/               # Demo datasets (for reproducing examples)
├── LICENSE             # GNU GPLv3 License text
└── README.md           # Project documentation
````

## System Requirements

  * **MATLAB:** R2021b or later (Recommended for best App Designer compatibility).
  * **Required Toolboxes:**
      * *Statistics and Machine Learning Toolbox* (for `pdist2`, distance metrics).
      * *Parallel Computing Toolbox* (for parallel acceleration).

## Installation

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/YourUsername/FASC.git](https://github.com/YourUsername/FASC.git)
    cd FASC
    ```

2.  **Setup MATLAB Path:**
    Open MATLAB, navigate to the `FASC` folder, and run the following command to add all subfolders to your path:

    ```matlab
    addpath(genpath(pwd));
    savepath;
    ```

## Usage

### 1\. Graphical User Interface (Recommended)

The GUI is the easiest way to explore your data and tune parameters.

```matlab
% In MATLAB Command Window:
FASC_GUI
```

  * **Load Data:** Supports `.mat` files (containing a numeric matrix) and `.csv` files.
  * **Parameters:** Adjust `Sim (intra)` and `Sim (inter)` directly in the sidebar.
  * **Run:** Click **RUN ANALYSIS**. The interface runs asynchronously, allowing you to monitor progress logs.

-----

### 2\. Command Line Interface (Scripting)

For batch processing or integration into existing pipelines, use the core function directly.

#### Example Script

```matlab
% 1. Load the demo dataset (included in the 'data' folder)
% Ensure you have added the folders to your path as described in Installation
load(fullfile('data', 'example_data.mat'), 'dataMatrix'); 

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

% 4. Visualize Results (Optional)
figure; 
histogram(labels);
title('Cluster Distribution');
```

#### Output Arguments

  * `centers`: The resulting cluster centroid matrix.
  * `counts`: The number of particles in each cluster.
  * `labels`: Vector containing the cluster assignment ID for each particle (0 = outlier).
  * `info`: Structure containing iteration history and convergence metrics.

## Data Availability

  * **Demo Data:** A subset of particle data is provided in the `data/` folder of this repository to verify the algorithm's functionality.
  * **Full Dataset:** The complete high-dimensional mass spectrometric dataset used in the manuscript is available at **[Zenodo/Figshare Repository Name]** under DOI: **[10.xxxx/xxxx]**.

## License & Patent Notice

### Copyright © 2025 Shao Shi, Southern University of Science and Technology (SUSTech).

#### 1\. Open Source License (GPLv3)

This software is provided under the **GNU General Public License v3.0 (GPLv3)**.

  * **You are free to:** Run, study, share, and modify the software for academic and non-commercial research.
  * **Conditions:** Any modifications or derived works must be open-sourced under the same GPLv3 license.

#### 2\. Patent Notice

> **The core algorithm implemented in this software is protected by China Invention Patent.**
>
>   * **Patent No.:** ZL 2025 1 0839139.3
>   * **Status:** Granted
>   * **Academic Use:** Free and unrestricted for non-commercial research and educational purposes.
>   * **Commercial Use:** Any commercial use (including integration into commercial software/hardware) is **strictly prohibited** without a separate commercial license. Please contact the author for licensing inquiries.

## Citation

If you use FASC or its GUI in your research, please cite our paper:

> **Flexible Adaptive Stable Clustering (FASC) for High-Dimensional Mass Spectrometry Data**
> Shao Shi, et al.
> *Manuscript Submitted / Under Review*, 2025.
> (The full citation and DOI will be updated upon publication.)

## Contact

**Shao Shi (石邵)**
School of Environmental Science and Engineering
Southern University of Science and Technology (SUSTech)
Shenzhen, China
Email: [Your.Email@sustech.edu.cn]
GitHub: [https://github.com/YourUsername](https://github.com/YourUsername)
