# FASC: Flexible Adaptive Stable Clustering

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17844844.svg)](https://doi.org/10.5281/zenodo.17844844)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Patent](https://img.shields.io/badge/Patent-Granted-brightgreen)](https://cpquery.cnipa.gov.cn/)

**FASC** (Flexible Adaptive Stable Clustering) is a novel, high-performance clustering algorithm designed specifically for high-dimensional mass spectrometric big data.

Developed at **Southern University of Science and Technology (SUSTech)**, FASC addresses the challenges of clustering single-particle mass spectrometry (SPMS) data by combining density-based initialization with an adaptive resonance theory-inspired learning process. It effectively handles million-scale datasets with high stability and speed.

## Key Features

* **Dual-Language Support:** Fully implemented in both **Python** (for seamless open-science integration) and **MATLAB** (for HPC-optimized matrix streaming).
* **Adaptive Clustering:** Automatically determines the number of clusters based on data density and similarity thresholds, acting as a high-pass density filter.
* **Metric Flexibility:** Supports arbitrary bounded symmetric kernels, including a specialized "Dual-Cosine" similarity optimized for dual-polarity atmospheric MS.
* **Deterministic Stability:** Guarantees mathematically reproducible convergence independent of data presentation order.

## Repository Structure

To facilitate reproducibility, the project is organized as follows:

```text
.
├── python/             # Python implementation (Core algorithm & Jupyter tutorials)
├── matlab/             # MATLAB implementation (Core algorithm, Utils, & GUI)
├── data/               # Demo datasets for reproducing examples
├── LICENSE             # GNU GPLv3 License text
└── README.md           # Project documentation
```

## System Requirements

  * **Hardware:**
      * **CPU:** Multi-core processor recommended (4+ cores) for parallel processing.
      * **RAM:** Depending on the size of the dataset, e.g., we use ~200 GB RAM for the clustering of 25 million 300-D vectors (fp64) on CentOS Linux 7.
      * **Non-standard Hardware:** None required.
  * **Operating System:** Windows 10 or above, macOS (10.15+), or Linux (Ubuntu 20.04+, CentOS Linux release 7.5.1804 (Core)+).

### Python Environment
* **Python:** 3.8 or higher.
* **Dependencies:** `numpy`, `scipy`, `matplotlib`, `jupyter`
* Install dependencies via: `pip install -r python/requirements.txt`

### MATLAB Environment
* **MATLAB:** R2021b or later (Windows), 2020b or later (Linux).
* **Toolboxes:** *Statistics and Machine Learning Toolbox*, *Parallel Computing Toolbox*.

## Data Availability

* **Demo Data:** A subset of the high-dimensional mass spectrometric dataset (dataMatrix.csv) is provided in the `data/` folder to verify the algorithm's functionality.
* **Full Dataset:** The complete 25-million spectra dataset used in the manuscript is available at **[Zenodo]** under DOI: **10.5281/zenodo.17788367**.

---

## Usage: Python (Recommended for Data Science Pipelines)

The Python implementation is optimized for sparse matrix operations and integrates seamlessly into modern data science workflows.

### Quick Start (Jupyter Notebook)
The easiest way to get started is to run the interactive tutorial:
```bash
cd python
jupyter notebook tutorial.ipynb
```

### Python Scripting Example
```python
import numpy as np
from fasc_core import run_fasc

# 1. Load the demo data
data_path = './dataMatrix.csv'

print(f"Loading dataset from {data_path}...")
data_matrix = np.loadtxt(data_path, delimiter=',')
print(f"Data loaded successfully! Shape: {data_matrix.shape} (Spectra x Features)")


# 2. Define Parameters

sim_inter = 0.70      # Merge threshold (Inter-cluster)
sim_intra = 0.70      # Assignment threshold (Intra-cluster)
max_clust = 50        # Capacity budget
strategy  = 'DASS'    # Density-Augmented Similarity Selection
# Here we use dual-cosine metric as an example, for other similarity algorithms, idx_pos and idx_neg is not required.
algo      = 'dual-cosine' 
idx_pos = np.arange(0, 300)   # Positive spectrum features
idx_neg = np.arange(300, 600) # Negative spectrum features

# 3. Run FASC
centers, counts, labels, iter_data, total_iters = run_fasc(
    data_matrix=data_matrix, 
    idx_pos=idx_pos, idx_neg=idx_neg, 
    sim_inter=sim_inter, sim_inner=sim_intra, 
    init_limit=10, max_clusters=max_clust, max_iter=200, 
    strategy=strategy, min_vol=10, algo=algo
)

print(f"FASC identified {len(counts)} pure clusters.")
```

---

## Usage: MATLAB (Recommended for HPC & GUI Users)

Assuming MATLAB is already installed.

1.  **Clone the repository:**
In terminal, Powershell or bash, run the following to get the codes:
    ```bash
    git clone https://github.com/s129136908794904/FASC.git
    cd FASC
    ```

2.  **Setup MATLAB Path:**
    Open MATLAB, navigate to the `FASC` folder, and run the following command to add all subfolders to your path:

    ```matlab
    addpath(genpath(pwd));
    savepath;
    ```

![FASC GUI Screenshot](assets/FASC_GUI_screenshot.png)

*Figure 1: The FASC GUI showing the main control panel (left) and the real-time visualization of clustering results (right) on the demo dataset.*

## Usage

### 1\. Graphical User Interface
FASC includes a user-friendly App Designer interface that allows researchers to visualize high-dimensional data clustering in real-time without writing code.

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
% In MATLAB:
% 1. Setup Environment
addpath(genpath(pwd)); % Ensure src/ and utils/ are in path
load(fullfile('data', 'dataMatrix.mat'), 'dataMatrix'); 

% 2. Define Parameters
% Define Column Indices (Required for 'dual-cosine', ignored for others)
idx_pos = 1:300;    % Positive spectrum features
idx_neg = 301:600;  % Negative spectrum features

% Thresholds & Constraints
sim_inter  = 0.70;      % Merge threshold (Inter-cluster)
sim_intra  = 0.70;      % Vigilance threshold (Intra-cluster)
seed_limit = 8;        % Max initial seeds
max_clust  = 50;        % Max allowed clusters
max_iter   = 200;       % Max iterations
min_vol    = 2;         % Minimum particles to form a valid cluster

% Strategy: 'DASS' (Density-Adaptive) or 'SF' (Similarity-First)
strategy   = 'DASS';    

% Metric: 'dual-cosine', 'cosine', 'euclidean', 'l1-norm', etc.
algorithm  = 'dual-cosine'; 

% 3. Run FASC Analysis
fprintf('Starting FASC Analysis...\n');
tic;
[centers, counts, labels, info] = FASC(...
    dataMatrix, ...           % dataMatrix
    idx_pos, idx_neg, ... % Ion mode indices
    sim_inter, ...      % Similarity Threshold (Inter)
    sim_intra, ...      % Similarity Threshold (Intra)
    seed_limit, ...     % Initial Seed Limit
    max_clust, ...      % Max Cluster Limit
    max_iter, ...       % Max Iteration Limit
    strategy, ...       % Optimization Strategy
    min_vol, ...        % Min Cluster Volume
    algorithm);         % Similarity Algorithm
toc;

% 4. Generate Visualizations (Using built-in helper functions)
% These functions generate the same plots as the GUI.

% Define output directory for figures
outDir = fullfile(pwd, 'results', filesep);
if ~exist(outDir, 'dir'), mkdir(outDir); end

% A. Plot Convergence & Outlier Stats
% Usage: Clustering_iterInfoAnalyzer(infoStruct, ifPlot, ifSave, FileName, SaveDir)
Clustering_iterInfoAnalyzer(info, true, true, 'Demo_Convergence', outDir);

% B. Plot Similarity Heatmap & Distribution
% Usage: Clustering_clusterListCountsHister(counts, centers, sim_threshold, algo, pos, neg, ifPlot, ifSave, FileName, SaveDir)
Clustering_clusterListCountsHister(...
    counts, centers, sim_inter, algorithm, idx_pos, idx_neg, ...
    true, true, 'Demo_Heatmap', outDir);

fprintf('Analysis complete. Results saved to %s\n', outDir);
```

#### Example Output

* **Console Output:**
    The script provides real-time feedback on the algorithm's stability and convergence progress:
    ```text
    =========FASC Start
    Iteration 1 start
    ...
    Iteration 22 start 
      Current cluster count: 50 
      No cluster merged 
      Outlier count: 10143
      Completed in 0.44s
      Inter iteration similarity: 99.99991%
    Converged!

    Total clock time: 1.199s
    CPU time: 3.862s

    Outliers count: 1043
    FASC End=========
    ```

* **Workspace Variables:**
    * `centers`: Matrix $(K \times D)$ containing the centroids of the $K$ identified clusters.
    * `counts`: Vector $(K \times 1)$ containing the number of particles in each cluster.
    * `labels`: Vector $(N \times 1)$ containing the cluster assignment ID for each particle ($0$ indicates an unassigned outlier).
    * `info`: Structure containing iteration history, convergence metrics, and CPU timing data.

* **Generated Figures (Saved in `results/`):**
    1.  **Convergence Plot:** A dual-axis graph showing the stabilization of "Inter-Iteration Similarity" and "Cluster Count" over time.
    2.  **Outlier Statistics:** A plot tracking the number of unassigned particles (outliers) per iteration.
    3.  **Similarity Heatmap:** A visual matrix showing the Dual-Cosine similarity between the final cluster centers.
    4.  **Cluster Distribution:** A bar chart displaying the fraction of total particles assigned to each cluster ID.
    

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

> **FASC: A Flexible Adaptive Stable Clustering Algorithm for Archive-Scale Mass Spectrometry**
> Shao Shi, et al.
> *Manuscript Submitted / Under Review*, 2025.
> (The full citation and DOI will be updated upon publication.)

## Contact

**Shao Shi (石邵)**
School of Environmental Science and Engineering
Southern University of Science and Technology (SUSTech)
Shenzhen, China
Email: [12231091@mail.sustech.edu.cn]
GitHub: [https://github.com/s129136908794904](https://github.com/s129136908794904)

***Prof. Xin Yang**
School of Environmental Science and Engineering
Southern University of Science and Technology (SUSTech)
Shenzhen, China
Email: [yangx@sustech.edu.cn]
