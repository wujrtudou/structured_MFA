# Module-Structured Mixture Factor Models (MS-MFM)

A statistical framework for discovering **disease subtypes and molecular signatures** from high-dimensional gene expression data.

---

## 📌 Overview

Modern transcriptomic datasets present several challenges:

* Extremely high dimensionality (p ≫ n)
* Strong gene–gene dependence
* Biological heterogeneity across samples

This project implements a **Module-Structured Mixture Factor Model (MS-MFM)** that addresses these issues by combining:

* Gene module structure (biological prior)
* Mixture models (clustering)
* Factor analysis (dimensionality reduction)

---

## 🧠 Core Idea

Instead of modeling each gene independently or using a full covariance matrix, we:

1. **Group genes into modules** (co-expression clusters)
2. **Model cluster differences at the module level**
3. **Capture dependencies using latent factors**

This results in a model that is both:

* ✅ Statistically efficient
* ✅ Biologically interpretable

---

## 🏗️ Model Formulation

Each sample ( x_i \in \mathbb{R}^p ) is modeled as:

[
x_i = \delta + M\alpha_k + B f_i + \epsilon_i
]

Where:

* ( \delta ): global gene-level baseline
* ( \alpha_k ): cluster-specific module effects
* ( f_i ): latent factors (low-dimensional)
* ( B ): structured loading matrix
* ( \epsilon_i ): gene-specific noise

---

## 🔍 Interpretation

The model decomposes variation into:

1. **Global effects** — baseline gene expression
2. **Module-level shifts** — disease/subtype differences
3. **Latent structure** — gene dependencies
4. **Residual noise** — measurement error

---

## ⚙️ Inference Algorithm

We use an **ECM (Expectation–Conditional Maximization)** algorithm:

### E-step

* Compute cluster responsibilities
* Estimate latent factor expectations

### CM-steps

* Update mixture weights
* Update module effects
* Update factor loadings
* Update gene-specific variances

---

## 📊 Data Processing Pipeline

### 1. Gene Filtering

* Select highly variable genes using **MAD (Median Absolute Deviation)**

### 2. Module Construction

* Hierarchical clustering
* Correlation-based distance
* Dendrogram cutting to define modules

### 3. Normalization

* Standardize gene expression (mean = 0, variance = 1)

---

## 📈 Model Selection

We choose:

* Number of clusters ( K )
* Latent dimension ( q )

Using:

👉 **Bayesian Information Criterion (BIC)**

---

## 🧪 Example Application

Applied to autoimmune diseases:

* Rheumatoid Arthritis (RA)
* Systemic Lupus Erythematosus (SLE)

### Key Results

* Identified **9 clusters (subtypes)**
* Achieved **clear disease separation without labels**
* Revealed **significant within-disease heterogeneity**
* Discovered **biologically meaningful module signatures**

---

## 💡 Why This Matters

Compared to traditional methods:

| Method  | Limitation                     |
| ------- | ------------------------------ |
| K-means | Ignores correlation            |
| GMM     | Unstable in high dimensions    |
| PCA     | Not interpretable biologically |

👉 **MS-MFM solves all three:**

* Handles high dimensions
* Models dependence
* Keeps biological interpretability

---

## 🚀 Getting Started

### Requirements

* Python / R (depending on implementation)
* Standard scientific libraries (numpy, scipy, sklearn)

### Basic Workflow

```bash
# 1. Preprocess data
# 2. Construct gene modules
# 3. Fit MS-MFM model
# 4. Select K and q via BIC
# 5. Analyze clusters and module effects
```

---

## 📂 Project Structure

```
├── data/               # Input gene expression data
├── preprocessing/      # Gene filtering & module construction
├── model/              # MS-MFM implementation
├── experiments/        # Case studies & results
├── utils/              # Helper functions
└── README.md
```

---

## 🔬 Applications

* Disease subtype discovery
* Biomarker identification
* Precision medicine
* Multi-omics integration (extendable)

---

## 📚 Reference

If you use this work, please cite:

```
Module-Structured Mixture Factor Models to Identify Outcome-Specific Signatures in Gene Expression Data
UQ Team, The University of Queensland
```

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repo
2. Create a feature branch
3. Submit a pull request

---

## 📜 License

MIT License (or specify your license here)

---

## ✨ Acknowledgements

* University of Queensland
* Open transcriptomic datasets (e.g., ADEx)

---

## 📬 Contact

For questions or collaborations:

* [Your Name / Email]
* GitHub Issues

---
