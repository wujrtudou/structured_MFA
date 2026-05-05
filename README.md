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
