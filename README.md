# 🧬 Gene Expression Analysis and Cancer Classification using Machine Learning

## 📘 Project Overview
This project analyzes gene expression data to distinguish between invasive and non-invasive cancer types. Using supervised and unsupervised machine learning techniques, the study performs dimensionality reduction, clustering, and classification to identify key gene markers and classify cancer samples accurately.

---

## 🗃️ Data Sources
- Gene expression dataset (`gene-expression-invasive-vs-noninvasive-cancer (1).csv`) with thousands of gene features.
- The dataset includes sample classifications into invasive (class 1) and non-invasive (class 0) cancer types.

---

## ⚙️ Methodology

### 1. Data Preprocessing
- Loaded and examined data structure and dimensions.
- Converted class labels (1 → 0 for non-invasive, 2 → 1 for invasive).
- Checked and imputed missing values using KNN imputation to clean data.
- Selected a subset of 2000 gene features for focused analysis.

### 2. Dimensionality Reduction
- Conducted two-sample t-tests for supervised feature selection, retaining genes with p-values < 0.05.
- Applied variance filtering to identify genes with variance > 0.25.
- Used PCA to summarize informative features and visualize variance explained.
- Applied t-SNE for nonlinear dimensionality reduction and visualization of sample groups.

### 3. Clustering Analysis
- Performed K-means clustering on both patient samples and gene sets to discover natural groupings.
- Utilized elbow method to select optimal cluster numbers.
- Visualized clusters with convex hulls using `factoextra`.
- Conducted hierarchical clustering with Euclidean distance and complete linkage, generating dendrograms.

### 4. Classification Modeling
- Split data into training (70%) and testing (30%) sets with controlled randomization.
- Implemented multiple classification algorithms:
  - Logistic Regression
  - Linear Discriminant Analysis (LDA)
  - Random Forests with hyperparameter tuning
  - Support Vector Machines (SVM) with radial kernel
  - Naive Bayes
  - K-Nearest Neighbors (KNN) with k selection via cross-validation
  - Quadratic Discriminant Analysis (QDA)
- Evaluated models using confusion matrices, accuracy, misclassification errors, and ROC metrics.
- Used K-fold cross-validation to assess model robustness and optimize hyperparameters.

---

## 📈 Key Results
- Identified key gene subsets with significant differential expression associated with cancer invasiveness.
- PCA and t-SNE visualization revealed clear separation between invasive and non-invasive groups.
- Random Forest and SVM models showed superior classification accuracy with tuned parameters.
- Hierarchical and K-means clustering captured meaningful groupings consistent with biological insights.
- Cross-validation confirmed model stability and predictive reliability.

---

## 💡 Practical Implications
- The study aids in understanding molecular differences in cancer types for better diagnosis.
- Demonstrates the utility of combining statistical tests, dimensionality reduction, clustering, and ML classification in bioinformatics.
- Provides a framework for integrating large-scale genomic data into clinical prognosis tools.

---

## 🚀 Future Work
- Incorporate additional clinical and genomic modalities to enrich models.
- Explore deep learning for feature extraction and classification improvements.
- Develop interpretable machine learning models for actionable insights in oncology.
- Validate findings on independent, external datasets for generalizability.

---


