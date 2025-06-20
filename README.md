# Genomic-Variant-Predictor

This repository contains a complete machine learning pipeline to classify human genetic variants (from the ClinVar database) as **Pathogenic** or **Benign** using biological features and mutation annotations.

🔍 The project combines:
- Data preprocessing & annotation
- Mutation classification
- Statistical testing (Chi-Square)
- Machine learning (XGBoost, Logistic Regression, Random Forest)
- Hyperparameter tuning
- Model explainability (SHAP)

---

## Dataset

- **Source**: [ClinVar Variant Summary](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
- Filtered for:
  - Human variants
  - Only `Benign` and `Pathogenic`
  - Reviewed (with clear `ReviewStatus`)
- Mutation types extracted: missense, nonsense, frameshift, splice-site, etc.

---

##  Methods

### Mutation Classification
- Regex-based parsing of the `Name` column to classify:
  - Missense
  - Nonsense
  - Frameshift
  - Splice-site
  - Others (e.g., synonymous, intronic)

### Statistical Analysis
- Chi-Square test showed significant association between **mutation type** and **pathogenicity** (p < 0.05)

### ML Models Trained
- **XGBoost (Best performing)**
- Random Forest
- Logistic Regression

### Hyperparameter Tuning
- Performed using `RandomizedSearchCV` and `GridSearchCV` to optimize XGBoost parameters like `max_depth`, `learning_rate`, `n_estimators`, etc.

---

## Final Results (XGBoost)

| Metric                | Value     |
|-----------------------|-----------|
| **Accuracy**          | 86%       |
| **ROC-AUC**           | `0.936`   |
| **Precision (Pathogenic)** | 0.89 |
| **Recall (Pathogenic)**    | 0.81 |
| **F1-score (Pathogenic)**  | 0.85 |
| **Total Samples**     | 155,185   |

### Confusion Matrix
|                            | **Predicted: Benign (0)** | **Predicted: Pathogenic (1)** |
| -------------------------- | ------------------------- | ----------------------------- |
| **Actual: Benign (0)**     | 76,085                    | 7,465                         |
| **Actual: Pathogenic (1)** | 13,582                    | 58,053                        |


## Model Explainability (SHAP)

We used SHAP (SHapley Additive exPlanations) to interpret the XGBoost model and identify which features most influenced predictions:

- **Top Features Identified**:
  - `mutation_type`
  - `OriginSimple`
  - `Chromosome`
  - `is_missense`, `is_nonsense`, etc.

## Visuals Included

- Bar plots of mutation types
- Pathogenic vs Benign distribution
- Heatmap of mutation type vs label
- SHAP summary and feature importance plots
- ROC curves and confusion matrix

