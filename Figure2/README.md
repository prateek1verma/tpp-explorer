# Figure 2 – Estimating Carrier Frequency

This folder contains R scripts used to generate Figure 2 of the gene drive monitoring manuscript. This figure visualizes how test sensitivity and specificity influence the accuracy (standard error) of estimating gene drive carrier frequency.

## 📊 Panels

- **Figure 2**: Standard error of gene drive carrier frequency estimate under various sensitivity/specificity configurations.

## 🛠️ Files

- `fig2_standard_error_estimation.R`

This script computes and visualizes the standard error of a corrected estimator for gene drive carrier frequency:
```math
p̂ = (y - f) / (s - f)
```
where `s` is test sensitivity and `f` is the false positive rate.

## 📎 Dependencies

- `ggplot2`
- `dplyr`
