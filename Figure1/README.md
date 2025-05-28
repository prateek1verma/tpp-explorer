# Figure 1 – Detection Accuracy and Pooling Effects

This folder contains R scripts used to generate the panels in Figure 1 of the gene drive monitoring manuscript. These figures explore how test sensitivity, specificity, and pooling strategies affect the probability of gene drive detection in mosquito populations.

## 📊 Panels

- **Figure 1A**: Required sample size vs. test sensitivity.
- **Figure 1B**: Required test specificity for controlling false positives.
- **Figure 1C**: Pooling effects on sensitivity.
- **Figure 1D**: Pooling effects on specificity.

## 🛠️ Files

- `fig1A_sample_size_vs_sensitivity.R`
- `fig1B_false_positive_vs_specificity.R`
- `fig1C_pooled_sensitivity.R`
- `fig1D_pooled_specificity.R`

Each script generates a standalone figure using base R or ggplot2.

## 🧠 Methodology

Figures are based on analytic expressions derived in the manuscript, such as:

```math
p_{GD} = 1 - (1 - xs)^n
p_{FP} = 1 - (1 - f)^n
```

These equations determine the probability of detection and false positives for given sample sizes.

## 📎 Dependencies

- `ggplot2`
- `dplyr`
