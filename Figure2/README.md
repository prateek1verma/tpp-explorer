# Figure 2 – Detection Accuracy and Pooling Effects

This folder contains R scripts used to generate the panels in Figure 1 of the gene drive monitoring manuscript. These figures explore how test sensitivity, specificity, and pooling strategies affect the probability of gene drive detection in mosquito populations.

## 📊 Panels

- **Figure 1A**: Required sample size vs. test sensitivity.
- **Figure 1B**: Required test specificity for controlling false positives.
- **Figure 1C**: Pooling effects on sensitivity.
- **Figure 1D**: Pooling effects on specificity.

## 🛠️ Files

- `Fig1_combine.R` generates data and combines individual panel plots (1A–1D) into a single multi-panel Figure 1 layout for the manuscript.
- `Fig1_Combined_MultiPanel.pdf` is the generated Figure in PDF format.
- `Fig1_Combined_MultiPanel.pdng` is the generated Figure in PNG format.


Each script generates a standalone figure using base R or ggplot2.

## 🧠 Methodology

Figures are based on analytic expressions derived in the manuscript, such as:

```math
x_{m} = 1 - (1 - x)^{m}
```

```math
p_{GD} = 1 - (1 - x_{m}s_{m})^{n/m}
```

```math
p_{FP} = 1 - (1 - f_{m})^{n/m}
```

These equations determine the probability of detection and false positives for given sample sizes.

## 📎 Dependencies

- `ggplot2`
- `dplyr`
