# Figure 2 – Estimating Carrier Frequency

This folder contains R scripts used to generate Figure 2 of the gene drive monitoring manuscript. This figure visualizes how test sensitivity and specificity influence the accuracy (standard error) of estimating gene drive carrier frequency.

## 📊 Panels

- **Figure 2**: Standard error of gene drive carrier frequency estimate under various sensitivity/specificity configurations.

## 🛠️ Files

| File                    | Description                                                                                                                                           |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Fig2_combined.m`       | Main script that generates Figure 2A and 2B showing standard error contours across sensitivity and specificity values for two gene drive frequencies. |
| `compute_sens_spec.m`   | Helper function used in `Fig2_combined.m` to compute required sensitivity when sensitivity = specificity for a target standard error.                 |
| `compute_sensitivity.m` | Computes required sensitivity given a fixed specificity to achieve a target standard error.                                                           |
| `brewermap.m`           | Utility script for generating colorblind-friendly colormaps used in the visualizations.                                                               |
| `Figure2_Combined.pdf`  | Exported final version of Figure 2 (A and B) as a publication-quality PDF.                                                                            |
| `Figure2_Combined.png`  | Exported PNG image of the final combined Figure 2 panels for quick preview or web use.                                                                |
| `README.md`             | Documentation describing the purpose and structure of the `Figure2` folder.                                                                           |


This script computes and visualizes the standard error of a corrected estimator for gene drive carrier frequency:
```math
p̂ = (y - f) / (s - f)
```
where `s` is test sensitivity and `f` is the false positive rate.

## 📎 Dependencies

- `ggplot2`
- `dplyr`
