# Figure 3 – Estimating Carrier Frequency

This folder contains R scripts used to generate Figure 3 of the gene drive monitoring manuscript. This figure visualizes how test sensitivity and specificity influence the accuracy (standard error) of estimating gene drive carrier frequency.

## 📊 Panels

- **Figure 3**: Standard error of gene drive carrier frequency estimate under various sensitivity/specificity configurations.

## 🛠️ Files

| File                    | Description                                                                                                                                           |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Fig2_combined.m`       | Main script that generates Figure 3A and 3B showing standard error contours across sensitivity and specificity values for two gene drive frequencies. |
| `compute_sens_spec.m`   | Helper function used in `Fig2_combined.m` to compute required sensitivity when sensitivity = specificity for a target standard error.                 |
| `compute_sensitivity.m` | Computes required sensitivity given a fixed specificity to achieve a target standard error.                                                           |
| `brewermap.m`           | Utility script for generating colorblind-friendly colormaps used in the visualizations.                                                               |
| `Figure3_Combined.pdf`  | Exported final version of Figure 3 (A and B) as a publication-quality PDF.                                                                            |
| `Figure3_Combined.png`  | Exported PNG image of the final combined Figure 3 panels for quick preview or web use.                                                                |
| `README.md`             | Documentation describing the purpose and structure of the `Figure3` folder.                                                                           |

## 🧠 Methodology

This script computes and visualizes the standard error of a corrected estimator for gene drive carrier frequency:
```math
p̂ = (y - f) / (s - f)
```
where `s` is test sensitivity and `f` is the false positive rate.

## 📎 Dependencies

The following MATLAB files and functions are required to generate Figure 3:

- `Fig3_combined.m`: Main script that generates Figure 3A and 3B; depends on the helper functions below.
- `compute_sens_spec.m`: Computes sensitivity when sensitivity equals specificity for a desired standard error.
- `compute_sensitivity.m`: Computes required sensitivity given fixed specificity to achieve target precision.
- `brewermap.m`: Provides ColorBrewer color palettes for consistent and accessible figure styling.

**MATLAB built-ins used**: `fzero`, `meshgrid`, `sqrt`, `imagesc`, `contour`, `plot`, `colorbar`, `exportgraphics`, etc.

No additional MATLAB toolboxes are required beyond the standard distribution.

