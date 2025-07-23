## TPP Explorer for Gene Drive Diagnostic Tools

This interactive tool allows users to explore test performance profiles for gene drive detection.

### Tabs Overview

1. **Sensitivity Explorer:**
   Plot of required sample size vs. test sensitivity for different pool sizes.

2. **Specificity Explorer:**
   Plot of required sample size vs. test specificity for different pool sizes.

3. **Precision Heatmap:**
   Heatmap of estimated standard error across sensitivity and specificity ranges for a fixed total sample and pool size.

4. **Pooling Design:**
   A simulation-based table of candidate pooling strategies, ranking designs by precision (standard deviation) of the prevalence estimate.

5. **Prevalence Estimator:**
   A simulation-based comparison of true vs. estimated prevalence across designs.



### How to Use

- Select a tab for the output you want to explore (Sensitivity, Specificity, Precision Heatmap, Pooling Design, or Prevalence Estimator).
- Adjust input parameters in the left‐hand sidebar. Tooltips and labels (with mathematical notation) guide you.
- Run the analysis by clicking the main action button (e.g., Run Simulation, Generate Plot) or by simply watching the live plots update for dynamic inputs.
- Reset to defaults at any time with the Default Parameter or Reset button in each panel.
- Interpret the results:
   * Outputs 1–3 give visual guidance on sample sizes or precision under varying test properties.
   * Output 4 helps you choose an optimal pooling strategy given resource constraints on sampled mosquitoes and total number of tests.
   * Output 5 shows how accurately each design recovers true prevalence across its full range.

---
