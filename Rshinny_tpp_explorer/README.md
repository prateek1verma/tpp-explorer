## Gene Drive TPP Explorer

This Shiny app helps researchers and developers explore diagnostic requirements for detecting gene drives in mosquito populations. It visualizes relationships between sensitivity, specificity, pool sizes, and required sample sizes based on Target Product Profiles (TPPs) for two use cases:

1. **Presence Detection** — Sample size required to detect gene drive presence at low frequency.
2. **Carrier Frequency Estimation** — Precision (standard error) in estimating gene drive carrier frequency.

📍 **Live App**: [TPP Explorer on shinyapps.io](https://pverma.shinyapps.io/tpp_explorer/)

### Features
- Dynamic sample size vs. sensitivity/specificity plots
- Pooled sample analysis with customizable pool sizes
- Interactive precision heatmap for diagnostic accuracy
- Reset buttons for quick parameter reconfiguration
- Markdown tutorial for onboarding new users

### Technologies
- R Shiny
- `ggplot2` and `image` for interactive visualization
- Bootstrap-based responsive UI

### Getting Started
To run locally:

```r
shiny::runApp()
