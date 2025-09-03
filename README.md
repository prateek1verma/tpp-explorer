# TPP Explorer – Gene Drive Diagnostic Tool

This repository hosts the code and data for the *TPP Explorer*, an R Shiny web application designed to support the development and assessment of Target Product Profiles (TPPs) for diagnostics detecting gene drive constructs in mosquito populations, and the codes for the plots used to generate panels in Figure 1 and 2 of the manuscript.

🧪 The app allows users to explore tradeoffs between test sensitivity, specificity, and sample size under different monitoring scenarios, as outlined in the accompanying research manuscript.

📍 **Live App**: [TPP Explorer on shinyapps.io](https://pverma.shinyapps.io/tpp_explorer/)

## 📁 Repository Structure

- `Rshinny_tpp_explorer/app.R` — Main R Shiny application script.
- `Rshinny_tpp_explorer/www/` — Folder containing assets such as logos or styling resources.
- `Figure1/` — Scripts and plots used to generate panels in Figure 1 of the manuscript.
- `Figure2/` — Scripts and plots used to generate panels in Figure 2 of the manuscript.

## 🧰 Features

- Computes minimum sample sizes for gene drive detection under different sensitivities/specificities.
- Visualizes false positive and false negative rates for both individual and pooled samples.
- Explores the standard error in estimating gene drive carrier frequency.
- Simulates and ranks candidate pooling designs by their precision (standard deviation), accounting for limits on mosquito samples and diagnostic tests.
- Simulates estimated prevalence across a range of true values to evaluate the performance of a given pooling design.

## 📖 Background

The methodology and motivation behind this tool are described in our working manuscript:
> *“Target Product Profiles for Gene Drive Monitoring”*

## 🚀 Getting Started (Local)

To run the app locally using RStudio:
1. Open the `app.R` file in RStudio.
2. Click the **Run App** button in the upper-right corner of the editor window.

Alternatively, you can run the app from the console with the following commands:

```R
install.packages("shiny")  # Install Shiny if not already installed
runApp("path/to/tpp-explorer")  # Replace with the actual path to the app directory

```

## 📜 License

This project is licensed under the MIT License.

## 👤 Author

Developed by [Prateek Verma](https://sites.google.com/view/prateekverma) and [​John M. Marshall](https://www.marshalllab.com/) .

