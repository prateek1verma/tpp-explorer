# install.packages("patchwork")   # If not already installed
# install.packages("ggplot2")     # For safety, ensure ggplot2 is installed
library(ggplot2)
library(patchwork)
library(ggplot2)
library("ggsci")
library(extrafont)

###############################################################################
##### FIGURE 1A #####
###############################################################################
# Define parameters
x <- 0.05  # Gene drive carrier frequency (5%)
pGD_values <- c(0.95, 0.99)  # Desired probabilities of detection

# Function to calculate required sample size given sensitivity 's' and pGD
required_sample_size <- function(s, x = 0.05, pGD) {
  ifelse(s <= 0, NA, log(1 - pGD) / log(1 - x * s))
}

# Generate data for both pGD values
sensitivity_values <- seq(0.2, 1, by = 0.01)
plot_data <- do.call(rbind, lapply(pGD_values, function(p) {
  data.frame(
    Sensitivity = sensitivity_values * 100,  # Convert to percentage
    SampleSize = sapply(sensitivity_values, required_sample_size, x = x, pGD = p),
    pGD = paste0("pGD = ", p * 100, "%")
  )
}))

# Choose colors from ColorBrewer Set1 palette
# colors <- brewer.pal(n = 5, "Greys")[3:4]  # Select first two distinct colors
colors <- c("#404040", "#CCCCCC") # ,"#969696"

# Plot
fig1A_plot <- ggplot(plot_data, aes(x = Sensitivity, y = SampleSize, color = pGD)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(100, 200), linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = colors,
    name = "Detection Probability",  # You can set to NULL or "" to remove legend title
    labels = c(
      expression(italic(p[GD]) == 95 * "%"),
      expression(italic(p[GD]) == 99 * "%")
    )
  ) +
  scale_x_continuous(name = "Test Sensitivity per Mosquito (%)", breaks = seq(20, 100, 20)) +
  scale_y_continuous(name = "Required Sample Size", limits = c(50, 260)) +
  theme_bw(base_size = 18, base_family = "Times") +
  theme(
    text = element_text(family = "Times"),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = alpha('white', 0.8)),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank()
  )


###############################################################################
##### FIGURE 1B #####
###############################################################################

# Define parameters
pFP <- 0.05  # Acceptable false positive rate at the sample level

# Function to calculate required sample size given specificity 'sp'
required_sample_size_specificity <- function(sp, pFP = 0.05) {
  f <- 1 - sp  # False positive rate per mosquito
  ifelse(f <= 0, NA, (log(1 - pFP) / log(1 - f)))
}

# Generate data for plot
specificity_values <- seq(0.9990, 1, by = 0.00001)  # 99.90% to 100% specificity
sample_sizes <- sapply(specificity_values, required_sample_size_specificity)

plot_data_B <- data.frame(
  Specificity = specificity_values * 100,  # Convert to percentage
  SampleSize = sample_sizes
)

# Plot using ggplot2
fig1B_plot <- ggplot(plot_data_B, aes(x = Specificity, y = SampleSize)) +
  geom_line(size = 1.2, color = "#404040") +  # Single curve, consistent black color
  geom_hline(yintercept = c(100, 200), linetype = "dashed", color = "grey40") +
  # annotate("text", x = 99.92, y = 107, label = "Sample Size = 100", hjust = 0, size = 5, family = "Times") +
  # annotate("text", x = 99.92, y = 207, label = "Sample Size = 200", hjust = 0, size = 5, family = "Times") +
  scale_x_continuous(name = "Test Specificity per Mosquito (%)", breaks = seq(99.90, 100, 0.02)) +
  scale_y_continuous(name = "Required Sample Size", limits = c(50, 260)) +
  theme_bw(base_size = 18, base_family = "Times") +
  theme(
    text = element_text(family = "Times"),
    legend.position = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.box.just = c("top"),
    legend.background = element_rect(fill = alpha('white', 0.8)),
    # axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank()
  )

###############################################################################
##### FIGURE 1C #####
###############################################################################
# Load necessary library
font_import()  # Run once to import system fonts
loadfonts(device = "pdf")  # Load fonts for PDF output

# Define parameters
x <- 0.05                 # Gene drive carrier frequency (5%)
pGD0 <- 0.95               # Desired detection probability
pool_sizes <- c(1, 5, 10, 20)  # Pool sizes to compare
sensitivity_values <- seq(0.2, 1, by = 0.01)  # Sensitivity range (40% to 100%)

# Function to calculate required sample size for pooled testing
required_sample_size_pooled <- function(s, m, x = 0.05, pGD = pGD0) {
  x_m <- 1 - (1 - x)^m  # Probability that at least one gene drive mosquito is in the pool
  if (x_m * s >= 1) return(NA)  # Avoid log(0) or negative values
  n <-  m *(log(1 - pGD) / log(1 - x_m * s))
  return(n)
}

# Create the plotting data
plot_data_C <- expand.grid(Sensitivity = sensitivity_values, PoolSize = pool_sizes)
plot_data_C$SampleSize <- mapply(required_sample_size_pooled, 
                                 s = plot_data_C$Sensitivity, 
                                 m = plot_data_C$PoolSize)

# Plot
fig1C_plot <- ggplot(plot_data_C, aes(x = Sensitivity * 100, y = SampleSize, color = factor(PoolSize))) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(100, 200), linetype = "dashed", color = "grey40") +
  # scale_color_npg(name = "Pool Size (m)") +  # Automatic NPG color palette
  scale_color_brewer(palette = "RdGy", direction = -1, name = "Pool Size (m)") +
  # scale_color_manual(name = "Pool Size (m)", values = c("black","#1f77b4", "#d62728", "#2ca02c")) +
  scale_x_continuous(expression("Test Sensitivity per Pool of " * italic(m) * " Mosquitoes (%)"), breaks = seq(20, 100, 20)) +
  scale_y_continuous(name = "Required Sample Size", limits = c(50, 260)) +
  theme_bw(base_size = 18, base_family = "Times") +
  theme(
    # legend.position = "inside",
    # legend.position.inside = c(0.85,0.8),
    text = element_text(family = "Times"),
    legend.position=c(0.99,0.65),legend.justification=c(1,1),
    legend.direction="vertical",
    legend.box="horizontal",
    legend.box.just = c("top"), 
    legend.background = element_rect(fill=alpha('white', 0.7)),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank()
  )

###############################################################################
##### FIGURE 1D #####
###############################################################################
# Define parameters
pFP <- 0.05                 # Acceptable sample-wide false positive rate (5%)
pool_sizes <- c(1, 5, 10, 20)  # Pool sizes to compare
specificity_values <- seq(0.994, 1, by = 0.00005)  # Specificity range (99.9% to 100%)

# Function to calculate required sample size for pooled testing based on specificity
required_sample_size_specificity <- function(sp, m, pFP = 0.05) {
  f <- 1 - sp  # False positive rate per mosquito
  f_m <- 1 - (1 - f)^m
  if (f >= 1) return(NA)  # Avoid log(0) or invalid values
  n <- m *(log(1 - pFP) / log(1 - f)) # assumes n/m to be integer
  # n <- ceiling(m * log(1 - pFP) / log(1 - f))
  return(n)
}

# Create the plotting data
plot_data_D <- expand.grid(Specificity = specificity_values, PoolSize = pool_sizes)
plot_data_D$SampleSize <- mapply(required_sample_size_specificity, 
                                 sp = plot_data_D$Specificity, 
                                 m = plot_data_D$PoolSize)

# Plot
# Plot using ggplot2 with updated styling matching Fig 1C
fig1D_plot <- ggplot(plot_data_D, aes(x = Specificity * 100, y = SampleSize, color = factor(PoolSize))) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(100, 200), linetype = "dashed", color = "grey40") +
  # scale_color_npg(name = "Pool Size (m)") +  # Use NPG color palette automatically
  scale_color_brewer(palette = "RdGy", direction = -1,name = "Pool Size (m)") +
  scale_x_continuous(expression("Test Specificity per Pool of " * italic(m) * " Mosquitoes (%)"), breaks = seq(99.4, 100, 0.1)) +
  scale_y_continuous(name = "Required Sample Size", limits = c(50, 260)) +
  theme_bw(base_size = 18, base_family = "Times") +
  theme(
    text = element_text(family = "Times"),
    legend.position = "none",  # Remove legend
    # legend.position = c(0.99, 0.99),
    # legend.justification = c(1, 1),
    # legend.direction = "vertical",
    # legend.box = "horizontal",
    # legend.box.just = c("top"),
    # legend.background = element_rect(fill = alpha('white', 0.6)),
    # axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_blank()  # No plot title as requested
  )


###############################################################################
##### COMBINE ALL FIGURES #####
###############################################################################

# Create combined figure with annotation tags A, B, C, D
combined_figure <- (fig1A_plot + fig1B_plot) / (fig1C_plot + fig1D_plot) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20, face = "bold"))

# Display the combined figure
print(combined_figure)

# Save in high resolution
ggsave("Fig1_Combined_MultiPanel.pdf", plot = combined_figure, width = 12, height = 10, units = "in")
ggsave("Fig1_Combined_MultiPanel.png", plot = combined_figure, width = 12, height = 10, units = "in", dpi = 600, bg = "white")
