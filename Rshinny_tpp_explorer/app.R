library(shiny)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(bslib)
library(shinyjs)
library(thematic)
library(markdown)
library(DT)

# --- Helper Functions ---
simulate_y_sum <- function(p_true, m_vec, k_vec, sens_vec, spec_vec, n_rep) {
  n_groups <- length(m_vec)
  x_m <- 1 - (1 - p_true)^m_vec
  theta <- sens_vec * x_m + (1 - spec_vec) * (1 - x_m)
  
  # Simulate n_rep binomial draws per group
  y_mat <- replicate(n_rep, {
    rbinom(n_groups, size = k_vec, prob = theta)
  })
  
  # Sum the total number of positives across replicates for each group
  y_sum <- rowSums(y_mat)  # total positive pools over all replicates
  return(y_sum)
}


log_likelihood_sum <- function(p, m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep) {
  if (p <= 0 || p >= 1) return(-Inf)
  # Compute expected positive test probability for each group
  x_m <- 1 - (1 - p)^m_vec
  theta <- sens_vec * x_m + (1 - spec_vec) * (1 - x_m)
  theta <- pmin(pmax(theta, 1e-10), 1 - 1e-10)  # avoid log(0)
  K_vec <- n_rep * k_vec # Total number of trials per group
  ll_terms <- y_vec_sum * log(theta) + (K_vec - y_vec_sum) * log(1 - theta)
  sum(ll_terms)
}


estimate_p_and_se <- function(m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep) {
  # Negative log-likelihood for optimization
  neg_logL <- function(p) -log_likelihood_sum(p, m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep)
  
  # Optimize to find MLE
  opt <- tryCatch({
    optimize(neg_logL, interval = c(0.0001, 0.9999))
  }, error = function(e) return(list(minimum = NA, objective = NA)))
  
  if (is.null(opt$minimum) || is.na(opt$minimum)) return(c(p_hat = NA, se = NA))
  p_hat <- opt$minimum
  
  # Estimate second derivative of log-likelihood at p_hat
  h <- 1e-5
  logL_plus <- log_likelihood_sum(p_hat + h, m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep)
  logL_minus <- log_likelihood_sum(p_hat - h, m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep)
  logL_0 <- log_likelihood_sum(p_hat, m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep)
  
  second_deriv <- (logL_plus - 2 * logL_0 + logL_minus) / (h^2)
  
  # Standard error: SE = sqrt(1 / observed information)
  if (second_deriv >= 0 || !is.finite(second_deriv)) {
    se <- NA  # curvature is non-negative or invalid; can't compute SE
  } else {
    se <- sqrt(-1 / second_deriv)
  }
  
  return(c(p_hat = p_hat, se = se))
}


parse_input <- function(x) as.numeric(unlist(strsplit(x, ",")))

ui <- page_fluid(  # replaces fluidPage
  # theme = bs_theme(version = 4, bootswatch = "flatly"),  # default light theme
  # theme = bs_theme(version = 4, bootswatch = "cosmo"),
  # theme = bs_theme(version = 4, bootswatch = "lumen"), 
  theme = bs_theme(version = 4, bootswatch = "yeti"), 
  tags$head(
    tags$title("TPP Explorer | Gene Drive"),
    tags$link(rel = "icon", type = "image/png", href = "favicon.png")  # or favicon.ico
  ),
  titlePanel(
    div(style = "display: flex; justify-content: space-between; align-items: center;",
        div(style = "display: flex; align-items: center;",
            img(src = "logo.png", height = "58px", style = "margin-right: 12px;"),
            div(
              span("Target Product Profile (TPP) Explorer", style = "font-size: 28px; font-weight: 600; color: #005B96;"),
              br(),
              span("For Gene Drive Diagnostic Tools", style = "font-size: 16px; color: #777;")
            )
        )
        # Optional toggle or dark mode checkbox can go here
    )
  ),
  div(style = "margin-top: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 6px; background-color: #ffffff;", 
      tabsetPanel(type = "pills",
                  tabPanel(title = tagList(icon("search-plus"),"Sensitivity Explorer"),
                           sidebarLayout(
                             sidebarPanel(
                               h5("Sampling Parameters"),
                               numericInput("x1",HTML( "Gene Drive Frequency, <i>x</i>"), value = 0.05, min = 0, max = 1,step = 0.01),
                               numericInput("pGD0", HTML("Desired Detection Probability, <i>p<sub>GD</sub></i>"), value = 0.95, min = 0, max = 1,step = 0.01),
                               h5("Pooling Options"),
                               textInput("pool_sizes1", HTML("Pool Sizes, <i>m</i> (comma-separated)"), value = "1,5,10,20"),
                               h5("Plot Settings"),
                               sliderInput("xlim1", "X-axis limits for Sensitivity (%)", min = 0, max = 100, value = c(20, 100),step = 0.05),
                               sliderInput("ylim1", "Y-axis limits for Sample Size", min = 0, max = 500, value = c(50, 260),step = 1),
                               # actionButton("reset1", "Reset Values"),
                               tags$hr(),
                               div(style = "margin-top: 10px;", 
                                   actionButton(
                                     inputId = "reset1", 
                                     label = "Default Parameter", 
                                     icon = icon("rotate-left"), 
                                     style = "margin-top:10px; font-size: 13px; padding: 6px 12px;",
                                     title = "Reset parameters to default values"
                                   ))
                             ),
                             mainPanel(plotOutput("fig1C",width = "600px", height = "400px"))
                           )
                  ),
                  tabPanel(
                    title = tagList(icon("bullseye"), "Specificity Explorer"),
                    sidebarLayout(
                      sidebarPanel(
                        h5("Sampling Parameters"),
                        numericInput("pFP", HTML("Acceptable Sample-Wide False Positive Rate, <i>p<sub>FP</sub></i>"), value = 0.05, min = 0, max = 1,step = 0.01),
                        textInput("pool_sizes2", HTML("Pool Sizes, <i>m</i> (comma-separated)"), value = "1,5,10,20"),
                        h5("Plot Settings"),
                        sliderInput("xlim2", "X-axis limits for Specificity (%)", min = 95, max = 100, value = c(99.4, 100), step = 0.01),
                        sliderInput("ylim2", "Y-axis limits for Sample Size", min = 0, max = 500, value = c(50, 260),step = 1,title("Adjust the Y-axis limits")),
                        tags$hr(),
                        div(style = "margin-top: 10px;", 
                            actionButton(
                              inputId = "reset2", 
                              label = "Default Parameter", 
                              icon = icon("rotate-left"), 
                              style = "margin-top:10px; font-size: 13px; padding: 6px 12px;", # "padding:4px 8px; font-size:12px; background-color:#f0f0f0; color:#333; border:1px solid #ccc; border-radius:4px;",
                              title = "Reset parameters to default values"
                            ))
                      ),
                      mainPanel(plotOutput("fig1D",width = "600px", height = "400px"))
                    )
                  ),
                  tabPanel(
                    title = tagList(icon("th-large"), "Precision Heatmap"),
                    sidebarLayout(
                      sidebarPanel(
                        h5("Sampling Parameters"),
                        numericInput("y_true", HTML("Estimated Gene Drive Frequency, <i>p&#770;</i>"), value = 0.5,step = 0.01),
                        numericInput("n", HTML("Pooled Sample Count, <i>n</i>"), value = 100, step = 1),
                        numericInput("m", HTML("Pool Size, <i>m</i>"), value = 1, step = 1),
                        textOutput("pool_count"),  # <-- New line here
                        h5("Plot Settings"),
                        sliderInput("xy_lim", "Sensitivity/Specificity Axis Limits (%)", 
                                    min = 50.1, max = 100, value = c(60, 100), step = 1),
                        tags$hr(),
                        div(style = "margin-top: 10px;", 
                            actionButton(
                              inputId = "reset3", 
                              label = "Default Parameter", 
                              icon = icon("rotate-left"), 
                              style = "margin-top:10px; font-size: 13px; padding: 6px 12px;", # "padding:4px 8px; font-size:12px; background-color:#f0f0f0; color:#333; border:1px solid #ccc; border-radius:4px;",
                              title = "Reset parameters to default values"
                            ))
                        
                      ),
                      mainPanel(plotOutput("heatmapPlot",width = "500px", height = "500px"))
                    )
                  ),
                  tabPanel(title = tagList(icon("cogs"), "Pooling design"),
                           sidebarLayout(
                             sidebarPanel(
                               numericInput("N", HTML("Max Mosquito Samples, <i>N</i> :"), value = 100, step = 1),
                               numericInput("max_tests", HTML("Total Pool Count (kits), <i>K</i> :"), value = 10, step = 1),
                               textInput("m_vals", HTML("Candidate Pool Sizes, <i>m<sub>i</sub></i> (comma-separated) :"), value = "1,5,20"),
                               textInput("sens_vec", HTML("Test Sensitivity, <i>s<sub>i</sub></i> (comma-separated) :"), value = "0.9,0.8,0.7"),
                               textInput("spec_vec", HTML("Test Specificity, <i>c<sub>i</sub></i> (comma-separated) :"), value = "0.9,0.8,0.7"),
                               numericInput("p_true", HTML("True GD Frequency, <i>p</i> :"), value = 0.1, min = 0.001, max = 0.999, step = 0.001),
                               # numericInput("sensitivity", "Test Sensitivity :", value = 0.95, min = 0.5, max = 1, step = 0.01),
                               # numericInput("specificity", "Test Specificity :", value = 0.99, min = 0.5, max = 1, step = 0.01),
                               numericInput("n_rep", HTML("No. of Simulations per Design, <i>n<sub>rep</sub></i> :"), value = 100),
                               actionButton("runSim", "Run Simulation"),
                               tags$hr(),
                               div(style = "margin-top: 10px;", 
                                   actionButton(
                                     inputId = "reset4", 
                                     label = "Default Parameter", 
                                     icon = icon("rotate-left"), 
                                     style = "margin-top:10px; font-size: 13px; padding: 6px 12px;", # "padding:4px 8px; font-size:12px; background-color:#f0f0f0; color:#333; border:1px solid #ccc; border-radius:4px;",
                                     title = "Reset parameters to default values"
                                   ))
                             ),
                             mainPanel(DTOutput("results_table"))
                           )
                  ),
                  tabPanel(
                    title = tagList(icon("book-open"), "Tutorial / Help"),
                    includeMarkdown("tutorial.md")
                  )
      )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$reset1, {
    updateNumericInput(session, "x1", value = 0.05)
    updateNumericInput(session, "pGD0", value = 0.95)
    updateTextInput(session, "pool_sizes1", value = "1,5,10,20")
    updateSliderInput(session, "xlim1", value = c(20, 100))
    updateSliderInput(session, "ylim1", value = c(50, 260))
  })
  
  observeEvent(input$reset2, {
    updateNumericInput(session, "pFP", value = 0.05)
    updateTextInput(session, "pool_sizes2", value = "1,5,10,20")
    updateSliderInput(session, "xlim2", value = c(99.4, 100))
    updateSliderInput(session, "ylim2", value = c(50, 260))
  })
  
  observeEvent(input$reset3, {
    updateNumericInput(session, "y_true", value = 0.5)
    updateNumericInput(session, "n", value = 100)
    updateNumericInput(session, "m", value = 1)
    updateSliderInput(session, "xy_lim", value = c(60, 100))
  })
  
  observeEvent(input$reset4, {
    updateNumericInput(session, "N", value = 100)
    updateNumericInput(session, "max_tests", value = 10)
    updateTextInput(session, "m_vals", value = "1,5,20")
    updateNumericInput(session, "p_true", value = 0.1)
    updateNumericInput(session, "sens_vec", value = "1,1,1")
    updateNumericInput(session, "spec_vec", value = "1,1,1")
    updateNumericInput(session, "n_rep", value = 50)
  })
  
  observeEvent(input$runSim, {
    m_vals <- parse_input(input$m_vals)
    N <- input$N
    max_tests <- input$max_tests
    sens_vec <- parse_input(input$sens_vec)
    spec_vec <- parse_input(input$spec_vec)
    k_vals <- 0:input$max_tests
    p_true <- input$p_true
    n_rep <- input$n_rep
    
    
    n <- length(m_vals)
    k_combinations <- expand.grid(replicate(n, k_vals, simplify = FALSE))
    results_list <- vector("list", nrow(k_combinations))
    row_idx <- 1
    
    if (any(is.na(m_vals)) || any(is.na(sens_vec)) || any(is.na(spec_vec))) {
      showNotification("Please ensure all inputs are valid numbers.", type = "error")
      return(NULL)
    }
    
    
    for (i in seq_len(nrow(k_combinations))) {
      k_vec <- as.integer(k_combinations[i, ])
      K <- sum(k_vec)
      N_total <- sum(k_vec * m_vals)
      
      if (K != input$max_tests || N_total > input$N) next
      m_vec <- m_vals
      y_vec_sum <- simulate_y_sum(p_true, m_vec, k_vec, sens_vec, spec_vec, n_rep)
      est <- estimate_p_and_se(m_vec, k_vec, y_vec_sum, sens_vec, spec_vec, n_rep)
      p_hats <- est["p_hat"]
      se_phat <- est["se"]
      
      results_list[[row_idx]] <- data.frame(
        t(setNames(as.list(rbind(m_vec, k_vec)), c(rbind(paste0("m", 1:n), paste0("k", 1:n))))),
        K = K, N = N_total,
        p_hat = p_hats,
        SE = se_phat
      )
      row_idx <- row_idx + 1
    }
    results <- do.call(rbind, results_list[1:(row_idx - 1)])
    
    
    # --- Reference Design ---
    m_vec_ref <- rep(1,length(m_vals))
    k_vec_ref <- rep(1,length(m_vals))
    k_vec_ref[1] <- input$N-length(m_vals)+1
    p_hats_ref <- numeric(0)
    
    y_vec_ref <- simulate_y_sum(p_true, m_vec_ref, k_vec_ref, sens_vec, spec_vec, n_rep)
    est_ref <- estimate_p_and_se(m_vec_ref, k_vec_ref, y_vec_ref, sens_vec, spec_vec, n_rep)
    p_hats_ref <- est_ref["p_hat"]
    se_hat_ref <-  est_ref["se"]
    
    # Create interleaved m1, k1, m2, k2, ...
    col_names <- as.vector(rbind(paste0("m", 1:n), paste0("k", 1:n)))
    ref_values <- as.vector(rbind(m_vec_ref, k_vec_ref))
    
    # Create 1-row data frame with atomic columns
    reference_row <- data.frame(
      matrix(as.numeric(ref_values), nrow = 1, dimnames = list(NULL, col_names)),
      K = input$N,
      N = input$N,
      p_hat = as.numeric(p_hats_ref),
      SE = as.numeric(se_hat_ref),
      stringsAsFactors = FALSE
    )
    
    
    output$results_table <- renderDT({
      # sorted_results <- results[order(as.numeric(results$SE)), ]
      results <- results[order(as.numeric(results$SE)), ]
      display <- head(results, 50)
      
      # Define the insertion position
      insert_pos <- 10
      
      common_cols <- intersect(names(reference_row), names(display))
      reference_row <- reference_row[, common_cols]
      display <- display[, common_cols]
      
      
      if (!is.data.frame(display)) {
        showNotification("Display is not a data frame", type = "error")
        return(NULL)
      }
      
      display_n <- nrow(display)
      
      if (display_n < insert_pos) {
        display <- rbind(display, reference_row)
      } else if (insert_pos <= 1) {
        display <- rbind(reference_row, display)
      } else {
        display <- rbind(
          display[1:(insert_pos - 1), ],
          reference_row,
          display[insert_pos:display_n, ]
        )
      }
      
      m_cols <- grep("^m[0-9]+$", names(display), value = TRUE)
      k_cols <- gsub("^m", "k", m_cols)
      
      ref_k_vals <- rep(1, length(k_cols))
      ref_k_vals[1] <- input$N - length(k_cols) + 1
      ref_m_vals <- rep(1, length(m_cols))
      
      is_ref_row <- rep(TRUE, nrow(display))
      for (i in seq_along(m_cols)) {
        is_ref_row <- is_ref_row & display[[m_cols[i]]] == ref_m_vals[i]
        is_ref_row <- is_ref_row & display[[k_cols[i]]] == ref_k_vals[i]
      }
      ref_row <- which(is_ref_row)
      
      
      # Create a logical vector to exclude reference row
      non_ref_indices <- setdiff(seq_len(nrow(display)), ref_row)
      
      # Rank only non-reference rows
      numeric_se <- as.numeric(display$SE)
      ranks <- rep(NA, length(numeric_se))
      ranks[non_ref_indices] <- rank(numeric_se[non_ref_indices], ties.method = "first")
      
      # Assign rank labels
      # rank_labels <- rep("", length(numeric_se))
      rank_labels <- character(nrow(display))
      rank_labels[ref_row] <- "📌 Reference"
      rank_labels[non_ref_indices] <- paste("Rank", ranks[non_ref_indices])
      rank_labels[ranks == 1] <- "🥇 Rank 1"
      rank_labels[ranks == 2] <- "🥈 Rank 2"
      rank_labels[ranks == 3] <- "🥉 Rank 3"
      
      # Reference label
      if (length(ref_row) == 1) {
        rank_labels[ref_row] <- "📌 Reference"
      }
      
      # Assign rank labels to display
      display$Rank <- rank_labels
      
      # Format numeric columns
      display$p_hat <- formatC(display$p_hat, format = "f", digits = 5)
      display$SE <- formatC(display$SE, format = "f", digits = 5)
      
      # Highlight only if Rank matches specific values in the row content
      js_callback <- JS("
    function(row, data, index) {
      if (data[0] === 'Rank 1') {
        $('td', row).css('background-color', '#76d69c');
      }
      if (data[0] === 'Rank 2') {
        $('td', row).css('background-color', '#a9e8bf');
      }
      if (data[0] === 'Rank 3') {
        $('td', row).css('background-color', '#d4f4dd');
      }
      if (data[0] === 'Reference') {
        $('td', row).css('background-color', 'lightyellow');
        $('td', row).css('font-weight', 'bold');
      }
    }
  ")
      
      datatable(
        display[, c("Rank", setdiff(names(display), "Rank"))],
        rownames = FALSE,
        options = list(pageLength = 10)
      ) %>%
        formatStyle(
          'Rank',
          target = 'cell',
          backgroundColor = styleEqual(
            c("🥇 Rank 1", "🥈 Rank 2", "🥉 Rank 3", "📌 Reference"),
            c("#76d69c", "#a9e8bf", "#d4f4dd", "lightyellow")
          ),
          fontWeight = styleEqual("📌 Reference", "bold")
        )
      
    })
  })
  
  
  
  thematic::thematic_shiny()
  
  output$fig1C <- renderPlot({
    x <- input$x1
    pGD0 <- input$pGD0
    pool_sizes <- as.numeric(strsplit(input$pool_sizes1, ",")[[1]])
    sens_range <- seq(input$xlim1[1] / 100, input$xlim1[2] / 100, by = 0.01)
    
    
    required_sample_size_pooled <- function(s, m) {
      x_m <- 1 - (1 - x)^m
      if (x_m * s >= 1) return(NA)
      m * (log(1 - pGD0) / log(1 - x_m * s))
    }
    
    df <- expand.grid(Sensitivity = sens_range, PoolSize = pool_sizes)
    df$SampleSize <- mapply(required_sample_size_pooled, df$Sensitivity, df$PoolSize)
    
    df <- na.omit(df)
    df <- df[df$SampleSize >= input$ylim1[1] & df$SampleSize <= input$ylim1[2], ]
    
    
    ggplot(df, aes(x = Sensitivity * 100, y = SampleSize, color = factor(PoolSize))) +
      geom_line(linewidth = 1.2) +
      scale_x_continuous(expression("Test Sensitivity per Pool of " * italic(m) * " Mosquitoes (%)"), limits = input$xlim1) +
      scale_y_continuous(name = "Required Sample Size", limits = input$ylim1) +
      scale_color_brewer(palette = "Set1", name = "Pool Size (m)") +
      theme_bw(base_size = 16) + 
      theme(
        axis.text = element_text(size = 14),        # Increase tick label size
        axis.title = element_text(size = 16),       # Axis title size (already covered by base_size)
        legend.text = element_text(size = 12),      # Optional: legend text
        legend.title = element_text(size = 14)      # Optional: legend title
      )
  })
  
  output$fig1D <- renderPlot({
    pFP <- input$pFP
    pool_sizes <- as.numeric(strsplit(input$pool_sizes2, ",")[[1]])
    spec_range <- seq(input$xlim2[1] / 100, input$xlim2[2] / 100, by = 0.00005)
    
    required_sample_size_specificity <- function(sp, m) {
      f <- 1 - sp
      if (f >= 1) return(NA)
      m * (log(1 - pFP) / log(1 - f))
    }
    
    df <- expand.grid(Specificity = spec_range, PoolSize = pool_sizes)
    df$SampleSize <- mapply(required_sample_size_specificity, df$Specificity, df$PoolSize)
    
    df <- na.omit(df)
    df <- df[df$SampleSize >= input$ylim2[1] & df$SampleSize <= input$ylim2[2], ]
    
    
    ggplot(df, aes(x = Specificity * 100, y = SampleSize, color = factor(PoolSize))) +
      geom_line(linewidth = 1.2) +
      scale_x_continuous(expression("Test Specificity per Pool of " * italic(m) * " Mosquitoes (%)"),, limits = input$xlim2) +
      scale_y_continuous(name = "Required Sample Size", limits = input$ylim2) +
      scale_color_brewer(palette = "Set1", name = "Pool Size (m)") +
      theme_bw(base_size = 16)+ 
      theme(
        axis.text = element_text(size = 14),        # Increase tick label size
        axis.title = element_text(size = 16),       # Axis title size (already covered by base_size)
        legend.text = element_text(size = 12),      # Optional: legend text
        legend.title = element_text(size = 14)      # Optional: legend title
      )
  })
  
  output$pool_count <- renderText({
    n <- input$n
    m <- input$m
    if (m > 0) {
      k <- round(100*n / m)/100
      paste("Number of Kits (Pools):", k)
    } else {
      "Invalid pool size (m must be > 0)"
    }
  })
  
  output$heatmapPlot <- renderPlot({
    p <- input$y_true
    n <- input$n
    m <- input$m
    k <- n/m    # pool count
    
    sens_range <- seq(input$xy_lim[1] / 100, input$xy_lim[2] / 100, length.out = 100)
    spec_range <- seq(input$xy_lim[1] / 100, input$xy_lim[2] / 100, length.out = 100)
    
    # compute_se <- function(p, s, f) {
    #   sqrt((f + (s - f) * p) * (1 - f - (s - f) * p) / (n * (s - f)^2))
    # }
    # 
    # compute_se_pooled <- function(p, s, f, n, m) {
    #   q <- (1 - p)^m                          # prob pool has no carriers
    #   y <- f + (s - f) * (1 - q)              # expected pool positive rate
    #   dy_dp <- (s - f) * m * (1 - p)^(m - 1)  # derivative dy/dp
    #   se_y <- sqrt(y * (1 - y) / n)           # SE from binomial(n, y)
    #   se_p <- se_y / abs(dy_dp)               # delta method SE
    #   return(se_p)
    # }
    
    compute_se_pooled_fixed_total <- function(p, s, f, n, m) {
      k <- n / m                                 # number of pools
      q <- (1 - p)^m                             # prob a pool has no carriers
      y <- f + (s - f) * (1 - q)                 # expected pool positive rate
      dy_dp <- (s - f) * m * (1 - p)^(m - 1)     # derivative dy/dp
      se_y <- sqrt(y * (1 - y) / k)              # SE from binomial(n, y)
      se_p <- se_y / abs(dy_dp)                  # delta method SE
      return(se_p)
    }
    
    # SE <- outer(sens_range, 1 - spec_range, Vectorize(function(s, f) 100 * compute_se(p, s, f)))
    # SE <- outer(sens_range, 1 - spec_range, Vectorize(function(s, f) 100 * compute_se_pooled(p, s, f, n, m)))
    SE <- outer(sens_range, 1 - spec_range, Vectorize(function(s, f) 100 * compute_se_pooled_fixed_total(p, s, f, n, m)))
    
    
    par(pty = "s", mar = c(5, 5, 4, 7))  # Extra margin on the right
    layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1.2))  # Adjust widths as needed
    
    # ---- Heatmap Panel ----
    par(pty = "s", mar = c(5, 5, 4, 2),
        cex.lab = 1.4,           # Axis label font size (xlab, ylab)
        cex.axis = 1.2           # Axis tick label font size (numbers))  # Square plot + space for labels
    )
    
    image(100 * sens_range, 100 * spec_range, SE,
          col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
          xlab = "Sensitivity (%)", 
          ylab = "Specificity (%)",
          asp = 1)
    
    # Automatically compute 6 pretty contour levels in percentage units
    contour_levels <- pretty(100 * seq(min(SE), max(SE), length.out = 6)) # 100 * c(0.04,0.03,0.02,0.01,0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2, 0.25)
    
    contour(100 * sens_range, 100 * spec_range, SE,
            levels = 100 * c(0.04,0.03,0.02,0.01,0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.20, 0.25,0.30,0.40,0.50),
            add = TRUE, col = "black", lty = 2)
    
    # Add red dotted diagonal line: Sensitivity = Specificity
    abline(a = 0, b = 1, col = "#ca0020", lty = 3, lwd = 1)
    
    # ---- Colorbar Panel ----
    par(pty = "m", mar = c(5, 0, 4, 5), bty = "n")  # Reset pty and set clean margins for colorbar
    # par(mar = c(5, 0, 4, 5), bty = "n")  # leave space only on right
    
    # Set up the z vector for colorbar scale
    z_vals <- seq(min(SE), max(SE), length.out = 100)
    
    # Plot the colorbar
    par(bty = "n")  # removes box around the plot area
    image(
      x = 1,  # <- prevent default "1" label from appearing
      y = z_vals,
      z = t(matrix(z_vals, ncol = 1)),
      col = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
      xaxt = "n", yaxt = "n",
      asp = NA,
      xlab = "", ylab = ""
    )
    
    # Axis and label tightly on the right
    axis_ticks <- pretty(range(z_vals), n = 5)
    axis(4, at = axis_ticks, labels = round(axis_ticks, 2), las = 1, cex.axis = 1.1)
    mtext("Standard Error (%)", side = 4, line = 3, cex = 1.2)
    
    
  })
}

shinyApp(ui = ui, server = server)
