library(shiny)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(bslib)
library(shinyjs)
library(thematic)
library(markdown)

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
          numericInput("n", HTML("Sample Size, <i>n</i>"), value = 100, step = 1),
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
    updateSliderInput(session, "xy_lim", value = c(60, 100))
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

  output$heatmapPlot <- renderPlot({
    p <- input$y_true
    n <- input$n
    sens_range <- seq(input$xy_lim[1] / 100, input$xy_lim[2] / 100, length.out = 100)
    spec_range <- seq(input$xy_lim[1] / 100, input$xy_lim[2] / 100, length.out = 100)

    compute_se <- function(p, s, f) {
      sqrt((f + (s - f) * p) * (1 - f - (s - f) * p) / (n * (s - f)^2))
    }

    SE <- outer(sens_range, 1 - spec_range, Vectorize(function(s, f) 100 * compute_se(p, s, f)))
    
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
