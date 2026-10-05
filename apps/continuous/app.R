# Launch from repository root: Rscript scripts/run_app.R continuous
library(shiny)
source("R/regression_power.R", local = TRUE)
source("R/training_data.R", local = TRUE)
source("R/regression_grant_text.R", local = TRUE)
training <- load_applet_training()

ui <- fluidPage(
  titlePanel("Power prototype: continuous outcome"),
  tags$head(tags$style(HTML("#estimates table {width:100%;table-layout:fixed;font-size:13px;} #estimates th,#estimates td {white-space:normal;overflow-wrap:anywhere;}"))),
  p("Continuous-outcome association • independent participants or complete participant pairs"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("mode", "Analysis", c("Cross-sectional: one visit" = "cross", "Longitudinal: repeated-measures association" = "mixed"), selected = "cross"),
      selectInput("dataset", "Training feature SDs", unique(training$key),
        selected = "Original: HMP1-2 taxonomy stool"),
      numericInput("n", "Independent participants / complete pairs", 200, min = 4, step = 1),
      numericInput("covariates", "Adjustment coefficients (excluding outcome, intercept and time)", 2, min = 0, step = 1),
      helpText("Count K−1 coefficients for a K-level categorical covariate. Longitudinal mode adds a visit/time term automatically."),
      sliderInput("target", "Target power", min = 0.5, max = 0.99, value = 0.8, step = 0.01),
      numericInput("family_alpha", "Family-wise significance level", 0.05, min = 0.0001, max = 0.2, step = 0.01),
      numericInput("features", "Number of features tested", 260, min = 1, step = 1),
      numericInput("analyses", "Tests per feature in the correction family", 1, min = 1, step = 1),
      helpText("Use 2 if correcting jointly across separate tests at two visits. This changes alpha only; it does not pool visits."),
      radioButtons("outcome_r2", "Assumed outcome variation explained by covariates",
        choices = c("0%" = "0", "5%" = "0.05", "10%" = "0.10", "15%" = "0.15", "20%" = "0.20"), selected = "0.05", inline = TRUE),
      textInput("outcome_name", "Outcome name for the summary", "outcome"),
      textInput("covariate_names", "Adjustment covariates for the summary", "age and sex"),
      helpText("These names label the paragraph. Set the matching number of covariate coefficients above."),
      helpText("Abundance examples are standardized to a starting relative abundance of 1%."),
      helpText("Use R² from clinical outcome ~ all adjustment covariates, without microbial features. Do not substitute a microbiome PERMANOVA R²."),
      helpText("Unknown? Compare the scenarios below. Zero assumes no overlap with covariates and is optimistic for this part of the calculation."),
      conditionalPanel("input.mode == 'mixed'",
        sliderInput("feature_rho", "Similarity of feature abundance between visits (ICC)", min = 0, max = 0.95, value = 0.5, step = 0.05),
        helpText("0 means no residual similarity; values near 1 mean strong similarity after adjustment. Default 0.5 is an assumption, not an estimate."),
        sliderInput("outcome_rho", "Similarity of outcome between visits after adjustment", min = -0.9, max = 1, value = 0.7, step = 0.05),
        helpText("This is a separate assumption about outcome, not feature abundance. Default 0.7 is illustrative. Both visits must measure outcome.")),
      numericInput("sd_multiplier", "Per-visit unexplained SD / training feature SD", 1, min = 0.01, step = 0.1),
      helpText("Default 1 uses the training SD as unexplained per-visit SD. In longitudinal mode this includes both the subject random intercept and visit-level error."),
      numericInput("beta", "Slope magnitude to evaluate (transformed units per predictor SD)", 0.005, min = 0, step = 0.001),
      downloadButton("download", "Download estimates")
    ),
    mainPanel(
      h3("Grant summary — select and copy"),
      uiOutput("grantSummary"),
      downloadButton("downloadSummary", "Download summary text"),
      p("The paragraph updates with the settings. outcome is expressed in standard deviations on the eventual analysis scale; you do not need to choose raw versus log outcome to explore these standardized scenarios. The abundance example is a back-transformed fitted value, not a universal raw-scale effect. In longitudinal mode, the example describes a difference in fitted abundance associated with higher outcome after accounting for repeated measurements, not a before-and-after change within each person."),
      h3("Model and assumptions"),
      uiOutput("model"),
      uiOutput("design"),
      p(strong("SD planning assumption: "), "training SD quartiles describe positive abundances only. They are not measured full-cohort regression residual SDs. These results are conditional on the assumed residual SD; they are not a validated zero-inclusive power estimate. Participant count must match the samples actually included in each fitted model."),
      h3("Detectable slopes and power"),
      tableOutput("estimates"),
      p("Slope is measured per one SD of the outcome in both modes; longitudinal outcome is standardized across the two visits, not as a change score. Each row uses its own SD quartile; 75% is the most variable of the three. No baseline-independent raw-abundance conversion is applied."),
      plotOutput("curves", height = "400px"),
      h3("If you do not know the covariate R²"),
      p("The following are illustrative assumptions, not estimates or bounds. Each row uses the 75th-percentile feature SD, with all other settings held fixed. Higher R² means less independent outcome variation remains, so a larger slope is needed."),
      tableOutput("r2Sensitivity"),
      tags$details(tags$summary("How to estimate the correct R² from your metadata"),
        p("For a cross-sectional analysis, use the same complete participants and covariate coding as the planned feature model. Fit the clinical outcome on all adjustment covariates together, without microbial features. Use ordinary Multiple R-squared, not Adjusted R-squared or a sum of separate covariate R² values."),
        tags$pre("summary(lm(clinical_outcome ~ age + sex + BMI, data = analysis_data))$r.squared"),
        p("Replace age, sex and BMI with your actual adjustment terms. In longitudinal mode, R² summarizes outcome variation explained by all fixed adjustment terms, including time. The separate outcome correlation input describes what remains after that adjustment."),
        p("A PERMANOVA on microbial community distances measures a different response. It also cannot supply the residual SD of each individual feature. The residual-SD multiplier is a separate planning assumption.")),
      p("Longitudinal mode uses a two-visit random-intercept model and a large-sample Wald power approximation. It assumes a common outcome coefficient across within-person and between-person associations; it does not isolate within-person change."),
      tags$details(tags$summary("Calculation details"),
        p("Cross-sectional: residual df = N − covariate coefficients − 2; residualized predictor sum of squares = (N − 1)(1 − outcome R²)."),
        p("Longitudinal: Y_it = intercept + beta*outcome_it + time + covariates + b_i + e_it. Marginal unexplained SD S = training SD × multiplier; Var(b) = ICC × S² and Var(e) = (1 − ICC) × S²."),
        p("For two balanced visits, coefficient information = [(2N − 1)/2](1 − R²)[(1 + outcome correlation)/(1 + feature ICC) + (1 − outcome correlation)/(1 − feature ICC)] / S². This is the GLS information after projecting out fixed nuisance terms under the assumed balanced predictor design."),
        p("Mixed-model power uses a two-sided normal/Wald test, treating variance components as known planning values. It is a large-sample approximation, not an exact small-sample t test. Covariate count does not add a separate finite-df penalty in this mode; adjustment affects information through R² and the assumed residual predictor correlation."),
        p("For slope β, noncentrality λ = β × sqrt(residualized predictor sum of squares) / residual SD. Power is the upper-tail probability of noncentral F(1, df, λ²) above the two-sided coefficient-test threshold."),
        p("Per-test alpha = family-wise alpha / tested features / tests per feature. Adjustment covariates consume degrees of freedom; they are not automatically counted as additional tested hypotheses."),
        p("Assumes independent participants, normal random intercepts independent of predictors, independent equal-variance visit-level errors, two complete visits and a full-rank fixed-effects design. Variance components and predictor moments are assumed. No random slopes, outcome-by-time interaction, or missing visits are modeled."),
        tags$a(href = "https://stat.ethz.ch/R-manual/R-devel/library/stats/html/TDist.html", "R noncentral t distribution documentation"))
    )
  )
)
server <- function(input, output, session) {
  selected <- reactive({ req(input$dataset); training[training$key == input$dataset, , drop = FALSE] })
  observeEvent(input$dataset, {
    updateNumericInput(session, "features", value = selected()$total_features[1])
  }, ignoreInit = FALSE)
  settings <- reactive({
    req(input$n, !is.null(input$covariates), input$features, input$analyses,
        input$family_alpha, input$target, !is.null(input$outcome_r2), input$sd_multiplier, !is.null(input$beta))
    r2 <- as.numeric(input$outcome_r2)
    whole <- function(x) is.finite(x) && x == floor(x)
    validate(need(whole(input$n) && input$n > input$covariates + 2, "Participants must exceed covariate coefficients + 2."),
      need(whole(input$covariates) && input$covariates >= 0, "Enter a nonnegative integer covariate count."),
      need(whole(input$features) && input$features >= 1, "Enter a positive integer feature count."),
      need(whole(input$analyses) && input$analyses >= 1, "Enter a positive integer test count."),
      need(is.finite(input$family_alpha) && input$family_alpha > 0 && input$family_alpha < 1, "Significance must be between 0 and 1."),
      need(is.finite(input$sd_multiplier) && input$sd_multiplier > 0, "Residual-SD multiplier must be positive."),
      need(is.finite(r2) && r2 %in% c(0,.05,.10,.15,.20), "Outcome R² must be in [0, 1)."),
      need(input$mode %in% c("cross", "mixed"), "Choose an analysis mode."),
      need(input$mode != "mixed" || (is.finite(input$feature_rho) && input$feature_rho >= 0 && input$feature_rho < 1), "Feature ICC must be between 0 and 1, excluding 1."),
      need(input$mode != "mixed" || (is.finite(input$outcome_rho) && abs(input$outcome_rho) <= 1), "outcome correlation must be between -1 and 1."),
      need(input$mode != "mixed" || input$n > input$covariates + 3, "Participants must exceed adjustment coefficients + 3."),
      need(is.finite(input$beta) && input$beta >= 0, "Slope magnitude must be nonnegative."))
    alpha <- input$family_alpha / input$features / input$analyses
    validate(need(input$target > alpha && input$target < 1, "Target power must exceed per-test alpha and be below 1."))
    list(n = input$n, covariates = input$covariates, alpha = alpha,
      outcome_r2 = r2, sd_multiplier = input$sd_multiplier, mode = input$mode,
      feature_rho = if (input$mode == "mixed") input$feature_rho else .5,
      outcome_rho = if (input$mode == "mixed") input$outcome_rho else .7)
  })
  estimates <- reactive({
    t <- selected(); s <- settings()
    do.call(rbind, lapply(seq_len(nrow(t)), function(i) {
      sd <- as.numeric(t[i, c("sd_25", "sd_50", "sd_75")])
      detectable <- do.call(analysis_detectable, c(list(feature_sd = sd, target = input$target), s))
      achieved <- do.call(analysis_power, c(list(beta = input$beta, feature_sd = sd), s))
      data.frame(Class = t$feature_class[i], Features = t$class_features[i],
        Quartile = c("25%", "50%", "75%"), Training_SD = sd,
        Assumed_residual_SD = sd * s$sd_multiplier, Detectable_slope = detectable,
        Power_at_entered_slope = achieved)
    }))
  })
  grant_text <- reactive({
    s <- settings(); t <- selected()
    req(input$outcome_name, input$covariate_names)
    make_regression_grant_text(t, s, input$target, input$features, input$analyses,
      input$family_alpha, input$outcome_name, input$covariate_names,
      1, input$mode, input$feature_rho)
  })
  output$grantSummary <- renderUI({
    tags$textarea(readonly = "readonly", style = "width:100%;min-height:280px;font-size:16px;line-height:1.5;padding:12px;",
      `aria-label` = "Grant summary", grant_text())
  })
  output$downloadSummary <- downloadHandler(filename = function() "regression_grant_summary.txt",
    content = function(file) writeLines(grant_text(), file, useBytes = TRUE))
  output$model <- renderUI({
    if (input$mode == "mixed") tagList(
      p(code(paste0("asin(sqrt(feature relative abundance)) ~ standardized ", input$outcome_name,
        " + time + ", input$covariate_names, " + (1 | subjectID)"))),
      p("Tests the overall association with the outcome while accounting for repeated measurements. It does not test feature change versus outcome change."),
      p("Uses two complete visits per participant, a fixed visit/time effect and a subject random intercept. Standardize the outcome across visits before fixed-effect adjustment. Both correlation settings are assumptions; compare plausible values.")) else tagList(
      p(code("asin(sqrt(feature relative abundance)) ~ standardized outcome + covariates")),
      p("The continuous outcome is the predictor here, following feature ~ outcome. Standardize it to SD 1 before covariate adjustment."))
  })
  output$design <- renderUI({
    t <- selected(); s <- settings()
    tagList(p(t$definition[1]),
      p("Training samples: ", t$samples[1], "; tested feature family in this training set: ", t$total_features[1],
        "; excluded input features: ", t$excluded[1], "."),
      if (input$mode == "mixed") p(strong(s$n), " independent participants, ", 2*s$n,
        " observations across two visits. Random-intercept model; large-sample power approximation.") else
        p(strong(s$n), " independent participants. Residual degrees of freedom: ", s$n - s$covariates - 2, "."),
      p("Per-test alpha: ", input$family_alpha, " / ", input$features, " / ", input$analyses,
        " = ", format(s$alpha, scientific = TRUE), ". Target power: ", input$target, "."))
  })
  output$estimates <- renderTable({
    display <- estimates()
    names(display) <- c("Class", "Features", "SD quartile", "Training SD", "Residual SD", "Detectable slope", "Power at slope")
    display
  }, digits = 6, striped = TRUE)
  r2_sensitivity <- reactive({
    t <- selected(); s <- settings()
    do.call(rbind, lapply(c(0, .05, .10, .15, .20), function(r2) {
      scenario <- s; scenario$outcome_r2 <- r2
      slope <- do.call(analysis_detectable,
        c(list(feature_sd = t$sd_75, target = input$target), scenario))
      data.frame(Class = t$feature_class, Assumed_R2 = r2,
        Outcome_variation_remaining = 1-r2, Detectable_slope_75_SD = slope)
    }))
  })
  output$r2Sensitivity <- renderTable({
    display <- r2_sensitivity()
    names(display) <- c("Class", "Assumed R²", "Outcome variance remaining", "Detectable slope (75% SD)")
    display
  }, digits = 6, striped = TRUE)
  output$curves <- renderPlot({
    e <- estimates(); s <- settings(); classes <- unique(e$Class)
    par(mfrow = c(1, length(classes)), mar = c(5, 4, 3, 1))
    for (class_name in classes) {
      rows <- e[e$Class == class_name, ]
      x <- seq(0, max(rows$Detectable_slope) * 1.3, length.out = 180)
      y <- sapply(rows$Training_SD, function(sd)
        do.call(analysis_power, c(list(beta = x, feature_sd = sd), s)))
      matplot(x, y, type = "l", col = c("#6689A6", "#176B83", "#A65326"),
        lty = c(2, 1, 3), lwd = 2, ylim = c(0, 1), main = class_name,
        xlab = "Slope per outcome SD", ylab = "Power")
      abline(h = input$target, col = "gray50", lty = 2)
      legend("bottomright", c("25% SD", "50% SD", "75% SD"),
        col = c("#6689A6", "#176B83", "#A65326"), lty = c(2, 1, 3), bty = "n")
    }
  })
  output$download <- downloadHandler(filename = function() "regression_power_estimates.csv",
    content = function(file) {
      s <- settings()
      data <- estimates()
      data$analysis <- input$mode
      data$feature_residual_correlation <- if (input$mode == "mixed") input$feature_rho else NA_real_
      data$outcome_residual_correlation <- if (input$mode == "mixed") input$outcome_rho else NA_real_
      data$visits_per_participant <- if (input$mode == "mixed") 2 else 1
      data$method <- if (input$mode == "mixed") "Two-visit random-intercept GLS/Wald approximation" else "Fixed-design noncentral t"
      data$dataset <- input$dataset
      data$participants <- s$n
      data$covariate_coefficients <- s$covariates
      data$outcome_covariate_R2 <- s$outcome_r2
      data$reference_abundance_percent <- 1
      data$per_visit_sd_multiplier <- input$sd_multiplier
      data$effective_sd_multiplier <- s$sd_multiplier
      data$family_alpha <- input$family_alpha
      data$tested_features <- input$features
      data$tests_per_feature <- input$analyses
      data$per_test_alpha <- s$alpha
      data$target_power <- input$target
      data$entered_slope <- input$beta
      data$sd_basis <- "Positive-abundance training SD used as assumed regression residual SD"
      write.csv(data, file, row.names = FALSE)
    })
}
shinyApp(ui, server)
