library(shiny)
library(plotly)
library(ggrepel)

# Calculations for simulation results
shinyServer(function(input, output, session){
  
  methods_selected <- reactive({input$methods})
  cond_selected    <- reactive({input$cond})
  aggr_selected    <- reactive({input$aggr})
  
  rescale_const <- 1
  
  downloadInput <- plotInput <- pngInput <- textInput <- plotInputHeight <- reactiveValues()
  
  observeEvent(input$do, {
    
    withProgress(message = "Building plot\u2026", value = 0, {
      
      # ── Step 1: load data ──────────────────────────────────────────────────
      incProgress(0.15, detail = "Loading scenario data\u2026")
      
      if(input$linlog == 1) {
        logregr = FALSE
        
        if (input$corr == 1) {
          # Realistic correlation structure
          if (input$scenario == 1){
            load("results/linear regression/results_sim_linear_main_preparation_for_shiny.RData")
            results <- results_sim_linear_main}
          if (input$scenario == 2){
            load("results/linear regression/results_sim_linear_lowR2_preparation_for_shiny.RData")
            results <- results_sim_linear_lowR2}
        } else {
          # No correlations
          if (input$scenario == 1){
            load("results/linear regression/results_sim_smpl3_linear_main_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_main}
          if (input$scenario == 2){
            load("results/linear regression/results_sim_smpl3_linear_lowR2_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_lowR2}
        }
      }
      
      if(input$linlog == 2) {
        logregr = TRUE
        
        if (input$corr == 1) {
          # Realistic correlation structure
          if (input$scenario == 1){
            load("results/logistic regression/results_sim_linear_rate03_main_preparation_for_shiny.RData")
            results <- results_sim_linear_rate03_main}
          if (input$scenario == 2){
            load("results/logistic regression/results_sim_linear_rate03_lowR2_preparation_for_shiny.RData")
            results <- results_sim_linear_rate03_lowR2}
          if (input$scenario == 3){
            load("results/logistic regression/results_sim_linear_rate005_main_preparation_for_shiny.RData")
            results <- results_sim_linear_rate005_main}
          if (input$scenario == 4){
            load("results/logistic regression/results_sim_linear_rate005_lowR2_preparation_for_shiny.RData")
            results <- results_sim_linear_rate005_lowR2}
        } else {
          # No correlations
          if (input$scenario == 1){
            load("results/logistic regression/results_sim_smpl3_linear_rate03_main_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_rate03_main}
          if (input$scenario == 2){
            load("results/logistic regression/results_sim_smpl3_linear_rate03_lowR2_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_rate03_lowR2}
          if (input$scenario == 3){
            load("results/logistic regression/results_sim_smpl3_linear_rate005_main_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_rate005_main}
          if (input$scenario == 4){
            load("results/logistic regression/results_sim_smpl3_linear_rate005_lowR2_preparation_for_shiny.RData")
            results <- results_sim_smpl3_linear_rate005_lowR2}
        }
      }
      
      # ── Step 2: compute performance measure & build plot ──────────────────
      incProgress(0.45, detail = "Computing performance measures\u2026")
      
      methods = methods_selected()
      
      # Measure values: 1=bias, 2=TPR/FPR, 3=true/over/under model, 4=model size
      if(input$measures == 1){
        plotInput$plot <- fun_plot_bias(results, methods, cond_selected(), aggr_selected(),
                                        estimand = 1, rescale_const)
      }
      if(input$measures == 2){
        plotInput$plot <- fun_plot_TPR_FPR(results, methods, aggr_selected())
      }
      if(input$measures == 3){
        plotInput$plot <- fun_plot_truemodel(results, methods)
      }
      if(input$measures == 4){
        plotInput$plot <- fun_plot_model_size(results, methods)
      }
      
      # ── Step 3: prepare labels & metadata ─────────────────────────────────
      incProgress(0.30, detail = "Preparing labels\u2026")
      
      pngInput$title <- c("Bias", "TPR_FPR", "Models", "Model_size")[as.numeric(input$measures)]
      
      allScenariosLin <- c(
        "Strong signal scenario",
        "Weak signal scenario")
      
      allScenariosLog <- c(
        "Strong signal scenario with event rate 0.3",
        "Weak signal scenario with event rate 0.3",
        "Strong signal scenario with event rate 0.05",
        "Weak signal scenario with event rate 0.05")
      
      allMeasures <- c(
        "Bias of estimated standardized regression coefficients",
        "Variable selection rates of predictors (true positive rate, TPR) and noise variables (false positive rate, FPR)",
        "Model selection rates of true, over-selection and under-selection models. True model: exactly the true predictors were selected. Over-selection: the true predictors and at least one noise variable were selected. Under-selection: at least one true predictor was not selected",
        "Model size (nr of selected variables)")
      
      allMethods <- c(
        "FU"            = "FU, Full model",
        "BE_BIC"        = "BE(BIC), Backward elimination with BIC",
        "BE_005"        = "BE(0.05), Backward elimination with alpha = 0.05",
        "BE_AIC"        = "BE(AIC), Backward elimination with AIC",
        "BE_05"         = "BE(0.50), Backward elimination with alpha = 0.50",
        "ABE_AIC"       = "ABE(AIC), Augmented backward elimination with AIC and tau=0.05",
        "FSel_AIC"      = "FSel(AIC), Forward selection with AIC",
        "Step_FSel_AIC" = "Step_FSel(AIC), Stepwise forward selection with AIC",
        "Lasso"         = "Lasso(CV), Least angle selection and shrinkage operator, lambda tuned with cross-validation",
        "RLasso"        = "RLasso(CV), Relaxed Lasso - OLS fit with variables selected by Lasso, lambda tuned with cross-validation",
        "RLasso_BIC"    = "RLasso(BIC), Relaxed Lasso - OLS fit with variables selected by Lasso, lambda tuned with BIC",
        "AdaLasso"      = "AdaLasso(CV), Adaptive Lasso, lambda tuned with cross-validation")
      
      corrLabel <- ifelse(input$corr == 1, "realistic correlation structure", "no correlations")
      
      textInput$LegendScenario <- ""
      if (input$linlog == 1) {
        textInput$LegendScenario <- paste0(allScenariosLin[match(input$scenario, 1:2)],
                                           ", ", corrLabel)
      } else if (input$linlog == 2) {
        textInput$LegendScenario <- paste0(allScenariosLog[match(input$scenario, 1:4)],
                                           ", ", corrLabel)
      }
      
      textInput$Measure <- allMeasures[as.numeric(input$measures)]
      if(cond_selected() & input$measures == 1)
        textInput$Measure <- paste0(textInput$Measure, ", conditional on selection")
      textInput$Methods <- allMethods[match(methods_selected(), names(allMethods))]
      
      textInput$colors <- ""
      if(input$measures %in% 1:2) {
        if(aggr_selected() == FALSE)
          textInput$colors <- "Predictors are represented by solid lines, and noise variables by dashed lines. The stronger the effect of a predictor, the thicker the line. The color of the line indicates the multiple R\u00b2 of a predictor or noise variable."
        else
          textInput$colors <- "Results are averaged over predictors (solid line) and over noise variables (dashed line)."
      }
      
      textInput$linebreak <- "-----"
      
      textInput$biasNote1 <- ""
      if(input$measures == 1 & aggr_selected() == FALSE)
        textInput$biasNote1 <- "For predictors, bias>0 denotes bias away from 0 and bias<0 denotes bias towards 0."
      if(input$measures == 1 & aggr_selected() == TRUE)
        textInput$biasNote1 <- "Aggregated bias is computed as mean (over all predictors or noise variables) of absolute values of bias for each individual predictor."
      
      # Plot height
      plotInputHeight$plot_height <- 300
      if(input$measures %in% 1:2) {
        if(length(methods_selected()) >= 4)  plotInputHeight$plot_height <- 500
        if(length(methods_selected()) >= 7)  plotInputHeight$plot_height <- 700
        if(length(methods_selected()) >= 10) plotInputHeight$plot_height <- 900
      }
      if(input$measures %in% 3:4) {
        if(length(methods_selected()) >= 4)  plotInputHeight$plot_height <- 450
        if(length(methods_selected()) >= 7)  plotInputHeight$plot_height <- 600
        if(length(methods_selected()) >= 10) plotInputHeight$plot_height <- 750
      }
      
      # ── Step 4: done ──────────────────────────────────────────────────────
      incProgress(0.10, detail = "Rendering\u2026")
      
    }) # end withProgress
  })
  
  
  
  output$measPlot <- renderPlot({
    plotInput$plot
  }, height = reactive({ ifelse(!is.null(plotInputHeight$plot_height),
                                plotInputHeight$plot_height, 400) }))
  
  output$mainLegend <- renderUI({
    req(textInput$LegendScenario)
    tags$p(paste0(textInput$LegendScenario, ": ", textInput$Measure, "."))
  })
  output$methodsLegend <- renderUI({
    req(textInput$Methods)
    tags$p(paste0(paste(textInput$Methods, collapse = "; "), "."))
  })
  output$linebreak <- renderText({ textInput$linebreak })
  output$biasNote  <- renderText({ textInput$biasNote1  })
  output$colors    <- renderText({ textInput$colors     })
  
  
  
  
})