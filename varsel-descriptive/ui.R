library(shiny)
library(shinycssloaders)
library(markdown)

# Helper: sidebar section heading with a left border accent
sh <- function(label) tagList(
  tags$h4(label, style = "font-size:15px; font-weight:bold; margin-top:14px;
                           margin-bottom:0; padding-left:8px;
                           border-left:3px solid #aaa;"),
  tags$div(style = "height:10px; line-height:10px; font-size:1px; overflow:hidden;", HTML("&nbsp;"))
)

navbarPage(
  title = "Comparison of methods for variable selection",
  header = tagList(
    withMathJax(),
    tags$head(
      tags$style(HTML("
      .form-group        { margin-bottom: 8px; }
      label              { font-size: 13px; margin-bottom: 2px; }
      .radio label,
      .checkbox label    { font-size: 13px; }
      .well              { padding: 12px; }
      #do, #downloadPlot { width: 100%; }
      h4 + .form-group,
      h4 + .shiny-input-container { margin-top: 10px; }
      .simdesign-content {
        max-width: 900px;
        margin: 0 auto;
        padding: 10px 20px 30px 20px;
        font-size: 14px;
        line-height: 1.6;
      }
      .simdesign-content h1 {
        font-size: 1.8em;
        margin-top: 16px;
        margin-bottom: 16px;
        border-bottom: 2px solid #ccc;
        padding-bottom: 6px;
      }
      .simdesign-content h2 {
        font-size: 1.3em;
        margin-top: 24px;
        border-bottom: 1px solid #ddd;
        padding-bottom: 4px;
      }
      .simdesign-content h3 {
        font-size: 1.05em;
        margin-top: 18px;
        color: #2c3e50;
      }
      .simdesign-content img {
        max-width: 100%;
        display: block;
        margin: 12px auto;
      }
      .simtbl { border-collapse: collapse; font-size: 13px; margin: 10px auto; }
      .simtbl th, .simtbl td { border: 1px solid #bbb; padding: 6px 12px; }
      .simtbl th { background-color: #f0f0f0; text-align: center; }
      .simtbl td:first-child { text-align: left; }
      .simtbl td:not(:first-child) { text-align: center; }
      .simgrid { border-collapse: collapse; font-size: 13px; margin: 10px auto; }
      .simgrid td { border: 1px solid #bbb; padding: 6px 10px; }
      .simgrid td:first-child { text-align: left; }
      .simgrid td:not(:first-child) { text-align: center; }
      .simgrid .rowhdr { vertical-align: middle; font-weight: 500; }
      .simgrid .italic { font-style: italic; }
      .intro-box {
        background-color: #f8f9fa;
        border-left: 4px solid #4a90d9;
        border-radius: 4px;
        padding: 12px 16px;
        margin-bottom: 14px;
        font-size: 14px;
        color: #333;
      }
      .intro-box h4 {
        margin-top: 0;
        margin-bottom: 6px;
        font-size: 17px;
        color: #2c3e50;
      }
      .howto-box {
        background-color: #fff;
        border: 1px solid #dde;
        border-radius: 4px;
        padding: 10px 16px;
        margin-bottom: 14px;
        font-size: 14px;
        color: #333;
      }
      .howto-box h4 {
        margin-top: 0;
        margin-bottom: 6px;
        font-size: 17px;
        color: #2c3e50;
      }
      .howto-box ol {
        margin: 0;
        padding-left: 20px;
      }
      .howto-box ol li {
        margin-bottom: 3px;
      }
    "))
    )
  ),
  
  # ── Tab 1: Simulation results ──────────────────────────────────────────────
  tabPanel(
    "Simulation results",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        # ── Regression type ──────────────────────────────────────────────────
        sh("Regression type"),
        radioButtons("linlog", label = NULL,
                     choices = list("Linear regression"  = 1,
                                    "Logistic regression" = 2)),
        
        # Scenario — linear (main, low R² only)
        conditionalPanel(
          condition = "input.linlog == 1",
          selectInput("scenario", label = "Scenario:",
                      choices = c(
                        "strong signal" = 1,
                        "weak signal"   = 2))
        ),
        
        # Scenario — logistic (main/low R² × rate 0.3/0.05 only)
        conditionalPanel(
          condition = "input.linlog == 2",
          selectInput("scenario", label = "Scenario:",
                      choices = c(
                        "strong signal with event rate 0.3"  = 1,
                        "weak signal with event rate 0.3"    = 2,
                        "strong signal with event rate 0.05" = 3,
                        "weak signal with event rate 0.05"  = 4))
        ),
        
        # Correlation structure
        radioButtons("corr", label = "Correlation structure:",
                     choices = list("Realistic correlation structure" = 1,
                                    "No correlations"                 = 2)),
        
        tags$hr(),
        
        # ── Methods (Uni variants excluded) ───────────────────────────────────
        sh("Methods"),
        checkboxGroupInput(
          "methods", label = NULL,
          choiceNames = list(
            "FU: full model",
            "BE(BIC): backward elimination, BIC",
            HTML("BE(0.05): backward elimination, &alpha; = 0.05"),
            "BE(AIC): backward elimination, AIC",
            HTML("BE(0.50): backward elimination, &alpha; = 0.50"),
            "ABE(AIC): augmented backward elimination, AIC",
            "FSel(AIC): forward selection, AIC",
            "Step_FSel(AIC): stepwise forward selection, AIC",
            "Lasso(CV): Lasso with cross-validation",
            "RLasso(CV): Relaxed Lasso with cross-validation",
            "RLasso(BIC): Relaxed Lasso with BIC",
            "AdaLasso(CV): Adaptive Lasso with cross-validation"
          ),
          choiceValues = list(
            "FU", "BE_BIC", "BE_005", "BE_AIC", "BE_05", "ABE_AIC",
            "FSel_AIC", "Step_FSel_AIC", "Lasso", "RLasso", "RLasso_BIC",
            "AdaLasso"
          ),
          selected = "FU"
        ),
        
        tags$hr(),
        
        # ── Performance measure (4 measures only) ─────────────────────────────
        sh("Performance measure"),
        radioButtons(
          "measures", label = NULL,
          choiceNames = list(
            "Bias of coefficients",
            "Variable selection rates (TPR and FPR)",
            "True / over- / under-selection model rate",
            "Model size (number of selected variables)"
          ),
          choiceValues = as.list(1:4)
        ),
        
        # Options section — only shown when at least one option is relevant
        conditionalPanel(
          condition = "input.measures == 1 | input.measures == 2",
          tags$hr(),
          sh("Options"),
          conditionalPanel(
            condition = "input.measures == 1",
            checkboxInput("cond", "Conditional on selection", value = FALSE)
          ),
          checkboxInput("aggr", "Average over predictors / noise variables", value = FALSE)
        ),
        
        tags$hr(),
        
        actionButton("do", "Update plot", class = "btn-primary")
        
      ),
      
      # ── Main panel ────────────────────────────────────────────────────────
      mainPanel(
        width = 9,
        
        # ── Introduction ──────────────────────────────────────────────────────
        tags$div(
          class = "intro-box",
          tags$h4("About this app"),
          tags$p(
            "This app presents results from a simulation study comparing variable selection methods
             for multivariable linear and logistic regression in a descriptive modelling context.
             The simulation comprised ",
            tags$strong("20 independent variables (10 true predictors and 10 noise variables)"),
            " whose marginal distributions were derived from data of the National Health and Nutrition
             Examination Survey (NHANES). Two settings for the correlation structure are available:
             a realistic correlation structure derived from NHANES, and a setting with no correlations
             between variables (both sharing the same marginal distributions). Predictor effects were
             chosen to represent a realistic mixture of stronger and weaker associations. Three
             regression types were evaluated: ",
            tags$strong("linear regression, and logistic regression with event rates of 0.3 and 0.05,"),
            " across a range of sample sizes. Two signal strength settings were considered:"
          ),
          tags$table(
            class = "table table-sm table-bordered",
            style = "width:auto; font-size:13px; margin-bottom:10px;",
            tags$thead(tags$tr(
              tags$th("Regression type"),
              tags$th("Strong signal"),
              tags$th("Weak signal")
            )),
            tags$tbody(
              tags$tr(
                tags$td("Linear regression"),
                tags$td(HTML("R&sup2; = 0.45")),
                tags$td(HTML("R&sup2; = 0.15"))
              ),
              tags$tr(
                tags$td("Logistic regression, event rate 0.3"),
                tags$td(HTML("Cox-Snell R&sup2; = 0.40")),
                tags$td(HTML("Cox-Snell R&sup2; = 0.13"))
              ),
              tags$tr(
                tags$td("Logistic regression, event rate 0.05"),
                tags$td(HTML("Cox-Snell R&sup2; = 0.16")),
                tags$td(HTML("Cox-Snell R&sup2; = 0.05"))
              )
            )
          ),
          tags$p(
            "For full details on the simulation design, see the ",
            tags$strong("Simulation Design"), " tab."
          )
        ),
        tags$div(
          class = "howto-box",
          tags$h4("How to use"),
          tags$ol(
            tags$li(HTML(
              "<strong>Choose a setting:</strong> Select a regression type, scenario (strong or
               weak signal), and correlation structure on the left."
            )),
            tags$li(HTML(
              "<strong>Select methods:</strong> Tick one or more variable selection methods to compare."
            )),
            tags$li(HTML(
              "<strong>Choose a performance measure:</strong> Select what aspect of the methods
               you want to evaluate &mdash; bias, variable selection rates, model selection rates,
               or model size."
            )),
            tags$li(HTML(
              "<strong>Update the plot:</strong> Click <strong>Update plot</strong> to display results."
            ))
          )
        ),
        uiOutput("mainLegend"),
        uiOutput("methodsLegend"),
        textOutput("colors"),
        textOutput("biasNote"),
        tags$hr(),
        withSpinner(plotOutput("measPlot"), type = 4)
      )
    )
  ),
  
  # ── Tab 2: Simulation design ───────────────────────────────────────────────
  tabPanel(
    "Simulation design",
    tags$div(class = "simdesign-content",
             withMathJax(includeMarkdown("www/design.md")))
  ),
  
  # ── Tab 3: About ──────────────────────────────────────────────────────────
  tabPanel(
    "About",
    fluidRow(
      column(10, offset = 1,
             tags$br(),
             includeHTML("www/about.html"))
    )
  )
)