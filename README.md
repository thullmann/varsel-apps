---
title: "README"
output: html_document
---

# varsel-apps

Two interactive Shiny apps accompanying the paper:

> Ullmann, T., Heinze, G., Kammer, M., Dunkler, D. for TG2 of the STRATOS
> initiative (2026). A preregistered simulation study provided evidence on 
> the appropriate use of data-driven variable selection in descriptive 
> multivariable regression modeling. *Currently under review.*

## Apps

| App | Description | Link |
|-----|-------------|------|
| Simulation results | Explore performance measures (bias, selection rates, model size, …) across methods, scenarios, and sample sizes | https://thullmann.shinyapps.io/varsel_descriptive |
| TPR/FPR trade-off | Visualize whether methods meet user-specified TPR and FPR thresholds | https://thullmann.shinyapps.io/varsel_TPR_FPR |

## Running locally

The simulation result files (`results/`) are not included in this repository
due to their size. They are bundled with the deployed apps. To obtain the results files,
please contact the authors.

```r
install.packages(c("shiny", "shinycssloaders", "ggplot2", "ggrepel",
                   "reshape2", "cowplot", "ggpubr", "RColorBrewer",
                   "pROC", "plotly"))

shiny::runApp("varsel-descriptive")
shiny::runApp("varsel-TPR-FPR")
```
