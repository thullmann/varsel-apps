library(pROC)

##### ##### ##### ##### ##### ##### ##### #####
##### Functions for calculating measures #####
##### ##### ##### ##### ##### ##### ##### #####

## For predictors, fun_bias should return negative values for 'bias towards 0' and positive values for 'bias away from 0'

fun_bias <- # on scale of standardized coefficients
  function(beta_approx, coefficients,  sdev, cond = F, rescale_const = 1){
    
    if(cond == T){
      TF <- bias <- NULL
      for(k in 1:length(beta_approx)){
        TF <- coefficients[,k] != 0
        bias[k] <- mean((coefficients[TF,k] - beta_approx[k]) * sdev[k])
      }
    }else{
      bias <- apply((coefficients -
                       matrix(beta_approx, 
                              ncol = ncol(coefficients), 
                              nrow = nrow(coefficients),
                              byrow = T) ) *
                      matrix(sdev, 
                             ncol = ncol(coefficients), 
                             nrow = nrow(coefficients),
                             byrow = T)
                    , 
                    2, function(x) mean(x))
    }
    ## start add lines to return bias away from (+) or towards (-) zero ##
    sign_beta_approx <- sign(beta_approx)
    sign_beta_approx[sign_beta_approx==0] <- 1
    bias <- bias * sign_beta_approx
    bias <- bias / rescale_const
    ## end add ##
    return(bias)
  }

fun_rmsen <- # on scale of standardized coefficients
  function(beta_approx, coefficients, sdev, n,  cond = F){
    if(cond == T){
      TF <- rmse <- NULL
      for(k in 1:length(beta_approx)){
        TF <- coefficients[,k] != 0
        rmse[k] <-  sqrt(mean(((beta_approx[k] - coefficients[TF,k])*sdev[k])^2)) * sqrt(n)
      }
    }else{
      rmse <- apply((coefficients -
                       matrix(beta_approx, 
                              ncol = ncol(coefficients), 
                              nrow = nrow(coefficients),
                              byrow = T)) *
                      matrix(sdev, 
                             ncol = ncol(coefficients), 
                             nrow = nrow(coefficients),
                             byrow = T), 
                    2, function(x) sqrt(mean(x^2))) * sqrt(n)
    }
    return(rmse)
    
  }


fun_cover <- 
  function(beta_approx, CI.lower, CI.upper, coefficients = NULL, cond = F){
    if(cond == T){
      TF <- cover <- NULL
      for(k in 1:length(beta_approx)){
        TF <- coefficients[,k] != 0
        cover[k] <- mean((beta_approx[k] >= CI.lower[TF,k]) & (beta_approx[k] <= CI.upper[TF,k]))
      }
    }else{
      beta_mat <- matrix(beta_approx, 
                         ncol = ncol(coefficients), 
                         nrow = nrow(coefficients),
                         byrow = T)
      cover <- apply((beta_mat >= CI.lower) & (beta_mat <= CI.upper),
                     2,
                     mean)
    }
    return(cover)
  }


fun_width <- # on scale of standardized coefficients
  function(CI.lower, CI.upper, sdev, n, coefficients = NULL, cond = F){
    if(cond == T){
      TF <- width <- NULL
      for(k in 1:ncol(CI.lower)){
        TF <- coefficients[,k] != 0
        width[k] <- mean((CI.upper[TF,k] - CI.lower[TF,k])* sdev[k]) * sqrt(n)
      }
    }else{
      width <- apply((CI.upper - CI.lower) * matrix(sdev, 
                                                    ncol = ncol(CI.upper), 
                                                    nrow = nrow(CI.upper),
                                                    byrow = T),
                     2, mean) * sqrt(n)
    }
    return(width)
  }



fun_power<- 
  function(CI.lower, CI.upper, coefficients = NULL, cond = F){
    if(cond == T){
      TF <- power <- NULL
      for(k in 1:ncol(CI.lower)){
        TF <- coefficients[,k] != 0
        power[k] <- mean(( 0 < CI.lower[TF,k]) | (0 > CI.upper[TF,k]))
      }
    }else{
      power <- apply((0 < CI.lower) | (0 > CI.upper),
                     2,
                     mean)
    }
    return(power)
  }


fun_TPR_FPR <- 
  function(coefficients){
    PR <-  apply(coefficients, 2, function(x) mean(x != 0))
    return(PR)
  }

fun_model_size <- 
  function(coefficients){
    model_size <-  apply(coefficients, 1, function(x) sum(x != 0))
    model_size <- mean(model_size) - 1 # subtract 1 for intercept 
    return(model_size)
  }


fun_rank <- # on scale of standardized coefficients
  function(beta_approx, coefficients, sdev){
    
    coef_stand <- coefficients  * matrix(sdev, 
                                         ncol = ncol(coefficients), 
                                         nrow = nrow(coefficients),
                                         byrow = T)
    beta_approx_stand = beta_approx * sdev
    
    kendalls_tau <- NULL
    for (j in 1:nrow(coefficients)) {
      kendalls_tau[j] <- cor(abs(coef_stand[j,-1]), abs(beta_approx_stand[-1]), method = "kendall")
    }
    kendalls_tau = mean(kendalls_tau) # mean over simulation repetitions
    
    return(kendalls_tau)
  }




fun_truemodel <- 
  function(coefficients, predictors){
    tm <-  mean(apply(coefficients, 1, function(x) identical(x != 0, predictors))) # true model
    om <-  mean(apply(coefficients, 1, function(x) !(0%in%(x[predictors] != 0)) & (sum(x[!predictors] != 0) > 0) )) # over-selection
    um <- mean(apply(coefficients, 1, function(x) 0%in%(x[predictors] != 0))) # under-selection 
    return(c(tm, om, um))
  }


fun_y_bias <- # on scale of linear predictor. For logistic regression: predictions should be predictions_link 
  function(truelp, predictions){
    localb <- apply(apply(predictions, 1, function(x) (x - truelp)),1, mean)
    globalb <- mean(localb)
    return(list("local_bias" = localb, "global_bias" = globalb))
  }


fun_y_rmse <- # on scale of linear predictor. For logistic regression: predictions should be predictions_link  
  function(truelp, predictions, n){
    localr <- apply(apply(predictions, 1, function(x) (x - truelp)^2), 1, function(y) sqrt(mean(y)) * sqrt(n)) 
    mean_over_test_obs <- apply(apply(predictions, 1, function(x) (x - truelp)^2), 2, function(y) mean(y)) # vector with nsim elements 
    globalr <- sqrt(mean(mean_over_test_obs)) * sqrt(n)
    return(list("local_rmse" = localr, "global_rmse" = globalr))
  }

fun_mae <- # on scale of linear predictor. For logistic regression: predictions should be predictions_link  
  function(truelp, predictions) {
    mean_over_test_obs <- apply(apply(predictions, 1, function(x) abs(x - truelp)), 2, function(y) median(y)) # vector with nsim elements 
    globalmae <- mean(mean_over_test_obs)
    return(globalmae)
  }


fun_AUC <-
  function(obsy, predictions_response){
    
    auc <- apply(predictions_response, 1, function(x) pROC::auc(obsy, x, levels = c(0,1), direction = "<")) # vector (length nsim) of AUC values for each simulation iteration
    return(mean(auc))
  }

fun_ICI <- # for logistic regression: predictions should be predictions_response
  function(obsy, predictions){
    # This code follows the implementation of:  
    # Austin, PC, Steyerberg, EW. The Integrated Calibration Index (ICI) and related metrics 
    # for quantifying the calibration of logistic regression models. 
    # Statistics in Medicine. 2019; 38: 4051– 4065. https://doi.org/10.1002/sim.8281 
    
    # Estimate loess-based smoothed calibration curve and determine points on calibration curve
    counter = numeric(nrow(predictions))
    y.calibrate <- sapply(1:nrow(predictions), function(x) if (length(unique(predictions[x,])) > 2) {
      predict(loess(obsy ~ predictions[x,]), newdata = as.data.frame(predictions[x,]))
    } else {
      counter[x] <<- length(unique(predictions[x,]))
      if (length(unique(predictions[x,])) == 1) {
        z = rep(mean(obsy), ncol(predictions))
      } else if ((length(unique(predictions[x,])) == 2)) {
        z = ave(obsy, predictions[x,], FUN = mean)
      }
      z
    })
    y.calibrate <- t(y.calibrate)
    # rows: simulation iterations
    # columns: test observations, points on the loess calibration curve corresponding to a given predicted y (lin. regr.) or predicted probability (log. regr.)
    
    ICI <- apply(abs(y.calibrate - predictions), 1, mean) # vector (length nsim) of ICI values for each simulation iteration
    
    return(list(ICI = mean(ICI), counter = counter))
  }

##### ##### ##### ##### ##### ##### ##### #####
##### Functions for plotting measures #####
##### ##### ##### ##### ##### ##### ##### #####

# libraries
library(reshape2)
library(ggplot2)
library(plotly)
library(cowplot)
library(ggpubr)
library(grid)
library(RColorBrewer)   # for color palettes

# ggplot theme
my_theme <-function(){
  theme_bw() +
    theme(panel.grid = element_blank(), text = element_text(size=20),
          legend.position = "none") +
    theme(axis.text.x = element_text(size = 15))
}

col_low = "blue"
col_mid = "#CC79A7"
col_high = "#E69F00"

# legend

fun_show_legend <- function(col = T){
  multi_corr <- 0:5
  
  pmat <- data.frame(
    "x" = rep(rep(1:10, 6),2),          
    "y" = runif(0,1, n=10*6*2),
    "Multiple correlation" = rep(rep(1:6, each = 10),2),
    "group" = rep(1:12, each = 10),
    "line" = factor(rep(c("Predictor" ,"Noise variable"), each = 6*10))
  )  
  
  
  if(col == T){
    p1 <- ggplot(pmat, aes(x,
                           y,
                           linetype = line,
                           color = Multiple.correlation,
                           group = group)) +
      geom_line(size = 1) +
      #scale_colour_gradient(low = "blue", high = "red", 
      #                      labels = c( "low","","",  "", "","high")) +
      scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 3.5,
                            labels = c( "low","","",  "", "","high")) +
      scale_linetype_manual(breaks=c("Predictor","Noise variable"), values=c("solid","dashed")) +
      theme_bw() +
      labs(linetype = "", color = "Multiple correlation") +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            legend.text = element_text(size = 12)) + 
      guides( color = guide_colorbar(barwidth = 10, 
                                     barheight = 0.6, 
                                     title.vjust = 1, 
                                     title.theme = element_text(size = 12))) 
  }else{
    p1 <- ggplot(pmat, aes(x,
                           y,
                           linetype = line,
                           #color = Multiple.correlation,
                           group = group)) +
      geom_line(size = 1) +
      scale_colour_gradient(low = "blue", high = "red", 
                            labels = c( "low","","",  "", "","high")) +
      scale_linetype_manual(breaks=c("Predictor","Noise variable"), values=c("solid","dashed")) +
      theme_bw() +
      labs(linetype = "") +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            legend.text = element_text(size = 12)) 
  }
  
  
  legend_p1 <- get_legend(p1)
  
  #return(legend_p1)
}


# color scale
gg_color_hue = function(k) {
  hues = seq(15, 375, length = k + 1)
  hcl(h = hues, l = 65, c = 100)[1:k]
}

# method labels
method_labels <- c(
  "FU" = "FU",
  "BE_BIC" = "BE(BIC)",
  "BE_005" = "BE(0.05)",
  "BE_AIC" = "BE(AIC)",
  "BE_05" = "BE(0.50)",
  "ABE_AIC" = "ABE(AIC)",
  "FSel_AIC" = "FSel(AIC)", 
  "Step_FSel_AIC" = "Step_FSel(AIC)",
  "Lasso" = "Lasso(CV)",
  "RLasso" = "RLasso(CV)",
  "RLasso_BIC" = "RLasso(BIC)",
  "AdaLasso" = "AdaLasso(CV)",
  "Uni_005" = "Uni(0.05)",
  "Uni_020" = "Uni(0.20)", 
  "Uni_020_BE_005" = "Uni(0.20)+BE(0.05)"
)

all_methods <- c(
  "FU", "BE(BIC)", "BE(0.05)", "BE(AIC)", "BE(0.50)", "ABE(AIC)",
  "Lasso(CV)", "RLasso(CV)", "RLasso(BIC)", "AdaLasso(CV)",  "FSel(AIC)", "Step_FSel(AIC)",
  "Uni(0.05)", "Uni(0.20)", "Uni(0.20)+BE(0.05)"
)

palette_15 <- c(
  brewer.pal(8, "Dark2"),              # strong base: 8 colors
  brewer.pal(4, "Set1")[c(1,2,3,4)],   # vibrant reds & blues
  brewer.pal(3, "Set2")[1:3]           # balanced mid-saturation tones
)
master_palette_15 <- setNames(palette_15, all_methods)

master_shapes_15 <- setNames(0:14, all_methods)

# bias
fun_plot_bias <- function(results, methods, cond = F, aggr = F, estimand = 1, rescale_const = 1) {
  n <- sapply(results, function(x) x$n)
  if (estimand == 1) {
    beta_approx <- results[[1]]$beta_approx1
  } else if (estimand == 2) {
    beta_approx <- results[[1]]$beta_approx2
  }
  sdev <- results[[1]]$sd
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  res_bias <- list()
  for(i in 1:length(methods)){
    res_bias[[i]] <- sapply(results, function(x) fun_bias(beta_approx, x[[methods[i]]]$coefficients, sdev, cond = cond, rescale_const = rescale_const))[-1,]
    colnames(res_bias[[i]]) <- n
  }
  names(res_bias) <- methods
  if(aggr == T){
    res_bias <- lapply(res_bias, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(abs(z), na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res_bias) <- methods
    res_bias_long <- cbind(melt(res_bias), 
                           "pred" = rep(pred, length(n)*length(methods)))
    res_bias_long$Var2 <- as.factor(res_bias_long$Var2)
    res_bias_long$L1 <- factor(res_bias_long$L1, levels = methods)
    
    return(
      plot_grid(
        fun_show_legend(col = F),
        ggplot(res_bias_long, aes(x = Var2, y = value, group = Var1, linetype = pred))+
          geom_hline(yintercept = 0 , col = "grey") +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"bias")) + 
          my_theme() +
          scale_colour_gradient(low = "blue", high = "red") +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)) ,
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))
    )
    
  }else{
    res_bias_long <- cbind(melt(res_bias), 
                           "pred" = rep(pred, length(n)*length(methods)), 
                           "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                           "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_bias_long$Var2 <- as.factor(res_bias_long$Var2)
    res_bias_long$L1 <- factor(res_bias_long$L1, levels = methods)
    
    res_bias_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T), 
        ggplot(res_bias_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          geom_hline(yintercept = 0 , col = "grey") +
          geom_line() +
          geom_text_repel(aes(label = effect_size),
                          data = res_bias_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_bias_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"bias")) + 
          #scale_y_continuous(limits = c(-0.85, 0.4), breaks=c(-0.8,-0.6,-0.4, -0.2, 0, 0.2, 0.4)) + # ADDED for unconditional bias 
          #scale_y_continuous(limits = c(-1.1, 1.5), breaks=seq(-1,1.5,by=0.5)) + # ADDED for conditional bias 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1)) 
    )
  } 
}


# rmse * sqrt(n)
fun_plot_rmsen <- function(results, methods, cond = F, aggr = F, estimand = 1){
  n <- sapply(results, function(x) x$n)
  if (estimand == 1) {
    beta_approx <- results[[1]]$beta_approx1
  } else if (estimand == 2) {
    beta_approx <- results[[1]]$beta_approx2
  }
  sdev <- results[[1]]$sd
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  
  res <- list()
  for(i in 1:length(methods)){
    res[[i]] <- sapply(results, function(x) fun_rmsen(beta_approx, x[[methods[i]]]$coefficients, sdev, x$n, cond = cond))[-1,]
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  
  if(aggr == T){
    res <- lapply(res, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(z, na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res) <- methods
    res_long <- cbind(melt(res), 
                      "pred" = rep(pred, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    
    return(
      plot_grid( 
        fun_show_legend(col = F),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred))+
          geom_hline(yintercept = 0 , col = "grey") +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab(bquote(.(ifelse(cond == T, "Conditional", "Unconditional"))~~ sqrt(n)~RMSE)) + 
          my_theme() +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)) ,
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1)) 
    )
    
  }else{
    res_long <- cbind(melt(res), "pred" = rep(pred, length(n)*length(methods)), 
                      "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                      "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    res_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T), 
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          geom_hline(yintercept = 0 , col = "grey") +
          geom_line() +
          geom_text_repel(aes(label = effect_size),
                          data = res_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          xlab("Sample size") + 
          ylab(bquote(.(ifelse(cond == T, "Conditional", "Unconditional"))~~ sqrt(n)~RMSE)) + 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
  }
}


# coverage
fun_plot_cover <- function(results, methods, cond = F, aggr = F, estimand = 1, robustSEs = F){
  n <- sapply(results, function(x) x$n)
  if (estimand == 1) {
    beta_approx <- results[[1]]$beta_approx1
  } else if (estimand == 2) {
    beta_approx <- results[[1]]$beta_approx2
  }
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  
  res <- list()
  for(i in 1:length(methods)){
    if (robustSEs == T & methods[i] == "FU") {
      res[[i]] <- sapply(results, function(x) fun_cover(beta_approx, 
                                                        x[[methods[i]]]$CI.lower.robust, 
                                                        x[[methods[i]]]$CI.upper.robust, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        cond = cond))[-1,]
    } else {
      res[[i]] <- sapply(results, function(x) fun_cover(beta_approx, 
                                                        x[[methods[i]]]$CI.lower, 
                                                        x[[methods[i]]]$CI.upper, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        cond = cond))[-1,]
    }
    
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  
  
  
  if(aggr == T){
    res <- lapply(res, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(z, na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res) <- methods
    res_long <- cbind(melt(res), 
                      "pred" = rep(pred, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    
    return(
      plot_grid(
        fun_show_legend(col = F),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred))+
          geom_hline(yintercept = 0.95 , col = "grey") +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"Coverage")) + 
          my_theme() +
          scale_colour_gradient(low = "blue", high = "red") +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
    
  }else{
    res_long <- cbind(melt(res), "pred" = rep(pred, length(n)*length(methods)), 
                      "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                      "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    res_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          geom_hline(yintercept = 0.95 , col = "grey") +
          geom_line() +
          geom_text_repel(aes(label = effect_size),
                          data = res_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"Coverage")) + 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
  }
}


# CI width
fun_plot_width <- function(results, methods, cond = F, aggr = F, robustSEs = F){
  n <- sapply(results, function(x) x$n)
  beta_approx <- results[[1]]$beta_approx1
  sdev <- results[[1]]$sd
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  
  res <- list()
  for(i in 1:length(methods)){
    if (robustSEs == T & methods[i] == "FU") {
      res[[i]] <- sapply(results, function(x) fun_width(x[[methods[i]]]$CI.lower.robust, 
                                                        x[[methods[i]]]$CI.upper.robust, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        sdev = sdev, 
                                                        n = x$n,
                                                        cond = cond))[-1,]
    } else {
      res[[i]] <- sapply(results, function(x) fun_width(x[[methods[i]]]$CI.lower, 
                                                        x[[methods[i]]]$CI.upper, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        sdev = sdev, 
                                                        n = x$n,
                                                        cond = cond))[-1,]
    }
    
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  
  
  
  if(aggr == T){
    res <- lapply(res, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(z, na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res) <- methods
    res_long <- cbind(melt(res), 
                      "pred" = rep(pred, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    
    return(
      plot_grid(
        fun_show_legend(col = F),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred)) +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab(bquote(.(ifelse(cond == T, "Conditional", "Unconditional"))~~ sqrt(n)~CI~width)) + 
          my_theme() +
          scale_colour_gradient(low = "blue", high = "red") +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
    
  }else{
    res_long <- cbind(melt(res), "pred" = rep(pred, length(n)*length(methods)), 
                      "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                      "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    res_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          #geom_hline(yintercept = 0 , col = "grey") +
          geom_line() +
          geom_text_repel(aes(label = effect_size),
                          data = res_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          xlab("Sample size") + 
          ylab(bquote(.(ifelse(cond == T, "Conditional", "Unconditional"))~~ sqrt(n)~CI~width)) + 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)) ,
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1)) 
    )
  }
}


# Power
fun_plot_power <- function(results, methods, cond = F, aggr = F, robustSEs = F){
  n <- sapply(results, function(x) x$n)
  beta_approx <- results[[1]]$beta_approx1
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  
  res <- list()
  for(i in 1:length(methods)){
    if (robustSEs == T & methods[i] == "FU") {
      res[[i]] <- sapply(results, function(x) fun_power(x[[methods[i]]]$CI.lower.robust, 
                                                        x[[methods[i]]]$CI.upper.robust, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        cond = cond))[-1,]
    } else {
      res[[i]] <- sapply(results, function(x) fun_power(x[[methods[i]]]$CI.lower, 
                                                        x[[methods[i]]]$CI.upper, 
                                                        coefficients =  x[[methods[i]]]$coefficients, 
                                                        cond = cond))[-1,]
    }
    
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  
  
  if(aggr == T){
    res <- lapply(res, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(z, na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res) <- methods
    res_long <- cbind(melt(res), 
                      "pred" = rep(pred, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    
    return(
      plot_grid(
        fun_show_legend(col = F),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred)) +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"power/type-1 error")) + 
          my_theme() +
          scale_colour_gradient(low = "blue", high = "red") +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
    
  }else{
    res_long <- cbind(melt(res), "pred" = rep(pred, length(n)*length(methods)), 
                      "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                      "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    res_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          #geom_hline(yintercept = 0 , col = "grey") +
          geom_line() +
          geom_text_repel(aes(label = effect_size),
                          data = res_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          xlab("Sample size") + 
          ylab(paste(ifelse(cond == T, "Conditional", "Unconditional"),"power/type-1 error")) + 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
  }
}


# TPR & FPR
fun_plot_TPR_FPR <- function(results, methods, aggr = F){
  n <- sapply(results, function(x) x$n)
  beta_approx <- results[[1]]$beta_approx1
  pred <- (beta_approx==0)[-1]
  beta_sd <- abs(results[[1]]$beta_sd)
  multi_corr <- results[[1]]$multi_corr
  multi_corr_scaled <- round(multi_corr/max(multi_corr)*5,0)
  
  res <- list()
  for(i in 1:length(methods)){
    res[[i]] <- sapply(results, function(x) fun_TPR_FPR( coefficients =  x[[methods[i]]]$coefficients))[-1,]
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  
  if(aggr == T){
    res <- lapply(res, function(x) apply(x,2, function(y) aggregate(y, by = list(pred), FUN = function(z) mean(z, na.rm = T))[,"x"]))
    pred <- c(F,T)
    names(res) <- methods
    res_long <- cbind(melt(res), 
                      "pred" = rep(pred, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    
    return(
      plot_grid(
        fun_show_legend(col = F),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred)) +
          geom_line(size = 1.5) +
          xlab("Sample size") + 
          ylab("Variable selection rate") + 
          my_theme() +
          scale_colour_gradient(low = "blue", high = "red") +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
    
  }else{
    res_long <- cbind(melt(res), "pred" = rep(pred, length(n)*length(methods)), 
                      "multi_corr" = rep(multi_corr_scaled, length(n)*length(methods)),
                      "size" = rep(beta_sd+1, length(n)*length(methods)))
    res_long$Var2 <- as.factor(res_long$Var2)
    res_long$L1 <- factor(res_long$L1, levels = methods)
    res_long$effect_size <- rep(1:20, length(n)*length(methods))
    
    return(
      plot_grid(
        fun_show_legend(col = T),
        ggplot(res_long, aes(x = Var2, y = value, group = Var1, linetype = pred, color = multi_corr, size = size))+
          #geom_hline(yintercept = 0 , col = "grey") +
          geom_text_repel(aes(label = effect_size),
                          data = res_long %>% filter(effect_size <= 10 & Var2 == min(as.numeric(as.character(res_long$Var2)))),
                          nudge_x = -0.2,
                          size = 4)  +
          geom_line() +
          xlab("Sample size") + 
          ylab("Variable selection rate") + 
          my_theme() +
          scale_size_continuous(range = c(0.3,1.5)) +
          # scale_colour_gradient(low = "grey20", high = "grey80") +
          #scale_colour_gradient(low = "blue", high = "red") +
          scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 2.5,
                                labels = c( "low","","",  "", "","high")) +
          facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                      4, 3), 2),1)),
        nrow = 2, 
        rel_heights = c(0.1,0.9)*ifelse(length(methods)>=4,ifelse(length(methods)>=7, 
                                                                  ifelse(length(methods)>=10, 4, 3), 2),1))  
    )
  }
}

# Model size
fun_plot_model_size <- function(results, methods){
  n <- sapply(results, function(x) x$n)
  
  res <- list()
  for(i in 1:length(methods)){
    res[[i]] <- sapply(results, function(x) fun_model_size(coefficients = x[[methods[i]]]$coefficients))
    res[[i]] = as.data.frame(t(res[[i]]))
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  res_long <- melt(res)
  res_long$variable <- factor(res_long$variable, levels = n[order(n)])
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      scale_x_discrete(expand = expansion(add = c(1, 0.6))) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.5, size = 4, max.overlaps = 3)  +
      xlab("Sample size") +
      ylab("Nr of selected variables") +
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
}

# Kendall's tau_B
fun_plot_rank <- function(results, methods, estimand = 1){
  n <- sapply(results, function(x) x$n)
  if (estimand == 1) {
    beta_approx <- results[[1]]$beta_approx1
  } else if (estimand == 2) {
    beta_approx <- results[[1]]$beta_approx2
  }
  sdev <- results[[1]]$sd
  
  res <- list()
  for(i in 1:length(methods)){
    res[[i]] <- sapply(results, function(x) fun_rank(beta_approx = beta_approx, coefficients = x[[methods[i]]]$coefficients, sdev = sdev))
    res[[i]] = as.data.frame(t(res[[i]]))
    colnames(res[[i]]) <- n
  }
  names(res) <- methods
  res_long <- melt(res)
  res_long$variable <- factor(res_long$variable, levels = n[order(n)])
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.2,
                      size = 4)  +
      xlab("Sample size") +
      ylab("Kendall's tau_B") +
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
}


# True model
fun_plot_truemodel <- function(results, methods){
  n <- sapply(results, function(x) x$n)
  beta_approx <- results[[1]]$beta_approx1
  predictors <- beta_approx != 0
  
  res <- list()
  for(i in 1:length(methods)){
    res[[i]] <- sapply(results, function(x) fun_truemodel(coefficients =  x[[methods[i]]]$coefficients, predictors))
    colnames(res[[i]]) <- n
    rownames(res[[i]]) <- c("True model", "Over-selection model", "Under-selection model")
    res[[i]] <- res[[i]][c("Under-selection model", "Over-selection model", "True model"),]
  }
  names(res) <- methods
  res_long <- melt(res)
  res_long$Var2 <- as.factor(res_long$Var2)
  res_long$L1 <- factor(res_long$L1, levels = methods)
  
  return(
    ggplot(res_long, aes(x = Var2, y = value, fill = Var1))+
      geom_bar(position = position_stack(reverse = F), stat = "identity") +
      xlab("Sample size") + 
      ylab("Model selection rate") + 
      my_theme() +
      theme(legend.title = element_blank(), legend.position = "top") +
      guides(fill = guide_legend(reverse = T)) +
      scale_fill_manual(
        breaks = rev(c("Under-selection model", "Over-selection model", "True model")),
        values = rev(c("Under-selection model" = "#A31515", "Over-selection model" = "#EEE163", "True model" = "#1A9640"))) +
      facet_wrap(.~ L1, labeller = as_labeller(method_labels), nrow = ifelse(length(methods)>=4,ifelse(length(methods)>=7, ifelse(length(methods)>=10, 
                                                                                                                                  4, 3), 2),1)) 
  )
}

# y bias, y rmse (both on scale of linear predictor) 
# based on pre-calculated y_bias and y_rmse,
# i.e. based on pre-generated res_long (see results_..._preparation_for_shiny.R)
fun_plot_y_bias_rmse_prep <- function(res_long, methods, bias = T, calibration = F, log = FALSE){ # if bias = F: RMSE 
  
  res_long <- res_long[res_long$L1 %in% methods, ] # L1 denotes the method name 
  n <- as.numeric(levels(res_long$Var2))
  n <- n[order(n)]
  
  res_long$Var2 <-  factor(paste0("n=", res_long$Var2), # Var2 denotes the sample size 
                           levels = paste0("n=", n) )
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  
  
  #res_long$Var2 <-  factor(paste0("n=", res_long$Var2), # Var2 denotes the sample size 
  #                         levels = unique(paste0("n=", res_long$Var2) ))
  
  if (log) {
    if(calibration) {
      p<- ggplot(res_long, aes(x=pred_prob, y=obsy, color=L1)) + geom_abline(intercept=0, slope=1, col="grey") 
    }  else {
      p<- ggplot(res_long, aes(x=truelp, y=value, color=L1)) + geom_hline(yintercept = 0 , col = "grey") 
      # "value" is y-bias (see fun_y_bias) or y-RMSE (see fun_y_rmse)
    }
    return(
      p + 
        geom_smooth(se = F, method="loess") +
        scale_color_manual(values = color_values) +
        #  geom_point() + 
        xlab(ifelse(calibration,"Predicted probability", "True linear predictor")) + 
        ylab(ifelse(bias == T, ifelse(calibration,"Observed proportion","Bias"), expression(RMSE*sqrt(n)))) + 
        my_theme() + 
        facet_grid( . ~ Var2) +
        #      geom_hline(yintercept = 0 , col = "grey") +
        theme(legend.position = "top", legend.title = element_blank()) +
        geom_rug(data = data.frame("x" = quantile(res_long[res_long$Var2 == levels(res_long$Var2)[1],ifelse(calibration, "pred_prob", "truelp")], seq(0,1,0.05))), 
                 aes(x = x), sides = "b", col = "black", inherit.aes = F) 
      
    )
  } else {
    if(calibration&bias) {
      p<- ggplot(res_long, aes(x=value+truelp, y=obsy, color=L1)) + geom_abline(intercept=0, slope=1, col="grey") 
      # "value" is y-bias (difference between truelp and predicted y, see fun_y_bias) 
      # value + truelp is thus equal to predicted y
    }  else {
      p<- ggplot(res_long, aes(x=truelp, y=value, color=L1)) + geom_hline(yintercept = 0 , col = "grey") #+
      # geom_text_repel(aes(label = L1),
      #                 data = res_long %>% group_by(Var2) %>% filter(abs(truelp - mean(truelp)) == min(abs(truelp - mean(truelp)))),
      #                 nudge_x = -0.5, size = 4, max.overlaps = 3)
      # "value" is y-bias (see fun_y_bias) or y-RMSE (see fun_y_rmse)
    }
    return(
      #      ggplot(res_long, aes(x = ifelse(calibration&bias, truey-value, truey), y = ifelse(calibration,truey,value), color = L1)) +
      p + 
        geom_smooth(se = F, method=ifelse(calibration&bias,"lm","loess")) +
        scale_color_manual(values = color_values) +
        #  geom_point() + 
        xlab(ifelse(calibration,"Predicted Y", "True linear predictor")) + 
        ylab(ifelse(bias == T, ifelse(calibration,"Observed Y","Bias"), expression(RMSE*sqrt(n)))) + 
        my_theme() + 
        facet_grid( . ~ Var2) +
        #      geom_hline(yintercept = 0 , col = "grey") +
        theme(legend.position = "top", legend.title = element_blank()) +
        geom_rug(data = data.frame("x" = quantile(res_long[res_long$Var2 == "n=100" ,"truelp"], seq(0,1,0.05))), 
                 aes(x = x), sides = "b", col = "black", inherit.aes = F) 
      
    )
  }
  
}

# Global RMSE w.r.t. true vs. estimated linear predictor
fun_plot_global_rmse_prep <- function(res_long, methods) {
  res_long <- res_long[res_long$L1 %in% methods, ] # L1 denotes the method name 
  
  n <- as.numeric(levels(res_long$variable))
  n <- n[order(n)]
  
  res_long$variable <-  factor(res_long$variable, # variable denotes the sample size 
                               levels = n )
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.5, size = 4, max.overlaps = 3) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == max(n)),
                      nudge_x = 0.5, size = 4, max.overlaps = 3) +
      xlab("Sample size") +
      ylab(expression(RMSE*sqrt(n))) +
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
}

# Global MAE w.r.t. true vs. estimated linear predictor
fun_plot_global_mae_prep <- function(res_long, methods) {
  res_long <- res_long[res_long$L1 %in% methods, ] # L1 denotes the method name 
  
  n <- as.numeric(levels(res_long$variable))
  n <- n[order(n)]
  
  res_long$variable <-  factor(res_long$variable, # variable denotes the sample size 
                               levels =  n )
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      scale_x_discrete(expand = expansion(add = c(1, 0.6))) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.5, size = 4, max.overlaps = 3)  +
      xlab("Sample size") +
      ylab(expression(MAE)) + 
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
}

# AUC

fun_plot_auc_prep <- function(res_long, methods){
  res_long <- res_long[res_long$L1 %in% methods, ] # L1 denotes the method name 
  
  n <- as.numeric(levels(res_long$variable))
  n <- n[order(n)]
  
  res_long$variable <-  factor(res_long$variable, # variable denotes the sample size 
                               levels = n)
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      scale_x_discrete(expand = expansion(add = c(1, 0.6))) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.5, size = 4, max.overlaps = 3)  +
      xlab("Sample size") +
      ylab(expression(AUC)) + 
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
  
}


# Integrated Calibration Index (ICI) 

fun_plot_ici_prep <- function(res_long, methods){
  res_long <- res_long[res_long$L1 %in% methods, ] # L1 denotes the method name 
  
  n <- as.numeric(levels(res_long$variable))
  n <- n[order(n)]
  
  res_long$variable <-  factor(res_long$variable, # variable denotes the sample size 
                               levels = n)
  
  res_long$L1 = factor(res_long$L1, levels = methods)
  levels(res_long$L1) = method_labels[methods]
  
  color_values <- master_palette_15[levels(res_long$L1)]
  shape_values <- master_shapes_15[levels(res_long$L1)]
  
  return(
    ggplot(res_long, aes(x = variable, y = value, group = L1)) +
      geom_line(aes(color = L1), size = 1.2) +
      geom_point(aes(color = L1, shape = L1), size = 2.5, stroke = 1.5) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = shape_values) +
      scale_x_discrete(expand = expansion(add = c(1, 0.6))) +
      geom_text_repel(aes(label = L1),
                      data = res_long %>% filter(variable == min(n)),
                      nudge_x = -0.5, size = 4, max.overlaps = 3)  +
      xlab("Sample size") +
      ylab(expression(ICI)) + 
      theme(text = element_text(size=20), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size=20),
            axis.text = element_text(size = 20)) +
      guides(linetype = guide_legend(override.aes = list(linewidth = 2))) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  )
  
}