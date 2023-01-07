# Script Author: Brenton Graham
# Last Edit: 12/15/2022

# Import libraries
require(tidyverse)
require(hdrm)
require(olsrr)
require(SignifReg)
require(MASS)
require(glmnet)
require(coefplot)
require(magrittr)
require(RcmdrMisc)
require(caret)
require(gmodels)

# Helper functions -------------------------------------------------------------
# Function to extract and clean model summary
modelSummary <- function(model){
  summary <- summary(model)$coefficients %>%
    as.data.frame() %>% 
    set_colnames(c("estimate", "stderr", "tval", "pval")) %>% 
    dplyr::select(-tval)
  return(summary)
}

# Function to merge results into template resultsDF
mergeResults <- function(resultsDF, simResults, modelType, sim){
  mergedDF <- resultsDF %>%
    merge(simResults, all.x = T, all.y = F, by = "row.names") %>% 
    mutate(modelType = modelType, sim = sim) %>%
    dplyr::select(modelType, sim, variable = Row.names, everything())
  return(mergedDF)
}

# Function to determine bias, coverage and power from model results
resultsCalc <- function(resultsDF, alpha){
  finalResults <- resultsDF %>%
    mutate(estimate_binary = case_when(is.na(estimate) ~ 0, TRUE ~ 1),
           true_beta_binary = case_when(true_beta == 0 ~ 0, TRUE ~ 1),
           bias = as.numeric(estimate - true_beta),
           coverage = as.numeric((estimate - 1.96*stderr <= true_beta) * (true_beta <= estimate + 1.96*stderr)),
           power = as.numeric(case_when((pval <= alpha) ~ 1, (pval > alpha) ~ 0)))
  return(finalResults)
}

# Function to determine mean estimate and CI for given column in DF
confInt <- function(DF, metric){
  value_to_report <- ci(DF %>% dplyr::select(metric) %>% drop_na() %>% pull) %>%
    as.matrix() %>% t() %>% 
    as.data.frame() %>% 
    set_colnames(c(paste(metric, "_estimate", sep = ""), paste(metric, "_ci_lower", sep = ""), paste(metric, "_ci_upper", sep = ""), "stderr")) %>% 
    #mutate(est_and_ci = paste(sprintf("%.3f", estimate), " (", sprintf("%.3f", ci_lower), ", ", sprintf("%.3f", ci_upper), ")", sep = "")) %>%
    dplyr::select(-stderr)
  return(value_to_report)
}

# Function to report TPR, FPR, bias, coverage and power estimates for model
report <- function(DF, n, rho, modelType){
  confusionMatrix <- confusionMatrix(data = factor(DF$estimate_binary), reference = factor(DF$true_beta_binary), positive = "1")
  tpr <- confusionMatrix$byClass["Sensitivity"]
  fpr <- 1 - confusionMatrix$byClass["Specificity"]
  bias <- confInt(DF, "bias")
  coverage <- confInt(DF, "coverage")
  power <- confInt(DF, "power")
  reportDF <- data.frame(modelType = modelType, n = n, rho = rho, tpr = tpr, fpr = fpr, bias, coverage, power) %>%
    set_rownames(NULL)
  return(reportDF)
}

selectionRates <- function(DF, nsim, modelType, n, rho){
  returnDF <- DF %>% 
    group_by(variable) %>%
    summarise(selection_rate = sum(!is.na(estimate))/nsim) %>%
    ungroup() %>%
    mutate(modelType = modelType, n = n, rho = rho) %>%
    dplyr::select(modelType, n, rho, everything())
  return(returnDF)
}

# Simulation -------------------------------------------------------------------

# Set up DF with changing parameter configurations
dataConfigurations <- data.frame(
  config = seq(1, 6),
  n = c(rep(250,  3), rep(500,  3)),
  rho = rep(c(0.0, 0.4, 0.8), 2))

# Set fixed parameters
nsim <- 100
true.beta <- c(seq(0.5/3, 2.5/3, 0.5/3), rep(0, 15))

# Set up empty DFs to store all metrics from every simulation
allResultsDF <- data.frame()
selectionRatesDF <- data.frame()

# For each data configuration
for (j in seq(1, nrow(dataConfigurations))){
  
  # Set configuration-specific parameter combination
  n <- dataConfigurations[j, ]$n
  rho <- dataConfigurations[j, ]$rho
  
  # Set up DF with all variables to merge variable selection results into
  resultsDF <- data.frame(
    variables = c("V01", "V02", "V03", "V04", "V05", "V06", "V07", "V08", "V09", "V10",
                  "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20"),
    true_beta = true.beta) %>%
    column_to_rownames("variables")
  
  # Set up empty DFs to store results from each model type in simulation
  resultsDF.pval <- data.frame()
  resultsDF.aic <- data.frame()
  resultsDF.bic <- data.frame()
  resultsDF.lasso_minLambda <- data.frame()
  resultsDF.lasso_1seLambda <- data.frame()
  resultsDF.elasticNet_minLambda <- data.frame()
  resultsDF.elasticNet_1seLambda <- data.frame()
  
  # Set seed for reproducibility
  set.seed(0)
  
  # Perform simulation
  for (i in seq(1, nsim)){
    
    # Generate Data ------------------------------------------------------------
    data <- genData(n = n, p = 20, p1 = 5, beta = true.beta, rho = rho)
    df <- cbind(y = data$y, data$X) %>% as.data.frame()
    
    # Models -------------------------------------------------------------------
    # Full linear regression model; need for backward selection methods
    full.model <- lm(y ~ ., data = df)
    
    # Backward selection variable selection - p-value, AIC, BIC
    ## p-value
    step_pval.model <- ols_step_backward_p(full.model, prem = 0.15)
    step_pval.summary <- modelSummary(step_pval.model$model)
    step_pval.mergedResults <- mergeResults(resultsDF, step_pval.summary, modelType = "step_pval", sim=i)
    resultsDF.pval <- rbind(resultsDF.pval, step_pval.mergedResults)
    
    ## AIC
    step_aic.model <- stepAIC(full.model, direction = "backward", trace = F)
    step_aic.summary <- modelSummary(step_aic.model)
    step_aic.mergedResults <- mergeResults(resultsDF, step_aic.summary, modelType = "step_aic", sim=i)
    resultsDF.aic <- rbind(resultsDF.aic, step_aic.mergedResults)
    
    ## BIC
    step_bic.model <- stepwise(full.model, direction = "backward", criterion = "BIC", trace = F)
    step_bic.summary <- modelSummary(step_bic.model)
    step_bic.mergedResults <- mergeResults(resultsDF, step_bic.summary, modelType = "step_bic", sim=i)
    resultsDF.bic <- rbind(resultsDF.bic, step_bic.mergedResults)
    
    # Model selection with LASSO and elastic net
    ## LASSO
    ### lambda.min
    lasso_minLambda <- cv.glmnet(
      x = as.matrix(df %>% dplyr::select(-"y")), y = df$y, family = "gaussian", 
      alpha = 1, type.measure = "mse", nfolds = 5)
    lasso_minLambda.coef <- extract.coef(lasso_minLambda, s = "lambda.min")
    lasso_minLambda.vars <- lasso_minLambda.coef %>% filter(Coefficient != "(Intercept)") %>% row.names()
    # Extract selected variables and fit new model to obtain estimates
    lasso_minLambda_DF <- df %>% dplyr::select(y, lasso_minLambda.vars)
    lasso_minLambda.lrmod <- lm(y ~ ., data = lasso_minLambda_DF)
    lasso_minLambda.summary <- modelSummary(lasso_minLambda.lrmod)
    lasso_minLambda.mergedResults <- mergeResults(resultsDF, lasso_minLambda.summary, modelType = "lasso_minLambda", sim=i)
    resultsDF.lasso_minLambda <- rbind(resultsDF.lasso_minLambda, lasso_minLambda.mergedResults)
    
    ### lambda.1se
    lasso_1seLambda <- cv.glmnet(
      x = as.matrix(df %>% dplyr::select(-"y")), y = df$y, family = "gaussian", 
      alpha = 1, type.measure = "mse", nfolds = 5)
    lasso_1seLambda.coef <- extract.coef(lasso_1seLambda, s = "lambda.1se")
    lasso_1seLambda.vars <- lasso_1seLambda.coef %>% filter(Coefficient != "(Intercept)") %>% row.names()
    # Extract selected variables and fit new model to obtain estimates
    lasso_1seLambda_DF <- df %>% dplyr::select(y, all_of(lasso_1seLambda.vars))
    lasso_1seLambda.lrmod <- lm(y ~ ., data = lasso_1seLambda_DF)
    lasso_1seLambda.summary <- modelSummary(lasso_1seLambda.lrmod)
    lasso_1seLambda.mergedResults <- mergeResults(resultsDF, lasso_1seLambda.summary, modelType = "lasso_1seLambda", sim=i)
    resultsDF.lasso_1seLambda <- rbind(resultsDF.lasso_1seLambda, lasso_1seLambda.mergedResults)
    
    ## Elastic net
    ### lambda.min
    elasticNet_minLambda <- cv.glmnet(
      x = as.matrix(df %>% dplyr::select(-"y")), y = df$y, family = "gaussian", 
      alpha = 0.5, type.measure = "mse", nfolds = 5)
    elasticNet_minLambda.coef <- extract.coef(elasticNet_minLambda, s = "lambda.min")
    elasticNet_minLambda.vars <- elasticNet_minLambda.coef %>% filter(Coefficient != "(Intercept)") %>% row.names()
    # Extract selected variables and fit new model to obtain estimates
    elasticNet_minLambda_DF <- df %>% dplyr::select(y, elasticNet_minLambda.vars)
    elasticNet_minLambda.lrmod <- lm(y ~ ., data = elasticNet_minLambda_DF)
    elasticNet_minLambda.summary <- modelSummary(elasticNet_minLambda.lrmod)
    elasticNet_minLambda.mergedResults <- mergeResults(resultsDF, elasticNet_minLambda.summary, modelType = "elasticNet_minLambda", sim=i)
    resultsDF.elasticNet_minLambda <- rbind(resultsDF.elasticNet_minLambda, elasticNet_minLambda.mergedResults)
    
    ### lambda.1se
    elasticNet_1seLambda <- cv.glmnet(
      x = as.matrix(df %>% dplyr::select(-"y")), y = df$y, family = "gaussian", 
      alpha = 0.5, type.measure = "mse", nfolds = 5)
    elasticNet_1seLambda.coef <- extract.coef(elasticNet_1seLambda, s = "lambda.1se")
    elasticNet_1seLambda.vars <- elasticNet_1seLambda.coef %>% filter(Coefficient != "(Intercept)") %>% row.names()
    # Extract selected variables and fit new model to obtain estimates
    elasticNet_1seLambda_DF <- df %>% dplyr::select(y, all_of(elasticNet_1seLambda.vars))
    elasticNet_1seLambda.lrmod <- lm(y ~ ., data = elasticNet_1seLambda_DF)
    elasticNet_1seLambda.summary <- modelSummary(elasticNet_1seLambda.lrmod)
    elasticNet_1seLambda.mergedResults <- mergeResults(resultsDF, elasticNet_1seLambda.summary, modelType = "elasticNet_1seLambda", sim=i)
    resultsDF.elasticNet_1seLambda <- rbind(resultsDF.elasticNet_1seLambda, elasticNet_1seLambda.mergedResults)
  }
  
  # Determine bias, coverage and power from model results
  reportDF.pval <- resultsCalc(resultsDF.pval, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "step_pval")
  reportDF.aic <- resultsCalc(resultsDF.aic, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "step_aic")
  reportDF.bic <- resultsCalc(resultsDF.bic, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "step_bic")
  reportDF.lasso_minLambda <- resultsCalc(resultsDF.lasso_minLambda, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "lasso_minLambda")
  reportDF.lasso_1seLambda <- resultsCalc(resultsDF.lasso_1seLambda, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "lasso_1seLambda")
  reportDF.elasticNet_minLambda <- resultsCalc(resultsDF.elasticNet_minLambda, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "elasticNet_minLambda")
  reportDF.elasticNet_1seLambda <- resultsCalc(resultsDF.elasticNet_1seLambda, alpha = 0.05) %>% report(n = n, rho = rho, modelType = "elasticNet_1seLambda")
  
  # Concatenate results into one DF
  allResultsDF.tmp <- rbind(
    reportDF.pval, reportDF.aic, reportDF.bic, reportDF.lasso_minLambda, reportDF.lasso_1seLambda, reportDF.elasticNet_minLambda, reportDF.elasticNet_1seLambda)
  
  # Concatenate configuration-specific results into final DF
  allResultsDF <- rbind(allResultsDF, allResultsDF.tmp)
  
  # Determine selection rates of each variable
  selectionRatesDF.pval <- selectionRates(resultsDF.pval, nsim = nsim, modelType = "step_pval", n = n, rho = rho)
  selectionRatesDF.aic <- selectionRates(resultsDF.aic, nsim = nsim, modelType = "step_aic", n = n, rho = rho)
  selectionRatesDF.bic <- selectionRates(resultsDF.bic, nsim = nsim, modelType = "step_bic", n = n, rho = rho)
  selectionRatesDF.lasso_minLambda <- selectionRates(resultsDF.lasso_minLambda, nsim = nsim, modelType = "lasso_minLambda", n = n, rho = rho)
  selectionRatesDF.lasso_1seLambda <- selectionRates(resultsDF.lasso_1seLambda, nsim = nsim, modelType = "lasso_1seLambda", n = n, rho = rho)
  selectionRatesDF.elasticNet_minLambda <- selectionRates(resultsDF.elasticNet_minLambda, nsim = nsim, modelType = "elasticNet_minLambda", n = n, rho = rho)
  selectionRatesDF.elasticNet_1seLambda <- selectionRates(resultsDF.elasticNet_1seLambda, nsim = nsim, modelType = "elasticNet_1seLambda", n = n, rho = rho)
  
  # Concatenate results into one DF
  selectionRatesDF.tmp <- rbind(
    selectionRatesDF.pval, selectionRatesDF.aic, selectionRatesDF.bic, selectionRatesDF.lasso_minLambda, selectionRatesDF.lasso_1seLambda, selectionRatesDF.elasticNet_minLambda, selectionRatesDF.elasticNet_1seLambda)
  
  # Concatenate configuration-specific results into final DF
  selectionRatesDF <- rbind(selectionRatesDF, selectionRatesDF.tmp)
}

write.table(allResultsDF, "simulation_results.csv", sep = ",", row.names = F)
write.table(selectionRatesDF, "var_selection_results.csv", sep = ",", row.names = F)

