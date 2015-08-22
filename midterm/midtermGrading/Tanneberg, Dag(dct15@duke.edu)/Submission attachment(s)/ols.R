rm(list=ls())
ols <- function(formula, data, impute = FALSE){
  ## Author: Dag Tanneberg
  ## This function implements OLS regression analysis for numeric
  ## variable objects. It computes regression coefficients, 
  ## standard errors, and corresponding test statistics. 
  ## Moreover, it calculates the R^2 statistic and conducts an 
  ## omnibus F-test on the model. 
  ## 
  ## The function takes a regression formula and data frame as 
  ## inputs and returns a list of regression results. A logical 
  ## flag 'impute' may be set to TRUE to estimate missing data
  ## with setting amelia(data, m=1). Default behavior is to 
  ## drop all observations with missing values.
  ## 
  ## formula ... An object of class formula specifying dv and iv's
  ## data ... An object of class data.frame
  ## impute ... A logical flag for imputing missing data. Default 
  ##  is FALSE.

  ## --- Data Input Objects -------------------------------------
  dta <- data[, all.vars(formula)]   ## subset data to regression
  # n.org <- nrow(dta)       ## test loop response to impute flag
  if(impute == TRUE){             ## imputation flag and response
    set.seed(6886)
    dta <-  Amelia::amelia(dta, m = 1)[['imputations']][['imp1']]
      ## Note: Process breaks if there are completely missing obs
  } else {
    dta <- na.omit(dta)            ## listwise deletion otherwise
  }
  X <- as.matrix(cbind(1, dta[, all.vars(formula)[-1]]))
  Y <- as.vector(dta[, all.vars(formula)[1]])
  n <- nrow(dta)
  # n.mis <- n.org-n         ## test loop response to impute flag
  degFree <- n-ncol(X)
  
  ## --- Define output object -----------------------------------
  output <- list(
    coefficients = matrix(NA, nrow = ncol(X), ncol = 6),  # betas
    varcov = matrix(NA, nrow = ncol(X), ncol = ncol(X)), # vcov b
    Rsq = vector(length = 1),                               # Rsq
    Fstat = vector(length = 1)                     # F-Statistic
    # data = dta ## return regression data, required for some plots
  #  , n.mis = n.mis         ## test loop response to impute flag
  )
  colnames(output[['coefficients']]) <- c(
    'Estimate', 'Std. Error', 'T-Statistic', 'P-Value', 
    'Lower 95% CI', 'Upper 95% CI'
  )
  rownames(output[['coefficients']]) <- c('(Intercept)', 
    all.vars(formula)[-1]
  )

  ## --- Compute Beta hat ---------------------------------------
  output[['coefficients']][, 1] <- solve(t(X)%*%X)%*%t(X)%*%Y

  ## --- Model fit basics ---------------------------------------
  fitted <- X%*%output[['coefficients']][, 1] ## yhat predictions
  residuals <- fitted-Y
  sigma2 <- 1/degFree*sum(diag(residuals%*%t(residuals)))
  
  ## --- Regressor Variance-Covariance Matrix -------------------
  output[['varcov']] <- sigma2*solve(t(X)%*%X)
  rownames(output[['varcov']]) <- c('(Intercept)', colnames(X)[-1])
  colnames(output[['varcov']]) <- c('(Intercept)', colnames(X)[-1])

  ## --- Define 1. Coefficient SE, 2. t-values, and 3. p-values -
  output[['coefficients']][, 2] <- sqrt(diag(output[['varcov']]))
  output[['coefficients']][, 3] <- output[['coefficients']][, 1] / 
    output[['coefficients']][, 2]
  output[['coefficients']][, 4] <- sapply( ## return vec of pvals
    output[['coefficients']][, 3], function(x) { 
      2*(1-pt(abs(x), df = degFree)) 
    }
  )

  ## --- .95 Confidence intervals -------------------------------
  output[['coefficients']][, 5] <- output[['coefficients']][, 1] - 
    qt(p = .975, df = degFree) * output[['coefficients']][, 2]
  output[['coefficients']][, 6] <- output[['coefficients']][, 1] + 
    qt(p = .975, df = degFree) * output[['coefficients']][, 2]
  
  ## --- Model fit sum of squares -------------------------------
  sst <- sum((Y-mean(Y))^2)               ## total sum of squares
  ssm <- sum((fitted-mean(Y))^2)          ## model sum of squares
  sse <- sum(residuals^2)              ## residual sum of squares 

  ## --- R-Squared statistic ------------------------------------
  output[['Rsq']] <- 1-(sse/sst)

  ## --- Omnibus F-Statistic ------------------------------------
  num <- ssm/ (ncol(X)-1)
  den <- sse/(n-(ncol(X)))
  fStat <- num/den
  pval <- 1-pf(abs(fStat), (ncol(X)-1), (n-(ncol(X))))
  output[['Fstat']] <- paste(
    'F-Statistic:', 
    format(round(fStat, digits = 3), nsmall = 3), 'on', 
    (ncol(X)-1), 'and', n-ncol(X), 'DF, p-value:', 
    format(round(pval, digits = 3), nsmall = 3), sep = ' '
  )

  ## --- Return results -----------------------------------------
  return(output)
}
## End

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

olsSM(form, data)
ols(form, data=data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))