# ---
# title: "Midterm Assignment"
# author: "Jaewon Chung"
# date: "March 17, 2015"
# output:
# pdf_document:
#   fig_caption: yes
# word_document: default
# header-includes:
# - \usepackage{multirow}
# - \usepackage{dcolumn}
# ---


# # Appendix: R codes for Midterm Assignment

# ## 0. Setting Up the Workspace  

# ```{r, message=FALSE, warning=FALSE}

# # Set up workspace
# rm(list=ls()) 
# setwd('C:/Users/Jaewon Chung/Desktop/Duke University/15 Spring/MLE/HWMID/')
# set.seed(6886)

# # Function to load packages
# loadPkg=function(toLoad){
#   for(lib in toLoad){
#     if(! lib %in% installed.packages()[,1])
#     { install.packages(lib, repos='http://cran.rstudio.com/') }
#     suppressMessages( library(lib, character.only=TRUE) ) }
# }

# # Load libraries
# packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop')
# loadPkg(packs)

# ```

# ## 1. Creating an OLS function
# ```{r, message=FALSE, warning=FALSE}

# # Load data
# load('midTermData.rda')

# Creating a function
rm(list=ls())
ols <- function(formula, data, impute=FALSE){
  
  
  if(impute){
    set.seed(6886)
    dataAmelia = amelia(x=data, m=1)
    names(dataAmelia$imp)
    data=dataAmelia$imp$imp1
    miss=FALSE
  }
  
  data=na.omit(data)
  
  # Retrieve vars from formula input
  dv = all.vars(form)[1]
  ivs = all.vars(form)[ 2:length(all.vars(form)) ]
  
  y = data[,dv]
  x = data.matrix(cbind('(Intercept)', data[,ivs]))
  i = diag(1, nrow=nrow(x), ncol=ncol(x))
  
  
  # General parameters
  n = nrow(x)                  # Number of observations
  p = length(ivs)              # Number of parameters
  
  
  # Calculate coefficient estimates
  xy = t(x) %*% y              # x`y
  xxi = solve(t(x) %*% x)      # (x`x)^-1
  h = x %*% xxi %*% t(x)       # hat mtrix of x
  i = diag(1, nrow=n, ncol=n)
  b = as.vector(xxi %*% xy)    # estimated coefficients 
  
  # Calculating standard errors
  
  yhat = as.vector(x %*% b)    # predicted values for y
  res = y-yhat   # model residuals
  
  sst = sum((y-mean(y))^2)     # total sum of squres 
  ssr = sum(res^2)             # residual sum of squares
  ssm = sst-ssr                # sum of squares for model
  
  df.e = (n-p-1)               # dgrees of freedom for error
  df.t = (n-1)                 # total degrees of freedom
  df.m = df.t-df.e             # degrees of freedom for model
  
  s2 = as.vector(ssr/df.e)     # (sigma hat)^2
  sigma2 = ssr/(n-p) 
  
  # Variance-covariance
  varcov = s2 * xxi
  serrors = sqrt(diag(varcov)) # coefficients standard errors
  
  # R-squared
  Rsq = 1-(ssr/sst)           
  
  # F-statisics
  f = (ssm/df.m)/(ssr/df.e) 
  fpvalue = pf(f, df.m, df.e, lower.tail =FALSE)
  
  # Calculate confidence intervals around regression estimates
  ## upper 95% CI and lower 95% CI
  
  up95=b+qt(.975, df.e)*serrors
  lo95=b-qt(.975, df.e)*serrors
  
  # t-statistics for st.errors
  tstats = b/serrors     
  
  # p-values
  pval = 2*pt(abs(tstats), df.e, lower.tail = F)
  
  # Coefficients matrix
  b.table = cbind(Estimate = b, 'Std. Error' = serrors, 'T-Statistic' = tstats, 'P-Value' = pval, 
                  'Lower 95% CI' = lo95, 'Upper 95% CI' = up95)
  
  rownames(b.table)[1]=c('(Intercept)')
  rownames(varcov)=colnames(varcov)=rownames(b.table)
  return(list(coefficients=b.table, varcov=varcov, Rsq=Rsq, 
              Fstat=paste("F-statistic:", f=round(f,3), "on", df.m,"and", df.e, "DF,", 
                          "p-value:", fpvalue=round(fpvalue,3) ))
  )
}
# 
# ```
library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)


expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
summary(lm(form, data=dataMiss))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))

# ## 2. Running models
# ```{r, message=FALSE, warning=FALSE}

# # Set up the model formula
# form = formula(gini_net_std ~ ELF_ethnic + polity2)

# # Run the various models by using the function
# model = ols(formula=form, data=data)
# modelListDel = ols(formula = form, data=dataMiss)
# modelAmelia = ols(formula = form, data=dataMiss, impute=TRUE)

# ```
# \
# \

# \begin{tabular}{ccccccc}
# \hline \hline
# & Estimate & Std.Error & T-Statistic & P-Value & Lower 95\% CI & Upper 95\% CI \\ 
# \hline \hline
# "(Intercept)" & $$-$0.220$ & $0.854$ & $$-$0.258$ & $0.798$ & $$-$8.808$ & $8.368$ \\ 
# ELF\_ethnic & $2.010$ & $0.631$ & $3.184$ & $0.003$ & $$-$6.577$ & $10.598$ \\ 
# polity2 & $$-$0.096$ & $0.080$ & $$-$1.204$ & $0.235$ & $$-$8.684$ & $8.492$ \\ 
# \hline \hline
# \end{tabular}

# \
# \

# \begin{tabular}{ccccccc}
# \hline \hline
# & Estimate & Std.Error & T-Statistic & P-Value & Lower 95\% CI & Upper 95\% CI \\ 
# \hline \hline
# "(Intercept)" & $0.669$ & $0.569$ & $1.175$ & $0.247$ & $$-$2.359$ & $3.698$ \\ 
# ELF\_ethnic & $0.805$ & $0.444$ & $1.816$ & $0.078$ & $$-$2.223$ & $3.834$ \\ 
# polity2 & $$-$0.165$ & $0.052$ & $$-$3.157$ & $0.003$ & $$-$3.193$ & $2.864$ \\ 
# \hline \hline
# \end{tabular}

# \
# \

# \begin{tabular}{ccccccc}
# \hline \hline
# & Estimate & Std.Error & T-Statistic & P-Value & Lower 95\% CI & Upper 95\% CI \\ 
# \hline \hline
# "(Intercept)" & $1.938$ & $0.460$ & $4.215$ & $0.0001$ & $$-$1.912$ & $5.788$ \\ 
# ELF\_ethnic & $0.334$ & $0.452$ & $0.739$ & $0.464$ & $$-$3.516$ & $4.185$ \\ 
# polity2 & $$-$0.287$ & $0.038$ & $$-$7.493$ & $0$ & $$-$4.137$ & $3.563$ \\  
# \hline \hline
# \end{tabular}
