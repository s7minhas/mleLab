# First set a seed
set.seed(6886)
rm(list=ls())
ols <- function(formula, data, impute=FALSE){

    if (impute==TRUE){

        #Another iteration- run regressions on m imputations for(i in i:m){}
        dvi = all.vars(formula)[1]
        ivsi = all.vars(formula)[ 2:length(all.vars(formula)) ]
                                        #y = data[,dv]
                                        #x = data.matrix(cbind(1, data[,ivs]))
     ## Impute data
     set.seed(6886)
        imputed <- amelia(x=data[,c(dvi, ivsi)], m=1)
        regdata <- as.data.frame(imputed$imp)
        colnames(regdata) <- c(dvi, ivsi)
        #rename columns so that it works with the ols function
        data <- as.matrix(regdata) #takes the imputed data, turns it into "data"
        print("Data has been imputed")
    }

    dta <- na.omit(data)
    check <- dim(dta)==dim(data)
    nrm <- dim(data)-dim(dta)
    if(check[1]==FALSE){print(paste0("Warning: ", nrm[1]," rows deleted during data prep.",
                "  Consider imputation"))}


    ##Step 2: create matrix with columns for intercept and IDVs:
    dv = all.vars(formula)[1]
    ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
        y = dta[,dv]
        x = data.matrix(cbind(1, dta[,ivs])) 
    colnames(x)[1] <- c('(Intercept)')

        ## here is where it becomes regular again, and I can wrap it into another function
    n = nrow(x) # Number of observations
    p = length(ivs) # Number of parameters
    df = n-p-1 # degrees of freedom

## coefficients:
xtx <-  t(x) %*%x
xty <- t(x) %*% y
## calculating Beta
betahat <-  solve(xtx) %*% xty
# Standard errors
yhat <- x %*% betahat
resid <- y- yhat
epe <- sum(resid^2)
sigmasq <- epe/df

varcov <- sigmasq*solve(xtx) #checks out so far

##elements of coefficients conf intervals, t-stat, pvals
stderrors <-  sqrt(diag(varcov))
upper95 = betahat + qt(.975,df)*stderrors
lower95 = betahat - qt(.975,df)*stderrors
tstats <- betahat/stderrors
pval <-  2 * pt(abs(tstats), df=n-p-1, lower.tail=F)

coefficients <- (cbind(betahat, stderrors,
                 tstats, pval,
                 lower95, upper95))
colnames(coefficients) <- c("Estimate", "Std. Error",
                                "T-Statistic", "P-Value",
                                "Lower 95% CI",
                                "Upper 95% CI")

# Rsquared:
SSreg <-  sum((yhat-mean(y))^2)
SStot <- sum((y-mean(y))^2)
Rsq <-  SSreg/SStot
# names(Rsq) <- "R-squared"

##F-statistic
    MSreg <-  sum(yhat^2)/ p
    MSresid <-  sum(resid^2)/df
    Fstat <-  MSreg/ MSresid  #Double check this- it is the
    #place that I'm not on the same page as the lm output
    pval2 <- pf(Fstat, p, df, lower.tail=FALSE)

    FSoutput <- paste0("F-statistic: ", format(Fstat, digits=5),
                          " on ", p ," and ", df, " DF,",
                  " p-value: ", round(pval2,3))

outList <- list(coefficients, varcov, Rsq, FSoutput)
names(outList) <- c('coefficients', 'varcov', 'Rsq', "Fstat")
return(outList)
}


library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

ols(form, data)
olsSM(form, data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))