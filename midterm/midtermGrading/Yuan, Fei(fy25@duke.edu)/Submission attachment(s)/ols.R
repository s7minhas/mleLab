rm(list=ls())
ols = function(formula, data, impute = FALSE){
  if(impute==TRUE){
  set.seed(6886)
  data=amelia(x=data,m=1)$imp$imp1}
  data=data[complete.cases(data),]
  dv = all.vars(formula)[1]
  ivs = all.vars(formula)[2:length(all.vars(formula))]
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  n = nrow(x) 
  p = length(ivs) 
  df = n-p-1 
  
  ## calculate beta 
xTx=t(x) %*% x
xTy=t(x)%*%y
beta=solve(xTx)%*%xTy
beta

  
 ## calculate standard error 
yhat=x%*%beta
e=y-yhat
eTe=t(e)%*%e
SSR=eTe
residual2=SSR/n-p-1
res = as.matrix(residual2)
n = nrow(x)
k = ncol(x)
varcov = 1/(n-k) * as.numeric(t(res)%*%res) * solve(t(x)%*%x)
se= sqrt(diag(varcov))

##tstat, p-val and confidence interval 
tstat = beta/se
pval=round(2*pt(abs(tstat),df,lower.tail=FALSE),digits=3)
up95=beta+qt(0.975,df)*se
lo95=beta-qt(0.975,df)*se
  

## Get R square
SSreg=sum((yhat-mean(y))^2)
SStot=sum((y-mean(y))^2)
Rsq=SSreg/SStot
                
## Get F stat
MSREG=sum(yhat^2)/p
MSRES=sum(e^2)/df
Fstat=round(MSREG/MSRES,3)
Fstatp = round(pf(Fstat,p,df,lower.tail=F),3)


## Get R-sq
SSreg=sum((yhat-mean(y))^2)
SStot=sum((y-mean(y))^2)
Rsq=SSreg/SStot

## list
coefficients = cbind(beta,se,tstat,pval,lo95,up95)
colnames(coefficients)=c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
rownames(coefficients)=c('(Intercept)',ivs)

Fstat=paste("F-statistic:",Fstat,"on",p,"and",df,"DF,","p-value:",Fstatp,sep=' ') 

list = list(coefficients=coefficients,varcov=varcov,Rsq=Rsq,Fstat=Fstat)
return(list)
  }

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

ols(form, data=data)
olsSM(form, data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))