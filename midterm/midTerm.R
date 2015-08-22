# rm(list=ls())
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
load('midTermData.rda')
library(Amelia)

olsSM = function(formula, data, impute=FALSE){
	# Organize iv and dvs
	dv=all.vars(formula)[1]
	ivs=all.vars(formula)[2:length(all.vars(formula))]

	# Missingness strategy
	if(impute){
		set.seed(6886)
		data = amelia(data[,c(dv,ivs)], m=1)$imp$imp1
	}
	data = na.omit(data)

	# Create matrix with column for intercept and
	## data from independent variables
	y = data[,dv]
	x = data.matrix(cbind(1, data[,ivs]))

	# General parameters
	n = nrow(x) # Number of observations
	p = length(ivs) # Number of parameters
	df = n-p-1 # degrees of freedom

	# Betas
	xTx = t(x) %*% x
	xTy = t(x) %*% y
	beta = solve(xTx) %*% xTy

	# Standard errors
	yhat = x %*% beta
	e = y - yhat
	eTe = sum( e^2 )
	sigma2 = eTe/(df)
	seBeta = sqrt( sigma2 * diag( solve(xTx) ) )

	# Calculate t values
	tStat = beta/seBeta

	# Calculate p-values
	pval = 2 * pt(abs(tStat), df = df, lower.tail = FALSE)

	# Varcovar
	varcovBeta = sigma2 * solve(xTx)
	colnames(varcovBeta)=rownames(varcovBeta)=c('(Intercept)', ivs)

	# R squared
	ssReg = sum( (yhat - mean(y))^2 )
	ssTot = sum( (y - mean(y))^2 )
	Rsq = ssReg/ssTot

	# Fstatistic statemtn
	muSqReg = sum(yhat^2)/p
	muSqRes = sum(e^2)/( df )
	Fstat = muSqReg/muSqRes
	Fstatp = pf(Fstat, p, df, lower.tail=FALSE)
	Fresult = paste0("F-statistic: ", round(Fstat,3), 
		" on ", p, " and ", df, " DF, p-value: ", round(Fstatp,3))

	# Org coef matrix
	lower  = beta - qt(.975, df) * seBeta
	upper  = beta + qt(.975, df) * seBeta
	coefMat = cbind(beta, seBeta, tStat, pval, lower, upper)
	colnames(coefMat)=c("Estimate", "Std. Error", "T-Statistic", 
		"P-Value", 'Lower 95% CI', 'Upper 95% CI')
	rownames(coefMat)=c('(Intercept)', ivs)

	# All results
	return( list(
		coefficients = coefMat,
		varcov = varcovBeta,
		Rsq = Rsq,
		Fstat = Fresult
		) )
}

# set.seed(6886)
# form = formula(gini_net_std ~ polity2 + ELF_ethnic)

# olsSM(form, data)$coefficients
# round(summary(lm(form, data=data))$coefficients,3)

# round(ols(form, dataMiss)$coefficients,3)
# round(summary(lm(form, data=dataMiss))$coefficients,3)

# round(ols(form, dataMiss, impute=TRUE)$coefficients,3)
# # dataAmelia = amelia(dataMiss, m=1)$imp$imp1
# # round(summary(lm(form, data=dataAmelia))$coefficients,3)