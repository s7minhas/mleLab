rm(list = ls())

# library(Amelia)
# library(ggplot2)

# setwd('~/OneDrive/MLE/Midterm/')

# load('midTermData.rda') #loads data & dataMiss

# form = formula(gini_net_std ~ ELF_ethnic + polity2)

# test = lm(gini_net_std ~ ELF_ethnic + polity2, data = data)

ols <- function(formula, data, impute=FALSE){
	if(impute){
		set.seed(6886)
		data = amelia(x = data, m = 1)
		data = data$imputations[[1]]
	}

	else{
		data = data[complete.cases(data), ]
	}


	dvs = (all.vars(formula)[1])
	ivs = (all.vars(formula)[2:length(all.vars(formula))])
	x = as.matrix(cbind(1, data[,ivs]))
	y = as.matrix(data[,dvs])
	feed = round((solve(t(x)%*%x)%*%t(x)%*%y), 5)

	#building the residulas
	#This will give me my resid but is not broadly applicable (work on that)
	res = as.matrix(data[,1] - feed[1] - feed[2]*data[,3] - feed[3]*data[,2])

	# establishing dimensions
	n = nrow(x)
	k = ncol(x)
	p2 = (k-1)
	df = n - p2 - 1

	#VCV
	varcov = 1/(n-k)*as.numeric(t(res)%*%res) * solve(t(x)%*%x)
	row.names(varcov) = c('Intercept', ivs)

	serror = round(sqrt(diag(varcov)), 5) #these are the standard errors

	t = round(rbind((feed[1]/serror[1]), (feed[2]/serror[2]), (feed[3]/serror[3])),3)

	p = round(rbind(2*pt(abs(t[1]), df = n-k, lower.tail = FALSE),
	2*pt(abs(t[2]), df = n-k, lower.tail = FALSE), 
	2*pt(abs(t[3]), df = n-k, lower.tail = FALSE)), 4)

	ci = NULL
	for(ii in 1:length(feed)){
		ci = rbind(ci, round(feed[ii] + c(-1,1)*serror[ii]*qt(.95, df),5))
	}

	#this is the final data, I guess i should make it a matrix before turning it in
	feed_f = as.data.frame(cbind(feed, serror, t, p, ci[,1], ci[,2]))
	colnames(feed_f) = c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
	row.names(feed_f) = c('Intercept', ivs)
	coefficients = as.matrix(feed_f)

	#R^2

	yhat = feed[1] + feed[2]*data[3] + feed[3]*data[2]

	y_mean = mean(y)

	ssreg = sum((yhat-y_mean)^2)
	sstot = sum((y - y_mean)^2)

	Rsq = round(ssreg/sstot, 4)

	adj_Rsq = round(1-(1-Rsq)*((n-1)/df), 4)

	#F-stat
	e = y - yhat

	msreg = sum(yhat^2)/p2
	msres = sum(e^2)/df

	Fstat = round(msreg/msres, 3)

	pwhole = round(pf(Fstat, p2, df, lower.tail = FALSE), 3)

	Fstat = paste("F-statistic:", Fstat, "on", p2, "and", df, "DF, p-value:", pwhole, sep = " ", collapse = NULL)

	#final
	final = list(coefficients = coefficients, varcov = varcov, Rsq = Rsq, Fstat = Fstat)
		return(final)
}

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

# #the models

# model = ols(formula=form, data=data)

# modelListDel = ols(formula = form, data=dataMiss)

# modelAmelia = ols(formula = form, data=dataMiss, impute=TRUE)

# #general plot

# gg = ggplot(data, aes(x = ELF_ethnic, y = gini_net_std)) + geom_point(col="firebrick") + theme_bw()
# gg = gg + geom_abline(aes(colour = "Original"), intercept = model$coefficients[1], 
# 	slope = model$coefficients[2], size = 1, show_guide = TRUE) #complete data
# gg = gg + geom_abline(aes(colour = "Listwise Delete"), intercept = modelListDel$coefficients[1], 
# 	slope = modelListDel$coefficients[2], size = 1)#listwise delete
# gg = gg + geom_abline(aes(colour = "Imputation"), intercept = modelAmelia$coefficients[1], 
# 	slope = modelAmelia$coefficients[2], size = 1) #Amelia
# gg = gg + scale_colour_manual(values = c('#00BA38','#F8766D', '#00bfc4'))
# gg = gg + xlab('ELF Ethnic') + ylab('GINI') + theme(legend.position="none")
# gg = gg + ggtitle('OLS Results')
# gg

# ggsave(gg, file = 'comp.png')

# #coef plot
# model_2 = as.data.frame(model$coefficients)
# modelListDel_2 = as.data.frame(modelListDel$coefficients)
# modelAmelia_2 = as.data.frame(modelAmelia$coefficients)

# row.names(model_2) = c('(Intercept)', 'ELF_ethnic', 'polity2')
# row.names(modelListDel_2) = c('(Intercept)', 'ELF_ethnic', 'polity2')
# row.names(modelAmelia_2) = c('(Intercept)', 'ELF_ethnic', 'polity2')

# colnames(model_2) = c('Coefficient', 'SE', 'T', 'P', 'Lower', 'Upper')
# colnames(modelListDel_2) = c('Coefficient', 'SE', 'T', 'P', 'Lower', 'Upper')
# colnames(modelAmelia_2) = c('Coefficient', 'SE', 'T', 'P', 'Lower', 'Upper')

# model_2frame = data.frame(Variable = rownames(model_2),
# 					Coefficient = model_2[,1],
# 					SE = model_2[,2],
# 					modelName = 'Original')

# modelListDel_2frame = data.frame(Variable = rownames(modelListDel_2),
# 					Coefficient = modelListDel_2[,1],
# 					SE = modelListDel_2[,2],
# 					modelName = 'Listwise Delete')

# modelAmelia_2frame = data.frame(Variable = rownames(modelAmelia_2),
# 					Coefficient = modelAmelia_2[,1],
# 					SE = modelAmelia_2[,2],
# 					modelName = 'Imputation')

# allModelFrame = data.frame(rbind(model_2frame, modelListDel_2frame, modelAmelia_2frame))

# interval1 <- -qnorm((1-0.9)/2)  # 90% multiplier
# interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# zp1 <- ggplot(allModelFrame, aes(colour = modelName))
# zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
# zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
#                                 ymax = Coefficient + SE*interval1),
#                             lwd = 1, position = position_dodge(width = 1/2))
# zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
#                                  ymax = Coefficient + SE*interval2),
#                              lwd = 1/2, position = position_dodge(width = 1/2),
#                              shape = 21, fill = "WHITE")
# zp1 <- zp1 + coord_flip() + theme_bw()
# zp1 <- zp1 + ggtitle("Comparing models") + scale_colour_manual(values = c('#F8766D', '#00bfc4', '#00BA38'))
# zp1

# ggsave(zp1, file = 'coef.png')