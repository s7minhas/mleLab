# Start with a clean workspace
rm(list=ls())

# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
	if(! lib %in% installed.packages()[,1])
	  { install.packages(lib, repos='http://cran.rstudio.com/') }
	suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c("ggplot2", 'lmtest', 'car', 'gap', 'sandwich')
loadPkg(packs)

# Set a theme for gg
theme_set(theme_bw())

# Functions that I use frequently
char = function(x){ as.character(x) }
num = function(x){ as.numeric(char(x)) }

# Relevant paths
labPath='~/Dropbox/Duke/Spring 2015/PS 733/lab3'

# Muller & Seligson Regression
# Load data
msrepPath=paste0(labPath, "/Msrepl87.asc")
msrep = read.table(msrepPath, header=TRUE)

# Create silly logged version of DV
msrep$deaths75ln = log(msrep$deaths75+1)

# Create logs of other things
msrep$deaths70ln = log(msrep$deaths70+1)
msrep$sanctions75ln = log(msrep$sanctions75+1)
msrep$sanctions70ln = log(msrep$sanctions70+1)
msrep$energypcln = log(msrep$energypc+1)

## Functional Form

### Omitted Variable Bias
msrepBiRegData=na.omit(msrep[,c('deaths75ln', 'upper20', 'sanctions75ln')])
cor(msrepBiRegData)

# Estimate coefficients
coeftest( lm(deaths75ln ~ upper20 + sanctions75ln, data=msrepBiRegData) )

# If we excluded sanctions75ln from the model what would you expect of the
# coefficient on upper20
coeftest( lm(deaths75ln ~ upper20, data=msrepBiRegData) )

# If we excluded upper20 from the model what would you expect of the
# coefficient on sanctions75ln
coeftest( lm(deaths75ln ~ sanctions75ln, data=msrepBiRegData) )

# More complicated model
ivs=c('upper20', 'energypcln', 'intensep', 
      'sanctions70ln', 'sanctions75ln', 'deaths70ln')
msrepRegData=na.omit(msrep[,c('deaths75ln', ivs)])
olsForm1=formula(paste0('deaths75ln ~ ', paste(ivs, collapse=' + ')))
olsForm2=formula(paste0('deaths75ln ~ ', paste(ivs[c(2:length(ivs))], collapse=' + ')))
mod1 = lm(olsForm1, data=msrepRegData)
mod2 = lm(olsForm2, data=msrepRegData)

# View model results
summary(mod1)
summary(mod2)

### F-test
# How would we do a partial F-test
ssr1=sum(resid(mod1)^2)
ssr2=sum(resid(mod2)^2)
chgReg=ncol(model.matrix(mod1)) - ncol(model.matrix(mod2))
# F statistic
Fstat=((ssr1-ssr2)/chgReg)/(ssr1/df.residual(mod1))
1-pf(abs(Fstat), chgReg, df.residual(mod1))

# Using anova function from base R
anova(mod1, mod2)

### Likelihood Ratio Test
# How would we run a likelihood ratio test
ltest=nrow(model.matrix(mod2))*log(ssr1/ssr2)
pchisq(abs(ltest), df=chgReg, lower.tail=FALSE)

# Using lrtest from lmtest library
lrtest(mod1, mod2)

### Partial Residual Plots
# How to test for non-linear effects in parameters
# Can use a the crPlots function from the car library
# to generate partial residual plots
crPlots(mod1)

# Incorporating a non-linear effect
olsForm3=formula('deaths75ln ~ poly(upper20,2) + energypcln + intensep + 
                 sanctions70ln + sanctions75ln + deaths70ln')
mod3=lm(olsForm3, data=msrepRegData)

# Examine coefficient results
coeftest(mod3)

# How to evaluate whether to keep non-linear term?
anova(mod1, mod3)
lrtest(mod1, mod3)

### Heteroskedasticity
# How to calculate the Breusch-Pagan test statistic
residMod1=resid(mod1)^2
bpForm=formula(paste0('residMod1 ~', paste(ivs, collapse=' + ')))
bpMod=lm(bpForm, data=msrepRegData)
bpStat=summary(bpMod)$r.squared*nrow(msrepRegData)
1-pchisq(bpStat, df=length(ivs))

# Breusch-Pagan test: using bptest from lmtest library
bptest(mod1)

#### Ways of testing: Visualization
# Residuals vs fitted
par(mfrow=c(1,2))
plot(mod1, 1)
plot(predict(mod1), resid(mod1))

# Scale-location
stdResid=(resid(mod1)-mean(resid(mod1)))/sd(resid(mod1))
plot(mod1, 3)
plot(predict(mod1), sqrt(abs(stdResid)))

# Other diagnostic plots
par(mfrow=c(2,3)); for(ii in 1:6){plot(mod1, ii)}

#### Dealing with heteroskedasticity: Robust standard errors
# Can use the sandwich package in combination with lmtest to help
# us calculate robust standard errors
coeftest(mod1, vcov=vcovHC(mod1, type='HC1'))
coeftest(mod1)