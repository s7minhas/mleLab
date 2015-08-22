# Nominal Lecture
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/nominalLecture")
set.seed(6886)

# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
	if(! lib %in% installed.packages()[,1])
	  { install.packages(lib, repos='http://cran.rstudio.com/') }
	suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c('foreign', 'MASS', 
	'ggplot2', 'reshape2', 'nnet')
loadPkg(packs)

# Helpful functions
char = function(x){ as.character(x) }
###################################################

###################################################
krain<-read.dta("isq05.dta",convert.dates = TRUE,
	convert.factors = TRUE, missing.type = FALSE,
	convert.underscore=TRUE, warn.missing.labels=TRUE)

mnlMod<-multinom(magnitud ~ intrvlag+icntglag+maglag
	+genyr+stfl+regtype+ethkrain+marg+coldwar,
	data=krain)

summary(mnlMod)
###################################################

###################################################
# no let’s make factors of the factor variables
krain$magnitud<-as.ordered(krain$magnitud)
z.out<-polr(magnitud ~
        intrvlag+icntglag+maglag+genyr+stfl+regtype+ethkrain+marg+coldwar,
        data=krain, method="logistic", Hess=T)
summary(z.out)
xtable(summary(z.out)$coefficients,digits=2)
###################################################