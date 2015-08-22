
#Replication data for:
#Berliner, Daniel. 2014. The Political Origins of Transparency. The Journal of Politics, 76(2): 479-491.

setwd("~/Desktop") 

data<-read.csv("FOI_passage_finaldata.csv", header=T)

#variable definitions
#CODE country code
#YEAR year
#DEVELOPING indicator for developing countries
#FOIPASS year of FOI law passage
#turnfreq frequency of turnover
#opp1vote vote share of largest opposition party
#oppvote vote share of all opposition parties
#opp1vote_dev change over time in opp1vote
#turnfreq1 turnover frequency in categories
#FOI_REGIONAL neighborhood context of FOI passage
#FOI_IGO inter-governmental organization context of FOI passage
#polity2 Polity2 democracy index
#NEWDEM new democracy indicator
#logINGO log INGO memberships
#CORRUPTION Transparency International CPI, average over time
#logtrade log trade per GDP
#loggdpcap log GDP per capita
#AID_DEPENDENCE log foreign aid per GDP
#IMFDUMMY use of IMF credit
#FHDEM Freedom House democracy measure
#checks veto players measure from DPI
#PRESSFREEDOM Freedom House freedom of the press
#PRES presidential systems from DPI
#FDISTOCK stock of inward FDI, from UNCTAD
#GDPGRTH GDP growth
#BATTLEDEATHS battle deaths
#INJUD judicial independence from CIRI
#yrcurnt years in office of current ruling party, from DPI


options(scipen=10)

#this function from Chris Adolph simcf package, available on his website
extractdata<-function (formula, data, extra = NULL, na.rm = FALSE) 
{
    subdata <- get_all_vars(formula, data)
    if (!is.null(extra)) {
        if (is.data.frame(extra)) {
            subdata <- cbind(subdata, extra)
        }
        else {
            if (any(class(extra) == "formula")) {
                subdata <- cbind(subdata, get_all_vars(extra, 
                  data))
            }
            else {
                stop("Extra must be a dataframe or formula object")
            }
        }
    }
    if (na.rm) 
        subdata <- na.omit(subdata)
    subdata
}


#data is all years all countries
#data_d is all years developing countries only
#data_s is post90 all countries
#data_ds is post90 developing only
data_d<-data[data$DEVELOPING==1,]
data_s<-data[data$YEAR>1989,]
data_ds<-data_d[data_d$YEAR>1989,]


#set which sample to use
data_x<-data_ds


###Models in main text of paper:

# .... removing code for models 1 through 4

#model 4: both
formula1<-FOIPASS ~ 
polity2 +
turnfreq +
opp1vote + 
NEWDEM +
logINGO +
CORRUPTION +
logtrade +
loggdpcap +
AID_DEPENDENCE +
IMFDUMMY  +
FOI_REGIONAL+
FOI_IGO

mdata1<-extractdata(formula1, data_x, extra=data_x[,c("CODE","YEAR")], na.rm = TRUE)
uc1<-unique(as.character(mdata1$CODE))
mdata1$yearsatrisk<-rep(NA,nrow(mdata1))
for(i in 1:length(uc1)){mdata1$yearsatrisk[mdata1$CODE==uc1[i]]<- 1:(length(mdata1$yearsatrisk[mdata1$CODE==uc1[i]])) }

formula1<-paste("FOIPASS","~", 
"polity2 +opp1vote + turnfreq +FOI_REGIONAL+ FOI_IGO+ NEWDEM  +logINGO+CORRUPTION +loggdpcap +logtrade + AID_DEPENDENCE +IMFDUMMY", sep="")

yearsatrisk_dummies <- model.matrix(~factor(mdata1$yearsatrisk)-1)
yearsatrisk_dummies<-yearsatrisk_dummies[,-1]
selectdata <- cbind(mdata1,yearsatrisk_dummies)
ddnames <- NULL
for (i in 1:ncol(yearsatrisk_dummies)) {
  formula1 <- paste(formula1,"+dd",i,sep="")
  ddnames <- c(ddnames,paste("dd",i,sep=""))
}
names(selectdata) <- c(names(mdata1),ddnames)

m4<-glm(as.formula(formula1), data=selectdata, family=binomial(link="logit"), control=glm.control(maxit=500))
summary(m4)

###############################################################
### My attempt to replicate figure 2

# Create scenario matrix for figure 2
data08 = selectdata[selectdata$YEAR==2008,] # Subset to 2008, so we can get mean 2008 levels
vars = names(coef(m4)); vars=vars[-1] # Pull out coefficients used in model 4, take out intercept
simData = apply(data08[,vars], 2, function(x){ mean(x, na.rm=TRUE) }) # Calculate mean of every variable in 2008
simData[paste0('dd',1:18)] = 0 # Set 
simData['dd18'] = 1

turnFreqRange = seq(0,4,.25) # Set up varying values for turnover frequency
scen = NULL ; for(ii in 1:length(turnFreqRange)){ scen=rbind(scen,simData) } # Expand scenario vector into matrix
scen[,'turnfreq'] = turnFreqRange # Replace mean turn over frequency value with varying values from 0 to 4
scen = cbind(1, scen) # Add intercept term
print(scen)

# Run sims
library(MASS) # load mass library for mvrnorm function
set.seed(6886) # set seed
sims=10000 # number of simulations
draws = mvrnorm(sims, coef(m4), vcov(m4)) # sample from multivariate normal
yVals = draws %*% t(scen) # Calculate predicted values
yProbs = 1/(1+exp(-yVals)) # Calculate predicted probabilities by applying logistic link function

# Get summary stats by scenario
info = function(x){ c( mean(x), quantile(x, probs=c(0.025, 0.975)) )  } # Function to calculate mean and 95% interval for every scenario
ySumm = t(apply(yProbs, 2, info)) # Apply function above to columns of yProbs matrix
colnames(ySumm)=c('mu','lo','hi') # Label columns

# Plot for gg
ggData = data.frame(turnFreqRange, ySumm, row.names=NULL) # Convert to dataframe for plotting

library(ggplot2) # load ggplot library
# Plot
tmp=ggplot(ggData, aes(x=turnFreqRange, y=mu, ymin=lo, ymax=hi))
tmp=tmp + geom_line() + geom_ribbon(alpha=.5)
tmp=tmp + ylab('Probability of FOI Passage') + xlab('Turnover Frequency')
tmp
###############################################################