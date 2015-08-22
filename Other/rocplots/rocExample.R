############################
### Here's how you use the new function
library(ggplot2)
# Load new function
roc = function(prediction, actual){
      library(ROCR)
      pred = prediction(prediction, actual)
      perf = performance(pred,"tpr","fpr")
      rocData = data.frame(attributes(perf)$x.values[[1]], attributes(perf)$y.values[[1]])
      names(rocData) = c('FPR', 'TPR')
      return(rocData)
}

getAUC = function(prediction, actual){
	library(ROCR)
	pred = prediction(prediction, actual)	
	attributes(performance(pred,"auc"))$y.values[[1]]
}

# Load your test file
load('~/Dropbox/Duke/Spring 2015/PS 733/Other/rocPlots/test.rda')
probs = cbindtest[,1]
dv = cbindtest[,2]

# Calculate roc curve using new function
rocCurve = roc(probs, dv)
rocPlot=ggplot(rocCurve, aes(x=FPR, y=TPR)) + geom_line() + geom_abline(intercept=0, slope=1)
rocPlot
############################

############################
## Comparison to old approach
# old method
rocOld = function(threshold){
	# Set up output matrix
	pr=matrix(NA, ncol=3, nrow=length(threshold), 
		dimnames=list(NULL, c('Threshold', 'FPR', 'TPR') ) )

	# Loop through thresholds
	for(ii in 1:length(threshold)){
		predRoc = as.numeric( probs > threshold[ii] )
		FPR = sum( (predRoc==1)*(dv==0) ) / sum(dv==0)
		TPR = sum( (predRoc==1)*(dv==1) ) / sum(dv==1)
		pr[ii,1] = threshold[ii]
		pr[ii,2] = FPR
		pr[ii,3] = TPR
	}

	# Return output
	return( data.frame( pr ) )
}

# Calculate roc curve using old function
thresh = seq(0,1,.0001)
rocCurveOld = rocOld(thresh)

# Merge data
ggData = rbind( 
	cbind( rocCurve, new=1 ),
	cbind( rocCurveOld[,2:3], new=0 )
 )

# Calculate AUCs
aucNew = getAUC(probs, dv)
i = 2:nrow(rocCurveOld)
aucOld = -1*(rocCurveOld$FPR[i] - rocCurveOld$FPR[i - 1]) %*% (rocCurveOld$TPR[i] + rocCurveOld$TPR[i - 1])/2

# Add labels
ggData$new[ggData$new==1] = paste0(
	'Updated \n Approach \n (AUC = ', round(aucNew,3), ')')
ggData$new[ggData$new==0] = paste0(
	'Old \n Approach \n (AUC = ', round(aucOld,3), ')')

# Plot old and new side by side
rocPlot = ggplot(ggData, aes(x=FPR, y=TPR))
rocPlot = rocPlot + geom_point() + geom_abline(intercept=0, slope=1)
rocPlot = rocPlot + facet_wrap(~new)
rocPlot
############################