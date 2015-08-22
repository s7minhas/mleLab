                 round(tstats, 3), pval,
                 lower95, upper95))
colnames(coefficients) <- c("Estimate", "Std. Error",
                                "T-Statistic", "P-Value",
                                "Lower 95% CI",
                                "Upper 95% CI")

# Rsquared:
SSreg <-  sum((yhat-mean(y))^2)
SStot <- sum((y-mean(y))^2)
Rsq <-  SSreg/SStot
names(Rsq) <- "R-squared"

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
