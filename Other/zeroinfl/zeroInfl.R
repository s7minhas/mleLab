rm(list=ls())

library(pscl) # loads zeroinfl library

load("~/Dropbox/Duke/Spring 2015/PS 733/Other/zeroinfl/zeroInfl.rda")

plot(density(datanew$Y))

# run a simple Poisson and get the predicted values for the actual observations in the dataset
m1 <- glm(Y ~ X + C1 + C4 + C5, data=datanew, family=poisson)
summary(m1)
XBm1 <- cbind(1, datanew$X, datanew$C1, datanew$C4, datanew$C5) %*% m1$coef
predm1 <- exp(XBm1)

# run a Poisson model correcting for possible overdispersion
m2 <- glm(Y ~ X + C1 + C4 + C5, data=datanew, family=quasipoisson)
summary(m2)

# run a zero-inflated Poisson with X and C1 for the count part and C4 and C5 for the zero-inflated part
m3 <- zeroinfl(Y ~ X + C1 | C4 + C5, data=datanew, dist="poisson", link="logit")
summary(m3)

# for each observation, what is the probability that it is zero-inflated?
ZG <- cbind(1,datanew$C4, datanew$C5) %*% m3$coef$zero
probzero <- 1/(1+exp(-ZG))
check <- predict(m3, type="zero")
cor(check,probzero)

# what is the predicted count, assuming it can have more than 0 (without the zero-inflated part)?
XB <- cbind(1,datanew$X, datanew$C1) %*% m3$coef$count
predcount <- exp(XB)
check <- predict(m3, type="count")
cor(check,predcount)

# what is the predicted count, all together?
pred <- (1-probzero)*predcount
check <- predict(m3, type="response")
cor(check,pred)

#plot predicted versus actual. For Poisson and Zero-Inflated. Which one does better? Compute some fit statistics
plot(datanew$Y, pred, pch=16)
abline(0,1)
plot(datanew$Y, predm1, pch=16)
abline(0,1)

# Absolute Prediction error
mean(abs(datanew$Y- pred))
mean(abs(datanew$Y- predm1))

# Squared Prediction error
mean((datanew$Y- pred)^2)
mean((datanew$Y- predm1)^2)

# what would the predicted count be if X1 was 0.2 lower for every observation (using the zero-inflated model)? Plot against the predicted values with the true X1
ZGexp <- cbind(1,datanew$C4, datanew$C5) %*% m3$coef$zero
probzeroexp <- 1/(1+exp(-ZGexp))

XBexp <- cbind(1,datanew$X-0.2, datanew$C1) %*% m3$coef$count
predcountexp <- exp(XBexp)

predexp <- (1-probzeroexp)*predcountexp

plot(pred, predexp, pch=16)
abline(0,1)