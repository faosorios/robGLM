# load Finney's dataset
load("../data/finney.rda")

# source .R files and load required packages
source("../code/robLogistic.R")
source("../code/qopt.R")
source("../code/envelope.R")
library("robustbase")

# fit logistic regression with q = 1
fm100 <- robLogistic(Response ~ log(Volume) + log(Rate), data = finney, qpar = 1.0)

# selection of distortion parameter
o <- qopt(Response ~ log(Volume) + log(Rate), data = finney)
qpar <- o$qpar

# fit logistic regression with optimal q = 0.79
fm079 <- robLogistic(Response ~ log(Volume) + log(Rate), data = finney, qpar = qpar)

# fit logistic regression removing observations 4 and 18
fm.rm <- robLogistic(Response ~ log(Volume) + log(Rate), data = finney, subset = -c(4,18), qpar = 1.0)

# fit logistic regression using robust procedures
fmCR <- glmrob(Response ~ log(Volume) + log(Rate), family = binomial, data = finney, method = "Mqle", control = glmrobMqle.control(tcc=3.5))
fmBY <- glmrob(Response ~ log(Volume) + log(Rate), family = binomial, data = finney, method= "BY")
fmWBY <- glmrob(Response ~ log(Volume) + log(Rate), family = binomial, data = finney, method= "WBY")

# Table 1
tab1 <- matrix(0, nrow = 5, ncol = 6)
tab1[1,] <- c(t(summary(fm100)$coef[,1:2]))
tab1[2,] <- c(t(summary(fm.rm)$coef[,1:2]))
tab1[3,] <- c(t(summary(fm079)$coef[,1:2]))
tab1[4,] <- c(t(summary(fmCR)$coef[,1:2]))
tab1[5,] <- c(t(summary(fmBY)$coef[,1:2]))
rownames(tab1) <- c("ML", "ML.remove", "MLq", "CR", "BY")
colnames(tab1) <- c("Intercept", "SE", "log(volume)", "SE", "log(rate)", "SE")

tab1
#          Intercept      SE  log(volume)       SE  log(rate)       SE
#ML          -2.8754   1.3207      5.1793   1.8647     4.5617   1.8379
#ML.remove  -24.5812  14.0200     39.5498  23.2444    31.9352  17.7581
#MLq         -5.1848   2.5631      8.2344   3.9197     7.2874   3.4546
#CR          -2.7534   1.3268      4.9739   1.8620     4.3881   1.8445
#BY          -6.8515  10.0403     10.7343  15.2948     9.3643  12.7701

# estimated probabilities
eta100 <- fm100$linear
eta079 <- fm079$linear
pred <- cbind(eta100, eta079, finney$Response)
colnames(pred) <- c("linear100", "linear079", "Response")
pred <- as.data.frame(pred)
oeta100 <- sort(eta100)
oeta100 <- c(-25, oeta100, 10)
prob100 <- exp(oeta100) / (1 + exp(oeta100))
oeta079 <- sort(eta079)
oeta079 <- c(-25, oeta079, 10)
prob079 <- exp(oeta079) / (1 + exp(oeta079))

# Fig 1.a
pdf()
par(pty = "s")
plot(oeta100, prob100, xlim = c(-25,10), ylim = c(0,1.05), type = "l", xlab = "Linear predictor", ylab = "Probability", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
plot(Response ~ linear100, data = pred, subset = (Response == 0), ylim = c(0,1.05), xlim = c(-25,10), pch = 1, axes = FALSE, xlab = "", ylab = "", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
plot(Response ~ linear100, data = pred, subset = (Response == 1), ylim = c(0,1.05), xlim = c(-25,10), pch = 16, axes = FALSE, xlab = "", ylab = "", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
text(eta100[c(4,24)], finney$Response[c(4,24)], as.character(c(4,24)), pos = 3, cex = 1.2)
text(eta100[18], finney$Response[18], as.character(18), pos = 1, cex = 1.2)
dev.off()

# Fig 1.b
ylim <- c(-7, 7)
n <- nrow(finney)
qq100 <- envelopeq.bin(fm100, reps = 500, conf = 0.95)
std <- qq100$std
pdf()
par(pty = "s")
qqnorm(std, ylim = ylim, main = "", ylab = "Studentized residuals", xlab = "Quantiles of standard normal", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
qqnorm(qq100$envelope[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(qq100$envelope[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
z <- qnorm(ppoints(n))[order(order(std))]
# order(std) reveals 'locations' of observations 4, 18 and 24
text(sort(z)[c(1,38,39)], sort(std)[c(1,38,39)], labels = as.character(c(24,18,4)), pos = 3, cex = 1.2)
dev.off()

# Fig 2
obs <- c(4,18,24)
pdf()
par(pty = "s")
plot(fm079$rob.weights, ylim = c(0,1), ylab = "Estimated weights", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
abline(h = 1, lwd = 2, lty = 2, col = "red")
text(obs, fm079$rob.weights[obs], labels = as.character(obs), pos = 3, cex = 1.2)
dev.off()

# Fig 3.a
pdf()
par(pty = "s")
plot(oeta079, prob079, xlim = c(-25,10), ylim = c(0,1.05), type = "l", xlab = "Linear predictor", ylab = "Probability", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
plot(Response ~ linear079, data = pred, subset = (Response == 0), ylim = c(0,1.05), xlim = c(-25,10), pch = 1, axes = FALSE, xlab = "", ylab = "", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
plot(Response ~ linear079, data = pred, subset = (Response == 1), ylim = c(0,1.05), xlim = c(-25,10), pch = 16, axes = FALSE, xlab = "", ylab = "", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
text(eta079[c(4,24)], finney$Response[c(4,24)], as.character(c(4,24)), pos = 3, cex = 1.2)
text(eta079[18], finney$Response[18], as.character(18), pos = 1, cex = 1.2)
dev.off()

# Fig 3.b
qq079 <- envelopeq.bin(fm079, reps = 500, conf = 0.95)
std <- qq079$std
pdf()
par(pty = "s")
qqnorm(std, ylim = ylim, main = "", ylab = "Studentized residuals", xlab = "Quantiles of standard normal", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
par(new = TRUE)
qqnorm(qq079$envelope[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
par(new = TRUE)
qqnorm(qq079$envelope[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
z <- qnorm(ppoints(n))[order(order(std))]
# order(std) reveals 'locations' of observations 4 and 18
text(sort(z)[c(38,39)], sort(std)[c(38,39)], labels = as.character(c(18,4)), pos = 3, cex = 1.2)
dev.off()
