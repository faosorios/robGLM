# load Leukemia dataset
load("../data/leukemia.rda")

# source .R files and load required packages
source("../code/robLogistic.R")
source("../code/qopt.R")
source("../code/envelope.R")
library("robust")
library("robustbase")

# fit logistic regression using ML 
f0 <- glm(y ~ WBC + AG, data = leukemia, family = binomial)

# fitted model removing observation 15
f1 <- glm(y ~ WBC + AG, data = leukemia, family = binomial, subset = -15)

# fit logistic regression with q = 1
fm100 <- robLogistic(y ~ WBC + AG, data = leukemia, qpar = 1.0)

# fit logistic regression with optimal q = 0.81
fm081 <- robLogistic(y ~ WBC + AG, data = leukemia, qpar = 0.81)

# fit logistic regression using robust procedures
fmCR <- glmrob(y ~ WBC + AG, data = leukemia, family = binomial, method = "Mqle")
fmCUBIF <- glmRob(y ~ WBC + AG, data = leukemia, family = binomial, method = "cubif")
fmBY <- glmrob(y ~ WBC + AG, data = leukemia, family = binomial, method = "BY")

# Table 2
CF <- matrix(0, nrow = 6, ncol = 6)
s0 <- summary(f0)
s1 <- summary(f1)
s081 <- summary(f081)
sCR <- summary(fmCR)
sCUBIF <- summary(fmCUBIF)
sBY <- summary(fmBY)

CF[1,] <- as.vector(t(s0$coef[,1:2]))
CF[2,] <- as.vector(t(s1$coef[,1:2]))
CF[3,] <- as.vector(t(s081$coef[,1:2]))
CF[4,] <- as.vector(t(sCR$coef[,1:2]))
CF[5,] <- as.vector(t(sCUBIF$coef[,1:2]))
CF[6,] <- as.vector(t(sBY$coef[,1:2]))
CF[,3:4] <- CF[,3:4] * 10^4
colnames(CF) <- c("Int","SE","WBC","SE","AG","SE")
rownames(CF) <- c("ML","ML*","MLq","CR","CUBIF","BY")

#CF
#             Int        SE        WBC        SE       AG        SE
#ML    -1.3073446 0.8144843 -0.3177170 0.1863111 2.261066 0.9522158
#ML*    0.2119170 1.0830472 -2.3544630 1.3540293 2.558064 1.2341075
#MLq   -0.9359840 0.8978418 -0.5475300 0.3667954 2.111350 1.0235823
#CR     0.1710991 1.0802041 -2.0417790 1.3613066 2.484917 1.2472440
#CUBIF -1.0450338 0.8503520 -0.5274136 0.3033654 2.220183 0.9795610
#BY    -1.3364274 1.2563473 -0.3247849 0.4834160 2.311365 1.0034179

# Figure 4.a
wts <- fm081$rob.weights
obs <- (1:33)[wts < .8]
obs <- obs[-4]
par(pty = "s")
plot(wts, ylim = c(0,1), ylab = "Estimated weights", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
abline(h = 1, lwd = 2, lty = 2, col = "red")
text(obs, wts[obs], labels = as.character(obs), pos = 3, cex = 1.2)
text(18, wts[18], labels = as.character(18), pos = 1, cex = 1.2)

# Figure 4.b
wr <- fmCR$w.r
obs <- (1:33)[wr < 1]
obs <- obs[-4]
par(pty = "s")
plot(wr, ylim = c(0,1), ylab = "Robust weights", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
abline(h = 1, lwd = 2, lty = 2, col = "red")
text(obs, wr[obs], labels = as.character(obs), pos = 3, cex = 1.2)
text(18, wr[18], labels = as.character(18), pos = 1, cex = 1.2)

# Figure 5.a
CD <- cooks.distance(f0)
obs <- (1:33)[wts < .8]
obs <- obs[-4]
par(pty = "s")
plot(CD, ylim = c(0,1.2), ylab = "Cook's distances", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
text(obs, CD[obs], labels = as.character(obs), pos = 3, cex = 1.2)
text(18, CD[18], labels = as.character(18), pos = 1, cex = 1.2)

# Figure 5.b
par(pty = "s")
env <- envelopeq.bin(fm081, reps = 5000)
