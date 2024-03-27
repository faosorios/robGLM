## Id: simul_Poi.R, last updated 2023-05-25, F.Osorio

simul.Poi <- function(Nsize = 1000, nobs = 100, eps = 0.05, nu = 5, trace = TRUE, msg = NULL)
{
  # setting covariates
  set.seed(43) # 1st prime number after 42
  cf  <- c(1,1,1)
  x1  <- runif(nobs)
  x2  <- runif(nobs)
  x3  <- runif(nobs)
  eta <- cf[1] * x1 + cf[2] * x2 + cf[3] * x3
  X <- cbind(x1, x2, x3)

  simul.resp <- function(nobs, eta, eps, nu) {
    # model building
    mu <- exp(eta)
    y  <- rpois(nobs, mu)
    if (eps > 0) { # if contamination is required
      obs <- sample(1:nobs, size = eps * nobs)
      y[obs] <- nu * y[obs]
    }
    y
  }

  if (is.null(msg))
    msg <- "Progress"

  # results containers
  cf100 <- matrix(0, nrow = Nsize, ncol = 3)
  cf099 <- matrix(0, nrow = Nsize, ncol = 3)
  cf097 <- matrix(0, nrow = Nsize, ncol = 3)
  cf095 <- matrix(0, nrow = Nsize, ncol = 3)
  cf093 <- matrix(0, nrow = Nsize, ncol = 3)
  cf091 <- matrix(0, nrow = Nsize, ncol = 3)
  cf089 <- matrix(0, nrow = Nsize, ncol = 3)
  cf087 <- matrix(0, nrow = Nsize, ncol = 3)
  cf085 <- matrix(0, nrow = Nsize, ncol = 3)
  cfrob <- matrix(0, nrow = Nsize, ncol = 3)
  qopt  <- rep(0, Nsize)
  now <- proc.time()

  if (trace) {
    cat(" ", paste(msg, ":", sep = ""), "\n")
    pb <- txtProgressBar(min = 0, max = Nsize, style = 3)
  }

  # Monte Carlo iterations
  for (i in 1:Nsize) {
    set.seed(97 + i) # 97

    # simulation
    y  <- simul.resp(nobs, eta, eps, nu)
    Xy <- cbind(X, y)
    Xy <- as.data.frame(Xy)

    # fitting model
    f100 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 1.00, control = glm.control(maxit = 100))
    f099 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.99, control = glm.control(maxit = 100))
    f097 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.97, control = glm.control(maxit = 100))
    f095 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.95, control = glm.control(maxit = 100))
    f093 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.93, control = glm.control(maxit = 100))
    f091 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.91, control = glm.control(maxit = 100))
    f089 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.89, control = glm.control(maxit = 100))
    f087 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.87, control = glm.control(maxit = 100))
    f085 <- robLoglin(y ~ 0 + ., data = Xy, qpar = 0.85, control = glm.control(maxit = 100))
    fQle <- glmrob(y ~ 0 + ., data = Xy, family = poisson, method = "Mqle")

    cf100[i,] <- f100$coef
    cf099[i,] <- f099$coef
    cf097[i,] <- f097$coef
    cf095[i,] <- f095$coef
    cf093[i,] <- f093$coef
    cf091[i,] <- f091$coef
    cf089[i,] <- f089$coef
    cf087[i,] <- f087$coef
    cf085[i,] <- f085$coef
    cfrob[i,] <- fQle$coef

    # update progress bar
    if (trace)
      setTxtProgressBar(pb, i)
  }
  if (trace)
    close(pb)
  
  # summaries
  z100 <- colMeans(sweep(cf100, 2, cf))
  z099 <- colMeans(sweep(cf099, 2, cf))
  z097 <- colMeans(sweep(cf097, 2, cf))
  z095 <- colMeans(sweep(cf095, 2, cf))
  z093 <- colMeans(sweep(cf093, 2, cf))
  z091 <- colMeans(sweep(cf091, 2, cf))
  z089 <- colMeans(sweep(cf089, 2, cf))
  z087 <- colMeans(sweep(cf087, 2, cf))
  z085 <- colMeans(sweep(cf085, 2, cf))
  zrob <- colMeans(sweep(cfrob, 2, cf))
  bias <- rep(0, 10)
  bias[1]  <- sqrt(sum(z100^2))
  bias[2]  <- sqrt(sum(z099^2))
  bias[3]  <- sqrt(sum(z097^2))
  bias[4]  <- sqrt(sum(z095^2))
  bias[5]  <- sqrt(sum(z093^2))
  bias[6]  <- sqrt(sum(z091^2))
  bias[7]  <- sqrt(sum(z089^2))
  bias[8]  <- sqrt(sum(z087^2))
  bias[9]  <- sqrt(sum(z085^2))
  bias[10] <- sqrt(sum(zrob^2))
  names(bias) <- c("ML","099","097","095","093","091","089","087","085","Rob")

  iqr <- rep(0, 10)
  iqr[1]  <- mean(apply(cf100, 2, IQR))
  iqr[2]  <- mean(apply(cf099, 2, IQR))
  iqr[3]  <- mean(apply(cf097, 2, IQR))
  iqr[4]  <- mean(apply(cf095, 2, IQR))
  iqr[5]  <- mean(apply(cf093, 2, IQR))
  iqr[6]  <- mean(apply(cf091, 2, IQR))
  iqr[7]  <- mean(apply(cf089, 2, IQR))
  iqr[8]  <- mean(apply(cf087, 2, IQR))
  iqr[9]  <- mean(apply(cf085, 2, IQR))
  iqr[10] <- mean(apply(cfrob, 2, IQR))
  names(iqr) <- c("ML","099","097","095","093","091","089","087","085","Rob")
  
  speed <- proc.time() - now

  out <- list(cf100 = cf100, cf099 = cf099, cf097 = cf097, cf095 = cf095, cf093 = cf093, 
              cf091 = cf091, cf089 = cf089, cf087 = cf087, cf085 = cf085, cfrob = cfrob, 
              bias = bias, iqr = iqr, speed = speed)
  out
}

