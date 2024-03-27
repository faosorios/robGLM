## ID: qopt.R, last updated 2024-03-27, F.Osorio

qopt <- function(formula, data, qstart = 1, qmin = 0.75, by = 0.01) {
  # local functions
  norm2 <- function(x, y) {
    diff <- x - y
    sqrt(sum(diff^2))
  }

  # defs
  incr <- -by
  grid <- seq(from = qstart, to = qmin, by = incr)
  ngrid <- length(grid)
  db <- data

  # initialize progress bar
  cat(" Progress:\n")
  pb <- txtProgressBar(min = 0, max = ngrid, style = 3)

  qpar <- grid[1] # qstart, initially equals 1.0
  fm0 <- robLogistic(formula, data = db, qpar = qpar, control = glm.control(maxit = 100), x = TRUE)
  n <- nrow(fm0$x)
  p <- ncol(fm0$x)
  cf0 <- fm0$coef
  SE0 <- summary(fm0)$coefficients[,2]
  z0 <- cf0 / SE0

  # container
  QV <- rep(0, ngrid)

  # sampling and estimation
  for (i in 2:ngrid) {
    #
    qpar <- grid[i]
    fm1 <- robLogistic(formula, data = db, qpar = qpar)
    cf1 <- fm1$coef
    SE1 <- summary(fm1)$coefficients[,2]
    z1 <- cf1 / SE1

    QV[i] <- norm2(cf0, cf1)
    z0 <- z1
    cf0 <- cf1

    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")

  delta <- sqrt(sum(cf1^2))

  names(QV) <- as.character(grid)

  ok <- QV < 0.05 * delta
  if (all(ok))
    qpar <- 1.0
  else {
    qgrid <- grid[-1]
    ok <- ok[-1]
    qpar <- min(qgrid[ok == TRUE])
  }

  cf <- robLogistic(formula, data = db, qpar = qpar)$coef

  o <- list(coef = cf, qpar = qpar, QV = QV, ok = ok)
  o
}
