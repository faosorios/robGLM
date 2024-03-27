## ID: envelope.R, last updated 2022-11-15, F.Osorio

## envelope.bin has been adapted from original sources by
## Gilberto A. Paula, https://www.ime.usp.br/~giapaula/envel_bino

envelope.bin <- function(object, reps = 50, conf = 0.95) {
    ## simulated envelope
    x <- model.matrix(object)
    n <- nrow(x)
    p <- ncol(x)
    w <- diag(object$weights)
    H <- solve(t(x) %*% w %*% x)
    H <- sqrt(w) %*% x %*% H %*% t(x) %*% sqrt(w)
    h <- diag(H)
    td <- resid(fm, type = "deviance") / sqrt(1 - h)

    #
    conf <- 1 - conf
    elims <- matrix(0, nrow = n, ncol = reps)
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    for(i in 1:reps) {
      dif <- runif(n) - fitted(object)
      dif[dif >= 0 ] <- 0
      dif[dif < 0] <- 1
      nresp <- dif
      fit <- glm(nresp ~ x, family = binomial)
      w <- diag(fit$weights)
      H <- solve(t(x) %*% w %*% x)
      H <- sqrt(w) %*% x %*% H %*% t(x) %*% sqrt(w)
      h <- diag(H)
      elims[,i] <- sort(resid(fit, type = "pearson") / sqrt(1-h))
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band

    ylim <- range(td, band)
    qqnorm(td, ylim = ylim, main = "", ylab = "Studentized residuals", xlab = "Quantiles of standard normal", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

    invisible(list(std.resid = td, envelope = band, ylim = ylim))
}

envelopeq.bin <- function(object, reps = 50, conf = 0.95) {
    ## simulated envelope
    x <- object$x
    n <- nrow(x)
    p <- ncol(x)
    WJ <- diag(object$weights)
    xx <- object$cov.unscaled
    M <- x %*% xx %*% crossprod(x, WJ)
    m <- diag(M)
    ti <- object$rob.weights * object$residuals / sqrt(object$weights)
    ti <- sqrt(2 - object$qpar) * ti / sqrt(1 - m)
    yfit <- object$fitted.values
    qpar <- object$qpar

    #
    conf <- 1 - conf
    elims <- matrix(0, nrow = n, ncol = reps)
    ok <- rep(0, reps)
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    for(i in 1:reps) {
      try <- 0
      repeat {
        # simulation and fitted model
        dif <- runif(n) - yfit
        dif[dif >= 0 ] <- 0
        dif[dif < 0] <- 1
        nresp <- dif
        fit <- tryCatch(robLogistic(nresp ~ -1 + x, qpar = qpar),
                        warning = function(w) list(),
                        error = function(e) list())
        try <- try + 1
        if (!is.null(fit$coef))
          break
      }
      ok[i] <- try
      WJ <- diag(fit$weights)
      xx <- fit$cov.unscaled
      M <- x %*% xx %*% crossprod(x, WJ)
      m <- diag(M)
      res <- fit$rob.weights * fit$residuals / sqrt(fit$weights)
      res <- sqrt(2 - qpar) * res / sqrt(1 - m)
      elims[,i] <- sort(res)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band

    ylim <- range(ti, band)
    qqnorm(ti, ylim = ylim, main = "", ylab = "Studentized residuals", xlab = "Quantiles of standard normal", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2, lwd = 2)
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")

    invisible(list(std.resid = ti, envelope = band, ylim = ylim, ok = ok))
}
