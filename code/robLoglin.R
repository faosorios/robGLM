## ID: robLoglin.R, last updated 2023-02-27, F.Osorio

robLoglin <-
function(formula, data, subset, na.action, qpar = 1, control = glm.control(),
  model = TRUE, x = TRUE, y = TRUE, contrasts = NULL, ...)
{
  # support functions
  linkfun <- function(mu) log(mu)
  linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
  kernel <- function(eta) exp(eta)
  variance <- function(mu) mu

  # extract model matrix and response variable
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$qpar <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  # mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  ynames <- names(y)
  xnames <- dimnames(x)[[2]]
  dims <- dim(x)
  nobs <- dims[1]

  # additional control parameters
  eps <- 1e-2

  # checking arguments
  if (any(y < 0))
	stop("negative values not allowed")
  if (is.null(control))
    control <- glm.control()

  now <- proc.time()
  ## initial estimates
  z   <- log(y + 0.1) - 0.1 / (y + 0.1)
  wts <- y + 0.1
  fit <- lsfit(x, z, wt = wts, intercept = FALSE)
  oldcoef <- fit$coefficients
  eta <- c(x %*% oldcoef)

  # IRLS iterations
  iter <- 0
  repeat {
    mu   <- linkinv(eta)
    rob.wts <- exp((1 - qpar) * (y * eta - kernel(eta) - lgamma(y + 1)))
    k    <- exp(qpar * kernel(eta) - kernel(qpar * eta))
    wts  <- k * variance(mu)
    z    <- eta + rob.wts * (y - mu) / wts
    fit  <- lsfit(x, z, wt = wts, intercept = FALSE)
    coef <- fit$coef

    iter <- iter + 1
    diff <- coef - oldcoef
    conv <- sum(diff^2) / (sum(oldcoef^2) + eps)

    if (conv < control$epsilon) break
    if (iter >= control$maxit)  break

    eta <- c(x %*% coef)
    oldcoef <- coef
  }

  ## calculate df
  nulldf <- nobs - sum(wts == 0)
  rank   <- fit$qr$rank
  resdf  <- nobs - rank

  ## Lq-likelihood
  if (qpar != 1)
    LqLik <- sum(rob.wts - 1) / (1 - qpar)
  else
    LqLik <- sum(y * eta - kernel(eta) - lgamma(y + 1))

  ## AIC
  p <- rank # must be length(coef)
  aic <- -2 * LqLik + 2 * p / (2 - qpar)

  speed <- proc.time() - now
  ## creating the output object
  mu <- linkinv(eta)
  R <- qr.R(qr(sqrt(wts) * x))
  ## unscaled covariance
  #R <- solve(R)
  cov.unscaled <- chol2inv(R) # tcrossprod(R)
  out <- list(call = Call,
              dims = dims,
              qpar = qpar,
              coefficients = qpar * coef,
              fitted.values = mu,
              linear.predictor = qpar * eta,
              residuals = y - mu,
              numIter = iter,
              control = control,
              weights = wts,
              rob.weights = rob.wts,
              rank = rank,
              df.null = nulldf,
              df.residual = resdf,
              LqLik = LqLik,
              AIC = aic,
              R = R,
              cov.unscaled = cov.unscaled,
              cov.trace = sum(diag(cov.unscaled)),
              speed = speed)
  out$terms <- Terms
  if (model)
    out$model <- mf
  if (ret.y)
    out$y <- y
  if (ret.x)
    out$x <- x
  class(out) <- "loglin"
  out
}

print.loglin <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  dput(x$call, control = NULL)
  cat("Converged in", x$numIter, "iterations\n\n")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else
    cat("No coefficients\n\n")
  cat("\nDistortion parameter:", x$qpar)
  cat("\nDegrees of Freedom:", x$df.null, "total;", x$df.residual, "residual\n")
  invisible(x)
}

summary.loglin <-
function (object, ...)
{
  z <- object
  p <- z$dims[2]
  #R <- solve(z$R)
  cov.unscaled <- z$cov.unscaled #R %*% t(R)
  qpar <- z$qpar
  se <- sqrt(diag(cov.unscaled) / (2 - qpar))
  est <- z$coefficients
  zval <- est / se
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$rank <- z$rank
  ans$qpar <- qpar
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients),
        c("Estimate", "Std.Error", "Z value", "p-value"))
  class(ans) <- "summary.loglin"
  ans
}

print.summary.loglin <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  dput(x$call, control = NULL)
  resid <- x$residuals
  nobs <- x$dims[1]
  rdf <- nobs - x$rank
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nDistortion parameter:", x$qpar)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
  invisible(x)
}
