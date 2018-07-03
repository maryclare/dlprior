#' Code for paper is given at https://github.com/debdeeptamu/Dirichlet_Laplace/blob/master/DL.m
#' @export
sym.sq.root <- function(A) {
  A.eig <- eigen((A + t(A))/2)
  crossprod(t(A.eig$vectors), tcrossprod(diag(sqrt(ifelse(A.eig$values > 0, A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A)), A.eig$vectors))
}

#' @export
dl.beta <- function(XtX, Xty, sig.sq, psi, lambda) {

  if (is.matrix(XtX)) {

    p <- nrow(XtX)

    phiphit <- tcrossprod(lambda*sqrt(psi))
    A <- (XtX/sig.sq)*(phiphit) + diag(p)
    A.inv <- solve(A)*phiphit
    mean <- crossprod(A.inv, Xty/sig.sq)
    var <- A.inv

    return(crossprod(sym.sq.root(var), rnorm(p)) + mean)
  } else if (is.vector(XtX)) {
    p <- length(XtX)

    A <- XtX/sig.sq + 1/(lambda^2*psi)
    A.inv <- 1/(A)
    mean <- A.inv*Xty/sig.sq
    var <- A.inv

    return(sqrt(var)*rnorm(p) + mean)
  }
}

# Generates inverse gaussian rv's based on
# https://www.jstor.org/stable/2683801?origin=crossref, as described on Wiki page for
# inverse gaussian distribution https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution#cite_note-2
#' @export
rinvgauss <- function(p, mu, lambda) {
  v <- rnorm(p)
  y <- v^2
  x <- mu + mu^2*y/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*y + mu^2*y^2)
  z <- runif(p)
  return(ifelse(z <= mu/(mu + x), x, mu^2/x))
}

#' @export
dl.psi <- function(beta, lambda) {
  p <- length(beta)
  psi.tilde <- rinvgauss(p, mu = lambda/abs(beta), lambda = rep(1, p))
  return(1/psi.tilde)
}

#' @export
dl.tau <- function(beta, phi, a) {
  p <- length(beta)
  return(GIGrvg::rgig(1, lambda = p*a - p, psi = 1, chi = 2*sum(abs(beta/(phi)))))
}

#' @export
dl.lambda <- function(beta, a) {
  p <- length(beta)
  # return(GIGrvg::rgig(p, lambda = a - 1, psi = 1, chi = 2*abs(beta)))
  return(GIGrvg::rgig(p, lambda = 1 - a, chi = 1, psi = 2*abs(beta)))
}

# Uses slice sampling algorithm of Damien, Wakefield and Walker (1999)
#' @export
dl.phi <- function(beta, a, s) {
  # - Draw slice variable
  u <- runif(p, 0, exp(-1/(2*s)))
  # - Convert slice variable to lower bound on s
  lb <- 1/(2*log(1/u))
  # - Find the gamma cdf value that corresponds to this lower bound
  Flb <- pgamma(lb, shape = 1 - a, rate = abs(beta))
  # - Use inverse CDF method to draw gamma random variable
  uu <- pmin(runif(p, Flb, 1), 1-(1e-16))
  s <- qgamma(uu,shape = 1-a,rate = abs(beta))
  # - Convert back to t and phi
  t <- 1/s
  phi <- t/sum(t)
  phi[phi <= (1e-20)] <- (1e-20)
  return(list("phi" = phi, "s" = s))
}

#' @export
dl.sampler <- function(y, X, a, sig.sq, num.samp = 10000,
                       burn.in = 0, thin = 1,
                       print.iter = FALSE, lambdapar = TRUE) {

  p <- ncol(X)
  n <- nrow(X)

  XtX <- crossprod(X)
  if (max(abs(XtX[lower.tri(XtX, diag = FALSE)])) == 0) {
    diagX <- TRUE
    XtX <- diag(crossprod(X))
  }
  Xty <- crossprod(X, y)

  betas <- psis <- matrix(nrow = num.samp, ncol = p)
  if (!lambdapar) {
    phis <- matrix(nrow = num.samp, ncol = p)
    taus <- numeric(num.samp)
  } else {
    lambdas <- matrix(nrow = num.samp, ncol = p)
  }

  # Starting values
  psi <- rep(1, p)
  if (!lambdapar) {
    t <- s <- rep(1, p)
    phi <- t/sum(t)
    tau <- 1
    lambda <- phi*tau
  } else {
    lambda <- rep(1, p)
  }


  for (i in 1:(num.samp*thin + burn.in)) {

    if (print.iter) {cat("i = ", i, "\n")}

    beta <- dl.beta(XtX = XtX, Xty = Xty, sig.sq = sig.sq, psi = psi,
                    lambda = lambda)
    psi <- dl.psi(beta = beta, lambda = lambda)

    if (!lambdapar) {
      tau <- dl.tau(beta = beta, phi = phi, a = a)
      samp.phi <- dl.phi(beta = beta, a = a, s = s)
      phi <- samp.phi$phi
      s <- samp.phi$s
      lambda <- phi*tau
    } else {
      lambda <- dl.lambda(beta = beta, a = a)
    }

    if (i > burn.in & (i - burn.in)%%thin == 0) {
      betas[(i - burn.in)/thin, ] <- beta
      psis[(i - burn.in)/thin, ] <- psi
      if (!lambdapar) {
        phis[(i - burn.in)/thin, ] <- phi
        taus[(i - burn.in)/thin] <- tau
      } else {
        lambdas[(i - burn.in)/thin, ] <- lambda
      }
    }
  }

  if (!lambdapar) {
    res <- list("betas" = betas, "psis" = psis, "phis" = phis, "taus" = taus)
  } else {
    res <- list("betas" = betas, "psis" = psis, "lambdas" = lambdas)
  }

  return(res)

}
