#' Generate Partially Observed Multivariate Functional Data Under Partial Separability
#'
#' @description
#' Generate Partially Observed Multivariate Functional Data Under Partial Separability assumption of the covariance operator.
#'
#'
#'
#' @param n the sample size.
#' @param a the lower bound of the full data domain. By default \code{a = 0}.
#' @param b the upper bound of the full data domain. By default \code{b = 1}.
#' @param d the number of full observed temporal points. By default \code{d = 60}.
#' @param w the number of missing temporal points. By default \code{w = 20}.
#' @param p the number of functional curves.
#' @param basis the type of basis that can be used. It is possible to choice \sQuote{\code{fourier}} or \sQuote{\code{bspline}} basis. By defaults \code{basis = "fourier"}.
#' @param nbasis the number of basis
#'
#'
#' @details
#' This function generate a set of partially observed multivariate functional data under partial separability assumption of the covariance operator. Other argument of the function \link{rPar} can be specified to generate a more or less sparse variance-covariance matrix.
#'
#'
#'
#' @return This function returns a list with the following objects:
#' \item{X}{ an array of the generated full observed multivariate functional data. By row the temporal points, by column the units, by the third dimension the functional curves.}
#' \item{Xpo}{an array of the generated partially observed multivariate functional data. By row the temporal points, by columns the units, by the third dimension the functional curves.}
#' \item{Xi}{an array of the true scores related to the full observed data. By row the units, by columns the functional curves, by the third dimension the basis.}
#' \item{U_mat}{ an array of the observed temporal grid. By row the temporal points, by columns the units, by the third dimension the functional covariates.}
#' \item{Tht}{ the true precision matrices under partial separability.}
#' \item{Sgm}{ the true var-cov matrices under partial separability.}
#'
#' @examples
#' set.seed(42)
#' rmfdpo(n = 100, p = 5)
#'
#' @export


rmfdpo <- function(n, a = 0, b = 1,  d = 60, w = 20, p, Sgm = NULL, basis = c("fourier", "bspline"), nbasis = 5, ...){

  if(missing(n)) stop("the argument n was not specified")
  if(missing(p)) stop("the argument p was not specified")
  if(p < 2) stop("the argument p must be greater than or equal to 2")
  if(is.null(Sgm)) Sgm <- rPar(p = p, K = nbasis, ...)$Sgm

  basis <- match.arg(basis)

  tp <- seq(a, b, length = d)
  Phi <- switch(basis,
    fourier = {
      phi.f <- create.fourier.basis(rangeval = c(a, b), nbasis = nbasis)
      Phi <- eval.basis(evalarg = tp, basisobj = phi.f)
      rownames(Phi) <- paste("tp = ", 1L:d)
      Phi
    },
    bspline = {
      phi.f <- create.fourier.basis(rangeval = c(a, b), nbasis = nbasis)
      Phi <- eval.basis(evalarg = tp, basisobj = phi.f)
      rownames(Phi) <- paste("tp =", 1L:d)
      Phi
    }
  )
  Xi <- array(0, dim = c(n, p, nbasis),
              dimnames = list(1:n,
                              paste0("Xi", seq_len(p)),
                              paste0("k", seq_len(nbasis))))
  X <- array(0, dim = c(d, n, p))
  U_mat <- array(NA, dim = c(d - w, n, p))
  Mask <- array(T, dim = c(d, n, p))
  for (k in seq_len(nbasis)) {
    Xi[, , k] <- MASS::mvrnorm(n = n,
                               mu = rep(0, p),
                               Sigma = Sgm[, , k])
  }
  for (j in seq_len(p)){
    X[, , j] <- Phi %*% t(Xi[, j, ])
  }
  Xpo <- X
  for (j in seq_len(p)) {
    ii <- sample(d - w + 1, 1)
    Mask[seq(from = ii, length = w), , j] <- FALSE
  }
  Xpo[!Mask] <- NA
  for(j in seq_len(p)) {
    for (i in seq_len(n)) {
      U_mat[,i,j]<- tp[Mask[, i, j]]
    }
  }
  out <- list("X" = X, "Xpo" = Xpo, "Xi" = Xi, "U_mat" = U_mat, "Sgm" = Sgm)
  class(out) <- 'pomfd'
  return(out)
}
