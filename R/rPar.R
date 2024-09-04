#' Generate variance-covariance matrix for multivariate functional data
#'
#'
#'
#' @param p the number of functional curves.
#' @param K the number of basis.
#' @param id.diag an argument that controls the sparsness of the matrices.
#' @param s1 an scale factor.
#' @param s2 a scale factor.
#'
#' @description
#' This function generates a variance-covariance matrices and the corresponding precision matrices from a multivariate functional data under partial separability assumption of the covariance operator
#'
#'
#' @export

rPar <- function(p, K, id.diag = 1, s1 = 3, s2 = -1.8) {
  # Argomenti
  # p: dimesione del vettore degli scores
  # K: numero di termini utilizzati per l'espansione multivariate di KL
  # id.diag =
  # s1, s2: fattore di scale

  #--- Setting Output
  Sgm <- Tht <- array(0, dim = c(p, p, K),
                      dimnames = list(paste0("X", seq_len(p)),
                                      paste0("X", seq_len(p)),
                                      paste0("K", seq_len(K))))
  for (k in seq_len(K)) {
    Tht[, , k] <- .rTht(p = p, id.diag = id.diag)
    sc <- s1 * k^s2
    Tht[, , k] <- Tht[, , k] / sc
    Sgm[, , k] <- solve(Tht[, , k])
  }
  return(list(Tht = Tht, Sgm = Sgm))
}
