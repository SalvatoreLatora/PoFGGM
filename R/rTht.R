#' Generate precision matrices for multivariate functional data

.rTht <- function(p, id.diag, det.min = 0.5) {
   # argomenti
   # p = dimensione della matrice
   # id.diag =
   # det.min = valore del determinate della matrice di precisione generata
   if (!(id.diag %in% 0:(p - 1)))
     stop("prova")
   rho <- optimise(\(rho) {
     outer(1:p, 1:p, function(i, j) {
       d <- abs(i - j)
       ifelse(d <= id.diag, rho^abs(i - j), 0)
     }) -> Tht
     (det(Tht) - det.min)^2
   }, c(0, 1))$min
   outer(1:p, 1:p, function(i, j) {
     d <- abs(i - j)
     ifelse(d <= id.diag, rho^abs(i - j), 0)
   }) -> Tht
   return(Tht)
 }
