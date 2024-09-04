plot.pomfd <- function(x, ...){
  if(missing(x)) stop("argument x was not specified")
  p <- dim(x$X)[3]
  if(p %% 2 == 0){
    par(mfrow(p/2, P/2))
  }
  else {
    par(mfrow = c(p/2-0.5, p/2+0.5))
  }
  for(i in seq_len(p)){
    matplot(x$X[, , i], type = "l", col = 'grey50', ...)
    matlines(x$Xpo[, , i], col = "red")


  }
}
