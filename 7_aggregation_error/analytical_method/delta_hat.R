
df_dxi <- function( f, x, i=1, dxi=1.E-2 ) {
  dimf   <- length(x)
  dx_i   <- rep(0,dimf) ; dx_i[i] <- dxi
  xiplus <- x + dx_i
  df_xi  <-( f(xiplus) - f(x) ) / dxi
  return(df_xi)
}

d2f_dxixj <- function( f, x, i=1, j=2, dxi=1.E-2, dxj=1.E-2 ) {
  dimf      <- length(x)
  dx_j      <- rep(0,dimf) ; dx_j[j] <- dxj
  xjplus    <- x + dx_j
  d2f_dxixj <- ( df_dxi(f,xjplus,i,dxi) - df_dxi(f,x,i,dxi) ) / dxj
  return(d2f_dxixj)
}

Hf <- function( f, x, dx=rep(1.E-2,length(x)) ) {
  dimf <- length(x)
  H    <- matrix(0,dimf,dimf)
  for(i in 1:dimf) {
    for(j in 1:dimf) {
      dxi<- dx[i] ; dxj <- dx[j]
      H[i,j] <- d2f_dxixj( f, x=x, i=i, j=j, dxi=dxi, dxj=dxj )
    }
  }
  return(H)
}

Delta_hat <- function( f, x, S, dx=rep(1.E-2,length(x)) ) {
  Delta_hat <- -0.5 * sum( S * Hf(f=f,x=x,dx=dx) )
  return(Delta_hat)
}
