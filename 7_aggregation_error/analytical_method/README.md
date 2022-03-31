Supplementary Information for the AGU paper on upscaling
<http://www.essoar.org/doi/10.1002/essoar.10509113.1>

## R-code for the Hessian matrix

The R code below gives the predicted upscaling error *Δ̂* by calculating
each block’s Hessian matrix (H). Variance-covariance matrices (S) were
calculated for each block from the input values of its constituent
cells, as described in:

Van Oijen, M., Cameron, D., Levy, P. E., & Preston, R. (2017).
Correcting errors from spatial upscaling of nonlinear greenhouse gas
flux models. Environmental Modelling & Software, 94, 157–165.
<https://doi.org/10.1016/j.envsoft.2017.03.023>

The code can be used with any model *f* to calculate its Hessian around
an input vector *x* (typically the grid-cell mean). The code consists of
three functions. The first function calculates first-order derivatives
(slopes). It is called by the second function, which calculates
second-order derivatives (curvatures) and is itself called by the third
function to populate the Hessian matrix.

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

## R-code for *Δ̂*

Direct implementation in R of Eq. (4) in van Oijen et al (2017) for the
estimation of *Δ* cannot be done with R’s base functions as they do not
include a function for calculating the trace of a matrix. To avoid
having to install an R-package (such as ‘matrixcalc’) that does contain
a trace function, our code implements the following formula which for
symmetric matrices such as S is mathematically equivalent to Eq. (4):

$$ \\widehat{\\Delta} = -\\frac{1}{2}\\sum\_{i,j}S(i,j)H(i,j). $$

The code that implements this formula is:

    Delta_hat <- function( f, x, S, dx=rep(1.E-2,length(x)) ) {
      Delta_hat <- -0.5 * sum( S * Hf(f=f,x=x,dx=dx) )
      return(Delta_hat)
    }
