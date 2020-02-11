""" Initializing custom Bingham functions from R source code """

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

directional = importr('Directional')

""" Initialize default rvmf, rbingham and custom rbingham functions for 
    input data and/or user-defined parameters """
def Bing_def():
    string_rbing = """
    rbingham <- function(n, A) {
       p <- ncol(A)  ## dimensionality of A
       eig <- eigen(A)
       lam <- eig$values  ## eigenvalues
       V <- eig$vectors  ## eigenvectors
       lam <- lam - lam[p]
       lam <- lam[-p]
       ### f.rbing part
       lam <- sort(lam, decreasing = TRUE)  ## sort the eigenvalues in desceding order
       nsamp <- 0
       X <- NULL
       lam.full <- c(lam, 0)
       qa <- length(lam.full)
       mu <- numeric(qa)
       sigacginv <- 1 + 2 * lam.full
       SigACG <- sqrt( 1 / ( 1 + 2 * lam.full ) )
       Ntry <- 0
    
       while (nsamp < n) {
         x.samp <- FALSE
         while ( !x.samp ) {
           yp <- rnorm(qa, mu, SigACG)
           y <- yp / sqrt( sum( yp^2 ) )
           lratio <-  - sum( y^2 * lam.full ) - qa/2 * log(qa) + 0.5 * (qa - 1) + qa/2 * log( sum(y^2 * sigacginv ) )
           if ( log(runif(1) ) < lratio) {
             X <- c(X, y)
             x.samp <- TRUE
             nsamp <- nsamp + 1
           }
           Ntry <- Ntry + 1
         }
       }
    
       x <- matrix(X, byrow = TRUE, ncol = qa)
       ## the avtry is the estimate of the M in rejection sampling
       ## 1/M is the probability of acceptance
       ## the x contains the simulated values
       tcrossprod(x, V) ## simulated data
     }
    """
    powerpack = SignatureTranslatedAnonymousPackage(string_rbing, "powerpack")
    
    return powerpack

def Bing_cust(lam1,lam2,lam3):
    string_rbing1 = """
    rbingham <- function(n, A) {
       p <- ncol(A)  ## dimensionality of A
       eig <- eigen(A)
       V <- eig$vectors  ## eigenvectors
       lam <- c(%f,%f,%f)
       lam <- lam - lam[p]
       lam <- lam[-p]
       ### f.rbing part
       lam <- sort(lam, decreasing = TRUE)  ## sort the eigenvalues in desceding order
       nsamp <- 0
       X <- NULL
       lam.full <- c(lam, 0)
       qa <- length(lam.full)
       mu <- numeric(qa)
       sigacginv <- 1 + 2 * lam.full
       SigACG <- sqrt( 1 / ( 1 + 2 * lam.full ) )
       Ntry <- 0
    
       while (nsamp < n) {
         x.samp <- FALSE
         while ( !x.samp ) {
           yp <- rnorm(qa, mu, SigACG)
           y <- yp / sqrt( sum( yp^2 ) )
           lratio <-  - sum( y^2 * lam.full ) - qa/2 * log(qa) + 0.5 * (qa - 1) + qa/2 * log( sum(y^2 * sigacginv ) )
           if ( log(runif(1) ) < lratio) {
             X <- c(X, y)
             x.samp <- TRUE
             nsamp <- nsamp + 1
           }
           Ntry <- Ntry + 1
         }
       }
    
       x <- matrix(X, byrow = TRUE, ncol = qa)
       ## the avtry is the estimate of the M in rejection sampling
       ## 1/M is the probability of acceptance
       ## the x contains the simulated values
       tcrossprod(x, V) ## simulated data
     }
    """ % (lam1,lam2,lam3) # 200,0.05
    
    powerpack1 = SignatureTranslatedAnonymousPackage(string_rbing1, "powerpack")
    
    return powerpack1
