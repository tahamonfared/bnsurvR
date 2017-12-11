
#################################Multi Variate Gamma

#### Random Multivariate Gamma Generator

#'@title Multivariate Gamma Distribution
#'@description Distribution function and random generator for the multivariate gamma distribution
#'with parameters \code{mu} for location, \code{alpha} for shape, and \code{beta} for scale. \code{alpha},
#'\code{x}, and \code{beta} should be in \code{R+}. \code{mu} resides in \code{R}.
#'
#'@param mu vector of location parameter.
#'@param alpha vector of shape parameter.
#'@param beta scale parameter.
#'@param x matrix of observations from a multivariate gamma distribution, columnar.
#'@param samp_size number of observations.
#'@param log logical; if \code{TRUE}, probabilities are returned as log(p).
#'
#'@retun \code{rmvgamma} returns a matrix of generated observations from multivariate
#'gamma distribution.
#'
#'\code{dmvgamma} returns a vector of distribution function of a given observation matrix,
#'compares the observations to the given
#'multivariate gamma distribution specified by the arguments in the function.
#'
#'@examples
#'samp<-rmvgamma(mu = c(-1,2), alpha = c(.1,2), beta = 3, samp_size = 10)
#'samp
#'dmvgamma(x = samp, mu = c(-1,2), alpha = c(1,1), beta = 3, log = FALSE)


rmvgamma<-function(mu,alpha,beta,samp_size){
  if (length(mu)!=length(alpha)){
    print(paste("Length of shape and location parameters",
                "\nshould be the same length as the dimensions",
                "\nof your required multivariate gamma distribution."))
    stop()
  }
  if(length(beta)!=1){
    print("This multivariate gamma has common scale parameter,
          \nso to get the same shape marginals.")
    stop()
  }
  mt<-cbind(mu,alpha,rep(beta,length(mu)))
  gamma_func<-function(vec){
    mu_i<-vec[1]
    alpha_i<-vec[2]
    beta_i<-vec[3]
    rgamma(n = samp_size,shape = alpha_i,rate = beta_i)+mu_i
  }
  v<-apply(X = mt,1,gamma_func)
  if(length(mu)==1){
    Z<-v
  }else{
    Z<-t(apply(v,1,cumsum))
  }
  return(Z)
}

###Multivariate gamma distribution

dmvgamma<-function(x,mu,alpha,beta,log=FALSE){
  ###Accepts a matrix X, n*m. A dataset of n observations from m dimensional gamma distribution
  ###Produces the probability of each observation in m dimensional gamma distribution with
  ###location parameter mu, an m length vector of alphas, a sigle beta. It can also create the
  ###log-probability if the log=TRUE

  if (length(mu)!=length(alpha)){
    print("Length of shape and location parameters
          should be the same length as the dimensions
          of your observed multivariate gamma distribution.")
    stop()
  }
  if(length(beta)!=1){
    print("This multivariate gamma has common scale parameter,
          so to get the same shape marginals.")
    stop()
  }

  calc_prob<-function(x){

    first<-x[1]
    x<-c(first,diff(x))
    if(any((x-mu)<=0)){
      return(0)
    }
    logd<-sum(log(x-mu)*(alpha-1))-
      (sum(x)-sum(mu))/beta-
      log(beta)*sum(alpha)-sum(lgamma(alpha))
    if (log == TRUE){return(logd)}
    else{return(exp(logd))}
  }
  x<-matrix(x,ncol=length(mu))
  res<-apply(x,1,calc_prob)
  return(res)
}

