####Stick Breaking Process Generator

#' @title Random Stick Breaking Process Generator
#' @description Random generator for stick breaking process, and posterior stick breaking process.
#'
#' @param alpha Concentration parameter in \code{R+}.
#' @param n Sample size (integer).
#' @param prior_probs a vector of prior probabilities.
#' @return \code{rsb} returns a list, a vector of sampled probabilities
#' and a vector of probabilities at each draw
#' @return \code{rsb_post} returns a list, a vector of sampled probabilities
#' and a vector of probability at each draw
#' @examples
#' samp<-rsb(n=5,alpha=5)
#' plot(1:5,samp$probs,type = "h")
#' probs <- c(.2,.2,.6)
#' samp<-rsb_post(alpha = 5, prior_probs = probs)
#' plot(1:3, probs,type = "h", col = alpha("red",.4), lwd = 5, xlim =c(0.5,3.5))
#' lines(1:3+.1, samp$probs,type = "h", col = alpha("blue", .4), lwd = 5)
#' samp$probs
#'
#' samp_size<-20
#' probs<-rsb(alpha = .2,n = samp_size)$probs
#' samp<-rsb_post(alpha = .5,prior_probs = probs)
#' plot(1:samp_size, probs,type = "h", col = alpha("red",.4), lwd = 5, xlim =c(0.5,samp_size+.5))
#' lines(1:samp_size+.1, samp$probs,type = "h", col = alpha("blue", .4), lwd = 5


rsb<-function(alpha, n){
  v<-rbeta(n-1,1,alpha)
  sb<-v*(c(1,cumprod(1-v))[-n])
  sb<-c(sb, 1-sum(sb))
  return(list(probs = sb,betas =  v))
}

###Stick Breaking Posterior Distribution

rsb_post <- function(alpha,prior_probs){
  p<-V<-numeric()
  n<-length(prior_probs)
  V[1]<-rbeta(1,(1 + prior_probs[1]),(alpha+sum(prior_probs[2:n])))
  p[1]<-V[1]
  if(V[1]==1){V[1]<-0.9999999999}

  for(i in 2:(n-1)){
    V[i]<-rbeta(1, (prior_probs[i] +1), (sum(prior_probs[(i+1):n]) + alpha))
    if(V[i] <= 0){V[i] <- 1e-308}
    if(V[i]==1){V[i]<-0.9999999999}
    p[i]<-V[i]*prod((1-V[1:(i-1)]))
  }

  p[n]<-1-sum(p[1:(n-1)])
  if(p[n] <= 0){p[n] <- 1e-308}
  return(list(probs = p,betas = V))
}
