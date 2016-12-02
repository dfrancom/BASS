#' @title BASS Plot Diagnostics
#'
#' @description Generate diagnostic plots for BASS model fit.
#' @param mod a \code{bass} object.
#' @param quants quantiles for intervals, if desired.  NULL if not desired.
#' @param ... graphical parameters.
#' @details The first two plots are trace plots for diagnosing convergence.  The third plot is posterior predicted vs observed, with intervals for predictions.  The fourth plot is a histogram of the residuals (of the posterior mean model), with a red curve showing the assumed Normal density (using posterior mean variance).
#' @keywords BMARS
#' @export
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # See examples in bass documentation.
#'
plot.bass<-function(mod,quants=c(.025,.975),...){
  op<-par(no.readonly=T)
  par(mfrow=c(2,2))
  plot(mod$nbasis,type='l',ylab='number of basis functions',xlab='MCMC iteration (post-burn)')
  plot(mod$s2,type='l',ylab='error variance',xlab='MCMC iteration (post-burn)')
  margin<-2
  if(mod$func)
    margin<-2:3
  s<-sqrt(mod$s2)
  if(!is.null(quants)){
    qq1<-apply(mod$yhat+qnorm(quants[2])*s,margin,quantile,probs=quants[2])
    qq2<-apply(mod$yhat+qnorm(quants[1])*s,margin,quantile,probs=quants[1])
    ylim=range(c(qq1,qq2))
    ylab='interval'
  } else{
    ylim=range(mod$yhat.mean)
    ylab='mean'
  }
  plot(mod$y,mod$yhat.mean,ylim=ylim,ylab=paste('posterior predictive',ylab),xlab='observed',main='Training Fit',type='n',...)
  if(!is.null(quants))
    segments(mod$y,qq1,mod$y,qq2,col='lightgrey')
  points(mod$y,mod$yhat.mean)
  abline(a=0,b=1,col=2)

  hist(mod$y-mod$yhat.mean,freq=F,main='Posterior mean residuals',xlab='residuals')
  curve(dnorm(x,sd=mean(s)),col=2,add=T)
  par(op)
}


#' @title Plot BASS sensitivity indices
#'
#' @description Generate plots for sensitivity analysis of BASS.
#' @param sens a \code{bassSob} object, returned from \code{sobol}.
#' @param ... graphical parameters.
#' @details If \code{func.var} in the call to \code{sobol} was \code{NULL}, this returns boxplots of sensitivity indices and total sensitivity indices.  If there were functional variables, they are labeled with integers falling after labels for the regular inputs.  Thus, if I fit a model with 4 regular inputs and 2 functional inputs, the functional inputs are labeled 5 and 6.  If \code{func.var} was not \code{NULL}, then posterior mean functional sensitivity indices are plotted, along with the functional partitioned variance.  Variables that are excluded did not explain any varaince.
#' @keywords BMARS
#' @export
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # See examples in bass documentation.
#'
plot.bassSob<-function(sens,...){
  op<-par(no.readonly=T)
  par(mfrow=c(1,2),xpd=T)
  if(sens$func){
    ord<-order(sens$xx)
    sens.mean<-apply(sens$S,2:3,mean)
    matplot(sens$xx[ord],t(apply(sens.mean,2,cumsum))[ord,],type='l',xlab='x',ylab='proportion variance',ylim=c(0,1),main='Sensitivity',...)
    lab.x<-apply(sens.mean,1,which.max)
    cs<-rbind(0,apply(sens.mean,2,cumsum))
    cs.diff<-apply(sens.mean,2,function(x) diff(cumsum(c(0,x))))
    text(x=sens$xx[lab.x],y=cs[cbind(1:length(lab.x),lab.x)] + (cs.diff/2)[cbind(1:length(lab.x),lab.x)],sens$names.ind,...)
    sens.mean.var<-apply(sens$S.var,2:3,mean)
    matplot(sens$xx[ord],t(apply(sens.mean.var,2,cumsum))[ord,],type='l',xlab='x',ylab='variance',main='Varaince Decomposition',...)
  } else{
    boxplot(sens$S,las=2,ylab='proportion varaince',main='Sensitivity',range=0,...)
    boxplot(sens$T,main='Total Sensitivity',range=0,...)
  }
  par(op)
}

