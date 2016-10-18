#' @title BASS Plot Diagnostics
#'
#' @description Generate diagnostic plots for BASS model fit.
#' @param mod a \code{bass} object.
#' @param quants quantiles for intervals, if desired.  NULL if not desired.
#' @param ... graphical parameters.
#' @details The first plot is posterior predicted vs observed, with intervals for predictions.  The second plot is a histogram of the residuals (of the posterior mean model), with a red curve showing the assumed Normal density (with posterior mean variance).
#' @keywords BMARS
#' @export
#' @seealso \link{bass}, \link{predict.bass}, \link{sobol}
#' @examples
#' # use example path/relative/to/packge/root instead
#'
plot.bass<-function(mod,quants=c(.025,.975),...){
  op<-par(no.readonly=T)
  par(mfrow=c(1,2))
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

# print.BMARS
# summary.BMARS

# plot.sobolBMARS, plot.predictBMARS
