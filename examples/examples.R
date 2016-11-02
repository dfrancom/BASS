### univariate example

## simulate data (Friedman function)
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma<-1 # noise sd
n<-500 # number of observations
x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y<-rnorm(n,f(x),sigma)

## fit BMARS
mod<-bass(x,y,maxInt=3,nmcmc=30000,nburn=20000,thin=10,npart=20)#,h2=1000,temp.ladder = c(1,.8,.65,.5,.4,.3,.23,.15,.1,.05),start.temper=1000)
plot(mod) # plot fit

## prediction
npred<-1000
xpred<-matrix(runif(npred*10),npred,10)
pred<-predict(mod,xpred,verbose = T) # posterior predictive samples
true.y<-f(xpred)
plot(true.y,colMeans(pred),xlab='true values',ylab='posterior predictive means'); abline(a=0,b=1,col=2)
quants.pred<-apply(pred,2,quantile,probs=c(.025,.975))
segments(true.y,quants.pred[1,],true.y,quants.pred[2,],col='lightgrey')
points(true.y,colMeans(pred)); abline(a=0,b=1,col=2)
mean(true.y>quants.pred[1,] & true.y<quants.pred[2,]) # empirical coverage of 95% interval - undercoverage indicates tempering could help


## sobol
sens<-sobol(mod)
plot(sens,cex.axis=.5)


### functional example

## simulate data (Friedman function with first variable as functional)
sigma<-1 # noise sd
n<-500 # number of observations
nfunc<-50
xfunc<-seq(0,1,length.out=nfunc)
x<-matrix(runif(n*9),n,9) #10 variables (1 functional), only first 5 matter
y<-matrix(f(cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x))),nrow=nfunc,byrow=T)+rnorm(n*nfunc,0,sigma)
mod<-bass(x,y,xx.func=xfunc,maxInt=3,nmcmc=20000,nburn=19000)
plot(mod)

## prediction
npred<-100
xpred<-matrix(runif(npred*9),npred,9)
ypred<-matrix(f(cbind(rep(xfunc,each=npred),kronecker(rep(1,nfunc),xpred))),nrow=nfunc,byrow=T)
pred<-predict(mod,xpred) # posterior predictive samples (each is a curve)
matplot(ypred,t(apply(pred,2:3,mean)),type='l',xlab='observed',ylab='mean prediction'); abline(a=0,b=1,col=2) # true values against posterior predictive means
matplot(ypred,type='l') # actual
matplot(t(apply(pred,2:3,mean)),type='l') # mean prediction


## sobol
sens<-sobol(mod,mcmc.use = 1:10) # speed this up by specifing a few samples (mcmc.use)
plot(sens) # note that the functional variable(s) are appended to the end of the list of variables (labeled 10 here)

sens.func<-sobol(mod,mcmc.use=1:10,func.var=1) # speed this up by specifing a few samples (mcmc.use)
plot(sens.func)
