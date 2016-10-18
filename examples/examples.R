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
bm<-BMARS(x,y,maxInt=3,nmcmc=20000,nburn=19000)
plot(bm) # plot fit

## prediction
npred<-100
xpred<-matrix(runif(npred*10),npred,10)
pred<-predictBMARS(bm,xpred) # posterior predictive samples
plot(f(xpred),colMeans(pred)); abline(a=0,b=1,col=2) # true values against posterior predictive means

## sobol
ss<-sobolBMARS(bm)
boxplot(ss$S)
boxplot(ss$T)


### functional example

## simulate data (Friedman function with first variable as functional)
sigma<-1 # noise sd
n<-500 # number of observations
nfunc<-50
xfunc<-seq(0,1,length.out=nfunc)
x<-matrix(runif(n*9),n,9) #10 variables (1 functional), only first 5 matter
y<-matrix(f(cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x))),nrow=nfunc,byrow=T)+rnorm(n*nfunc,0,sigma)
bm<-BMARS(x,y,xx.func=xfunc,maxInt=3,nmcmc=20000,nburn=19000)
plot(bm)

## prediction
npred<-100
xpred<-matrix(runif(npred*9),npred,9)
ypred<-matrix(f(cbind(rep(xfunc,each=npred),kronecker(rep(1,nfunc),xpred))),nrow=nfunc,byrow=T)
pred<-predictBMARS(bm,xpred) # posterior predictive samples (each is a curve)
matplot(ypred,t(apply(pred,2:3,mean)),type='l',xlab='observed',ylab='mean prediction'); abline(a=0,b=1,col=2) # true values against posterior predictive means
matplot(ypred,type='l') # actual
matplot(t(apply(pred,2:3,mean)),type='l') # mean prediction


## sobol
ss<-sobolBMARS(bm,mcmc.use = 1:10) # speed this up by specifing a few samples (mcmc.use)
boxplot(ss$S)
boxplot(ss$T) # note that the functional variable(s) are appended to the end of the list of variables (labeled 10 here)

ss<-sobolBMARS(bm,mcmc.use=1:10,func.var=1) # speed this up by specifing a few samples (mcmc.use)
dim(ss$S)
matplot(t(apply(ss$S[1,,],2,cumsum)),type='l') # functional sensitivity indices for 1st posterior draw
ss.mean<-apply(ss$S,2:3,mean) # posterior mean functional sensitivity indices
matplot(t(apply(ss.mean,2,cumsum)),type='l')

ss.mean.scale<-apply(ss$S.var,2:3,mean)
matplot(t(apply(ss.mean.scale,2,cumsum)),type='l') # functional partitioning of variance (posterior mean)
