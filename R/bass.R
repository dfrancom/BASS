
########################################################################
## make basis functions
########################################################################

pos<-function(vec){ # makes negative values 0
  replace(vec,vec<0,0)
}
const<-function(signs,knots,degree){ # largest value of basis function, assuming x's in [0,1], used for scaling
  cc<-prod(((signs+1)/2 - signs*knots))^degree
  if(cc==0)
    return(1)
  return(cc)
} # since a product, can find for functional & categorical pieces separately, take product
makeBasis<-function(signs,vars,knots,datat,degree){ #faster than apply
  cc<-const(signs,knots,degree)
  temp1<-pos(signs*(datat[vars,]-knots))^degree
  if(length(vars)==1){
    return(temp1/cc)
  } else{
    temp2<-1
    for(pp in 1:length(vars)){
      temp2<-temp2*temp1[pp,]
    }
    return(temp2/cc)
  }
}

makeBasisCat<-function(vars,sub,data){
  temp<-1
  for(ii in 1:length(vars)){
    temp<-temp*as.numeric(data[,vars[ii]] %in% sub[[ii]])
  }
  #browser()
  return(temp)
}

########################################################################
## functions used in MCMC
########################################################################

# CHANGE FOR FUNC
logProbChangeMod<-function(n.int,vars,I.vec,z.vec,p,vars.len,maxInt,miC){ # reversibility term
  if(n.int==1){ #for acceptance ratio
    out<-log(I.vec[n.int+miC])-log(2*p*vars.len[vars]) + #proposal
      log(2*p*vars.len[vars])+log(maxInt) # prior
  } else{
    # perms<-permutations(n.int,n.int,vars)
    # sum.perm<-sum(apply(perms,1,function(row){1/prod(1-cumsum(z.vec[row][-n.int]))}))
    # lprob.vars.noReplace<-sum(log(z.vec[vars]))+log(sum.perm)
    #require(BiasedUrn) # this is much faster than above (esp for large maxInt)
    x<-rep(0,p)
    x[vars]<-1
    #lprob.vars.noReplace<-log(BiasedUrn::dMWNCHypergeo(x,rep(1,p),n.int,z.vec)) - do this in combination with imports: BiasedUrn in DESCRIPTION file, but that has a limit to MAXCOLORS
    lprob.vars.noReplace<-log(dMWNCHypergeo(x,rep(1,p),n.int,z.vec))
    out<-log(I.vec[n.int+miC])+lprob.vars.noReplace-n.int*log(2)-sum(log(vars.len[vars])) + # proposal
      +n.int*log(2)+sum(log(vars.len[vars]))+lchoose(p,n.int)+log(maxInt) # prior
  }
  return(out)
}


logProbChangeModCat<-function(n.int,vars,I.vec,z.vec,p,nlevels,sub.size,maxInt,miC){
  if(n.int==1){ #for acceptance ratio
    out<-log(I.vec[n.int+miC])-log(p*(nlevels[vars]-1))-lchoose(nlevels[vars],sub.size[1:n.int]) + # proposal
      log(p*(nlevels[vars]-1))+lchoose(nlevels[vars],sub.size[1:n.int])+log(maxInt) # prior
  } else{
    x<-rep(0,p)
    x[vars]<-1
    lprob.vars.noReplace<-log(dMWNCHypergeo(x,rep(1,p),n.int,z.vec))
    out<-log(I.vec[n.int+miC])+lprob.vars.noReplace-n.int*sum(log(nlevels[vars]-1))-sum(lchoose(nlevels[vars],sub.size[1:n.int])) + # proposal
      n.int*sum(log(nlevels[vars]-1))+sum(lchoose(nlevels[vars],sub.size[1:n.int]))+lchoose(p,n.int)+log(maxInt) # prior
  }
  if(length(out)>1)
    browser()
  if(is.na(out))
    browser()
  return(out)
}

# CHANGE FOR FUNC
lp<-function(curr,prior,data){ # log posterior
  if(curr$nbasis==0)
    return(NA)
  tt<-(
    - (curr$s2.rate+prior$g2)/curr$s2
    -(data$n/2+1+(curr$nbasis+1)/2 -prior$g1)*log(curr$s2)
    + sum(log(abs(diag(curr$R))))
    + (prior$a.beta.prec+(curr$nbasis+1)/2-1)*log(curr$beta.prec) - prior$b.beta.prec*curr$beta.prec
    - (curr$nbasis+1)/2*log(2*pi)
    + (prior$h1+curr$nbasis-1)*log(curr$lam) - curr$lam*(prior$h2+1)
    )
    #priors for basis parameters
  if(data$des){
    tt<-tt+(
      - sum(curr$n.int.des)*log(2)
      - sum(lchoose(data$pdes,curr$n.int.des))
      - sum(log(data$vars.len.des[na.omit(curr$vars.des)]))
      - curr$nbasis*log(prior$maxInt.des)
    )
  }
  if(data$cat){
    tt<-tt+(
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CHECK THIS
      - sum(sapply(1:curr$nbasis,function(i) curr$n.int.cat[i]*sum(log(data$nlevels[na.omit(curr$vars.cat[i,])]-1))))
      - sum(sapply(1:curr$nbasis,function(i) sum(lchoose(data$nlevels[na.omit(curr$vars.cat[i,])],curr$sub.size[i,1:curr$n.int.cat[i]]))))
      - sum(lchoose(data$pcat,curr$n.int.cat))
      - curr$nbasis*log(prior$maxInt.cat)
    )
  }
  if(data$func){
    tt<-tt+(
      - sum(curr$n.int.func)*log(2)
      - sum(lchoose(data$pfunc,curr$n.int.func))
      - sum(log(data$vars.len.func[na.omit(curr$vars.func)]))
      - curr$nbasis*log(prior$maxInt.func)
    )
  }

  return(tt)
}

getQf<-function(XtX,Xty){ # get quadratic from that shows up in acceptance probability
  R<-tryCatch(chol(XtX), error=function(e) matrix(F))
  if(R[1,1]){
    dr<-diag(R)
    #if(min(dr)<1e-8) ## hack, otherwise can add the exact same basis sometimes (not sure why)
    if(length(dr)>1){
      if(max(dr[-1])/min(dr)>1e3) # TODO: this is a hack, otherwise we get huge variance inflation in beta
        return(NULL)
    }
    bhat<-backsolve(R,forwardsolve(R,Xty,transpose=T,upper.tri=T))
    qf<-crossprod(bhat,Xty)# same as sum((R%*%bhat)^2)
    return(list(R=R,bhat=bhat,qf=qf))
  } else{
    return(NULL)
  }
}

rgammaTemper<-function(n,shape,rate,temper){ # sample a tempered gamma
  rgamma(n,temper*(shape-1)+1,temper*rate)
}
rigammaTemper<-function(n,shape,scale,temper){ # sample a tempered IG
  1/rgamma(n,temper*(shape+1)-1,rate=temper*scale)
}

########################################################################
## MCMC update
########################################################################


updateMCMC<-function(curr,prior,data,funcs=funcs){

  ## RJMCMC update

  u<-sample(1:3,size=1)
  if(curr$nbasis==0){
    u<-1 # birth for sure
  }
  if(curr$nbasis==prior$maxBasis){
    u<-sample(2:3,size=1) # no birth
  }
  if(u==1){ # birth
    curr<-funcs$birth(curr,prior,data)
  } else if(u==2){ # death
    curr<-funcs$death(curr,prior,data)
  } else{ # change
    curr<-funcs$change(curr,prior,data)
  }

  #print(c(curr$qf))
  ## Gibbs updates

  # beta
  curr$beta<-curr$bhat/(1+curr$beta.prec)+curr$R.inv.t%*%rnorm(curr$nc)*sqrt(curr$s2/(1+curr$beta.prec)/data$temp.ladder[curr$temp.ind])

  # lambda
  lam.a<-prior$h1+curr$nbasis
  lam.b<-prior$h2+1
  curr$lam<-rgammaTemper(1,lam.a,lam.b,data$temp.ladder[curr$temp.ind])

  # s2
  qf2<-crossprod(curr$R%*%curr$beta)
  curr$s2.rate<-(data$ssy + (1+curr$beta.prec)*qf2 - 2*crossprod(curr$beta,curr$Xty[1:curr$nc]))/2
  s2.a<-prior$g1+(data$n+curr$nbasis+1)/2
  s2.b<-prior$g2+curr$s2.rate
  curr$s2<-rigammaTemper(1,s2.a,s2.b,data$temp.ladder[curr$temp.ind])
  if(is.nan(curr$s2) | is.na(curr$s2)) # major variance inflation, get huge betas from curr$R.inv.t, everything becomes unstable
    browser()
  if(curr$s2==0){ # tempering instability, this temperature too small
    curr$s2<-runif(1,0,1e6)
    #browser()
  }


  # beta.prec
  beta.prec.a<-prior$a.beta.prec+(curr$nbasis+1)/2
  beta.prec.b<-prior$b.beta.prec+1/(2*curr$s2)*qf2
  curr$beta.prec<-rgammaTemper(1,beta.prec.a,beta.prec.b,data$temp.ladder[curr$temp.ind])

  ## save log posterior
  curr$lpost<-lp(curr,prior,data) # doesn't include cat yet

  return(curr)
}
scale.range<-function(x,r=NULL){ # x is a vector
  if(is.null(r))
    r<-range(x)
  (x-r[1])/(r[2]-r[1])
}
unscale.range<-function(x,r){
  x*(r[2]-r[1])+r[1]
}

getYhat_des<-function(curr,nb){
  curr$des.basis%*%curr$beta
}
getYhat_cat<-function(curr,nb){
  curr$cat.basis%*%curr$beta
}
getYhat_des_cat<-function(curr,nb){
  curr$dc.basis%*%curr$beta
}
getYhat_des_func<-function(curr,nb){
  tcrossprod(curr$des.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_cat_func<-function(curr,nb){
  tcrossprod(curr$cat.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}
getYhat_des_cat_func<-function(curr,nb){
  tcrossprod(curr$dc.basis%*%diag(c(curr$beta),nb+1),curr$func.basis)
}

########################################################################
## bass function
########################################################################
#' @title Bayesian Adaptive Spline Surfaces (BASS)
#'
#' @description Fits a BASS model using RJMCMC.  Optionally uses parallel tempering to improve mixing.  Can be used with scalar or functional response.  Also can use categorical inputs.
#' @param xx a data frame or matrix of predictors.  Categorical predictors should be included as factors.
#' @param y  a response vector (scalar response) or matrix (functional response).
#' @param maxInt integer for maximum degree of interaction in spline basis functions.  Defaults to the number of predictors, which could result in overfitting.
#' @param maxInt.func (functional response only) integer for maximum degree of interaction in spline basis functions describing the functional response.
#' @param maxInt.func (categorical input only) integer for maximum degree of interaction of categorical inputs.
#' @param xx.func a vector, matrix or data frame of functional variables.
#' @param degree degree of splines.  Stability should be examined for anything other than 1.
#' @param maxBasis maximum number of basis functions.
#' @param npart minimum number of non-zero points in a basis function.  If the response is functional, this refers only to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the number of observations, whichever is smaller.
#' @param npart.func same as npart, but for functional portion of basis function.
#' @param nmcmc number of RJMCMC iterations.
#' @param nburn number of the \code{nmcmc} iterations to disregard.
#' @param thin keep every \code{thin} samples
#' @param g1 shape for IG prior on \eqn{\sigma^2}.
#' @param g2 scale for IG prior on \eqn{\sigma^2}.
#' @param h1 shape for gamma prior on \eqn{\lambda}.
#' @param h2 rate for gamma prior on \eqn{\lambda}.  This is the primary way to control overfitting.  A large value of \code{h2} favors fewer basis functions.
#' @param a.beta.prec shape for gamma prior on \eqn{\tau}.
#' @param b.beta.prec rate for gamma prior on \eqn{\tau}. Defaults to one over the number of observations, which is the unit information prior.
#' @param w1 nominal weight for degree of interaction, used in generating candidate basis functions.
#' @param w2 nominal weight for variables, used in generating candidate basis functions.
#' @param temp.ladder temperature ladder used for parallel tempering.  The first value should be 1 and the values should decrease.
#' @param start.temper when to start tempering (after how many MCMC iterations).
#' @param ncores currently disregarded.
#' @param curr.list list of starting models (one element for each temperature), could be output from a previous run under the same model setup.
#' @param save.yhat logical; should predictions of training data be saved?
#' @param verbose logical; should progress be displayed?
#' @details Explores BASS model space by RJMCMC.  The BASS model has \deqn{y = f(x) + \epsilon,  \epsilon ~ N(0,\sigma^2)} \deqn{f(x) = a_0 + \sum_{m=1}^M a_m B_m(x)} and \eqn{B_m(x)} is a BASS basis function (tensor product of spline basis functions). We use priors \deqn{a ~ N(0,\sigma^2/\tau (B'B)^{-1})} \deqn{M ~ Poisson(\lambda)} as well as the priors mentioned in the arguments above.
#' @return An object of class 'bass'.  The other output will only be useful to the advanced user.  Rather, users may be interested in prediction and sensitivity analysis, which are obtained by passing the entire object to the predict.bass or sobol functions.
#' @keywords BMARS
#' @seealso \link{predict.bass} for prediction and \link{sobol} for sensitivity analysis.
#' @export
#' @import stats
#' @import utils
#' @example examples/examples.R
#'
bass<-function(xx,y,maxInt=3,maxInt.func=3,maxInt.cat=3,xx.func=NULL,degree=1,maxBasis=1000,npart=NULL,npart.func=NULL,nmcmc=10000,nburn=9000,thin=1,g1=0,g2=0,h1=10,h2=10,a.beta.prec=1,b.beta.prec=NULL,w1=5,w2=5,temp.ladder=NULL,start.temper=NULL,ncores=1,curr.list=NULL,save.yhat=TRUE,verbose=TRUE){

  ########################################################################
  ## setup



  # there always has to be an xx in this version.  Version that allows for multiple functional variables could have no xx if xx.func (list of 2) included.


  xx<-as.data.frame(xx)
  dx<-dim(xx)
  dxf<-dim(xx.func)
  dy<-dim(y) # need to change this for array input
  if(any(dy==1))
    y<-c(y)
  dy<-dim(y)

  if(is.null(dy)){
    func<-F
    pfunc<-0
    if(!is.null(xx.func))
      warning('xx.func ignored because there is no functional variable')
    if(length(y)!=dx[1])
      stop('dimension mismatch between xx and y')
  } else {
    func<-T
    if(is.null(xx.func))
      stop('missing xx.func')
    xx.func<-as.matrix(xx.func)
    dxf<-dim(xx.func)
    if(dy[1]!=dx[1]){
      y<-t(y)
      dy<-dim(y)
    }
    if(dy[1]!=dx[1])
      stop('dimension mismatch between xx and y')
    if(dy[2]!=dxf[1])
      xx.func<-t(xx.func)
    dxf<-dim(xx.func)
    if(dy[2]!=dxf[1])
      stop('dimension mismatch between xx.func and y')
    pfunc<-dxf[2]
    range.func<-apply(xx.func,2,range)
    xx.func<-apply(xx.func,2,scale.range)
  }



  des<-T
  cx<-sapply(xx,class)
  cx.factor<- cx == 'factor'
  if(any(cx.factor)){
    cat<-T
    if(all(cx.factor))
      des<-F
    xx.des<-as.matrix(xx[,!cx.factor,drop=F])
    xx.cat<-xx[,cx.factor,drop=F]
  } else{
    cat<-F
    xx.des<-as.matrix(xx)
    xx.cat<-NULL
  }
  if(des){
    range.des<-apply(xx.des,2,range)
    xx.des<-apply(xx.des,2,scale.range)
  }
  des.vars<-which(!cx.factor)
  cat.vars<-which(cx.factor)
  pdes<-length(des.vars)
  pcat<-length(cat.vars)
  # do we even need xx.des and xx.cat if we have these indexes?

  type<-''
  if(des)
    type<-paste(type,'des',sep='_')
  if(cat)
    type<-paste(type,'cat',sep='_')
  if(func)
    type<-paste(type,'func',sep='_')

  # so cases are des,des_cat,des_cat_func,cat,cat_func

  if(is.null(temp.ladder)){
    temp.ladder<-1
  }
  if(min(temp.ladder)<(2/dx[1])){
    temp.ladder<-temp.ladder[temp.ladder>(2/dx[1])]
    if(length(temp.ladder)==0)
      stop('invalid temp.ladder (temperatures too small)')
  }
  ntemps<-length(temp.ladder)
  if(ntemps==1){
    start.temper<-nmcmc
  }
  temp.val<-matrix(nrow=nmcmc,ncol=ntemps)







  # data object
  data<-list()
  data$y<-y
  if(des){
    data$xxt.des<-t(xx.des)
    data$vars.len.des<-NA
    data$xxt.des.unique<-list()
    data$unique.ind.des<-list()
    for(i in 1:pdes){
      data$xxt.des.unique[[i]]<-unique(data$xxt.des[i,])
      data$unique.ind.des[[i]]<-which(!duplicated(data$xxt.des[i,])) # gets the first instance of each unique value
      data$vars.len.des[i]<-length(data$xxt.des.unique[[i]])
    }
  }
  if(func){
    data$xxt.func<-t(xx.func)
    data$vars.len.func<-NA
    data$xxt.func.unique<-list()
    data$unique.ind.func<-list()
    for(i in 1:pfunc){
      data$xxt.func.unique[[i]]<-unique(data$xxt.func[i,])
      data$unique.ind.func[[i]]<-which(!duplicated(data$xxt.func[i,])) # gets the first instance of each unique value
      data$vars.len.func[i]<-length(data$xxt.func.unique[[i]])
    }
  }
  if(cat){
    data$levels<-lapply(xx.cat,levels)
    data$nlevels<-sapply(data$levels,length)
    data$xx.cat<-xx.cat
  }
  data$pdes<-pdes
  data$pfunc<-pfunc
  data$pcat<-pcat
  data$p<-dx[2]
  data$ndes<-dx[1]
  data$nfunc<-dxf[1]
  data$des<-des
  data$func<-func
  data$cat<-cat

  data$n<-prod(data$ndes,data$nfunc)
  data$ssy<-sum(data$y^2)
  data$death.prob.next<-1/3
  data$birth.prob<-1/3
  data$birth.prob.last<-1/3
  data$death.prob<-1/3
  data$temp.ladder<-temp.ladder

  npart.des<-npart
  if(is.null(npart.des)){
    npart.des<-min(20,.1*data$n)
  }
  if(is.null(npart.func) & func){
    npart.func<-min(20,.1*data$nfunc)
  }

    

  maxBasis<-min(maxBasis,data$n) # can't have more basis functions than data points
  maxInt.des<-min(maxInt,pdes) # can't have more interactions than variables #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ may want to change this...
  maxInt.cat<-min(maxInt.cat,pcat)
  maxInt.func<-min(maxInt.func,pfunc)

  # prior object
  prior<-list()
  prior$maxInt.des<-maxInt.des
  prior$maxInt.cat<-maxInt.cat
  prior$maxInt.func<-maxInt.func
  prior$q<-degree
  prior$npart.des<-npart.des
  prior$npart.func<-npart.func
  prior$h1<-h1
  prior$h2<-h2
  prior$g1<-g1
  prior$g2<-g2
  prior$a.beta.prec<-a.beta.prec
  if(is.null(b.beta.prec)){
    prior$b.beta.prec<-1/data$n
  } else{
  prior$b.beta.prec<-b.beta.prec
  }
  prior$maxBasis<-maxBasis
  prior$minInt<-0
  if(des+cat+func==1) # if there is only one part, can't have minInt of 0
    prior$minInt<-1
  prior$miC<-abs(prior$minInt-1)

  # current mcmc state (can use from previous)
  if(is.null(curr.list)){
    curr.list<-list()
    for(i in 1:ntemps){
      curr.list[[i]]<-list()

      if(des){
        curr.list[[i]]$I.star.des<-rep(w1,prior$maxInt.des+prior$miC)
        curr.list[[i]]$I.vec.des<-curr.list[[i]]$I.star.des/sum(curr.list[[i]]$I.star.des)
        curr.list[[i]]$z.star.des<-rep(w2,data$pdes)
        curr.list[[i]]$z.vec.des<-curr.list[[i]]$z.star.des/sum(curr.list[[i]]$z.star.des)
        curr.list[[i]]$des.basis<-matrix(rep(1,data$ndes))
      }
      if(cat){
        curr.list[[i]]$I.star.cat<-rep(w1,prior$maxInt.cat+prior$miC)
        curr.list[[i]]$I.vec.cat<-curr.list[[i]]$I.star.cat/sum(curr.list[[i]]$I.star.cat)
        curr.list[[i]]$z.star.cat<-rep(w2,data$pcat)
        curr.list[[i]]$z.vec.cat<-curr.list[[i]]$z.star.cat/sum(curr.list[[i]]$z.star.cat)
        curr.list[[i]]$cat.basis<-matrix(rep(1,data$ndes))
      }
      if(func){
        curr.list[[i]]$I.star.func<-rep(w1,prior$maxInt.func+prior$miC)
        curr.list[[i]]$I.vec.func<-curr.list[[i]]$I.star.func/sum(curr.list[[i]]$I.star.func)
        curr.list[[i]]$z.star.func<-rep(w2,data$pfunc)
        curr.list[[i]]$z.vec.func<-curr.list[[i]]$z.star.func/sum(curr.list[[i]]$z.star.func)
        curr.list[[i]]$func.basis<-matrix(rep(1,data$nfunc))
      }

      if(des & cat)
        curr.list[[i]]$dc.basis<-curr.list[[i]]$des.basis*curr.list[[i]]$cat.basis

      curr.list[[i]]$s2<-1
      curr.list[[i]]$lam<-1
      curr.list[[i]]$beta.prec<-1
      curr.list[[i]]$nbasis<-0
      curr.list[[i]]$nc<-1

      curr.list[[i]]$knots.des<-matrix(numeric(0),ncol=maxInt.des)
      curr.list[[i]]$knotInd.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$signs.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$vars.des<-matrix(integer(0),ncol=maxInt.des)
      curr.list[[i]]$n.int.des<-NA

      curr.list[[i]]$sub.list<-list()
      curr.list[[i]]$sub.size<-matrix(integer(0),ncol=maxInt.cat)
      curr.list[[i]]$vars.cat<-matrix(integer(0),ncol=maxInt.cat)
      curr.list[[i]]$n.int.cat<-NA

      curr.list[[i]]$knots.func<-matrix(numeric(0),ncol=maxInt.func)
      curr.list[[i]]$knotInd.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$signs.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$vars.func<-matrix(integer(0),ncol=maxInt.func)
      curr.list[[i]]$n.int.func<-NA

      curr.list[[i]]$Xty<-rep(NA,maxBasis+2)
      curr.list[[i]]$Xty[1]<-sum(data$y)
      curr.list[[i]]$XtX<-matrix(NA,nrow=maxBasis+2,ncol=maxBasis+2)
      curr.list[[i]]$XtX[1,1]<-data$n
      curr.list[[i]]$R<-chol(curr.list[[i]]$XtX[1,1])
      curr.list[[i]]$R.inv.t<-t(solve(curr.list[[i]]$R))
      curr.list[[i]]$bhat<-mean(data$y)
      curr.list[[i]]$qf<-crossprod(curr.list[[i]]$R%*%curr.list[[i]]$bhat)
      curr.list[[i]]$count<-rep(0,3)
      curr.list[[i]]$cmod<-F
      curr.list[[i]]$step<-NA
      curr.list[[i]]$temp.ind<-i
      curr.list[[i]]$type<-type
    }
  }

  # define functions according to type.  Doing eval parse every time the functions were used is slow.
  funcs<-list()
  funcs$birth<-eval(parse(text=paste('birth',type,sep='')))
  funcs$death<-eval(parse(text=paste('death',type,sep='')))
  funcs$change<-eval(parse(text=paste('change',type,sep='')))
  funcs$getYhat<-eval(parse(text=paste('getYhat',type,sep='')))
  
  # CHANGE FOR FUNC
  # do something like maxInt.tot<-maxInt.des+maxInt.func+maxInt.cat
  # then have an index, int.des<-1:maxInt.des; int.func<-maxInt.des + 1:maxInt.func; int.cat<-maxInt.des+maxInt.func
  nmod.max<-(nmcmc-nburn)/thin # max number of models (models don't necessarily change every iteration)
  if(des){
    signs.des<-knotInd.des<-vars.des<-array(dim=c(nmod.max,maxBasis,maxInt.des)) # truncate when returning at end of function
    n.int.des<-matrix(nrow=nmod.max,ncol=maxBasis) # degree of interaction
  }
  if(cat){
    sub.list<-list() # this is big...
    sub.size<-vars.cat<-array(dim=c(nmod.max,maxBasis,maxInt.cat))
    n.int.cat<-matrix(nrow=nmod.max,ncol=maxBasis)
  }
  if(func){
    signs.func<-knotInd.func<-vars.func<-array(dim=c(nmod.max,maxBasis,maxInt.func))
  # arrays use less space, esp integer arrays
    n.int.func<-matrix(nrow=nmod.max,ncol=maxBasis)
  }

  beta<-matrix(nrow=nmod.max,ncol=maxBasis+1) # +1 for intercept, nmcmc-nburn instead of nmod because beta updates every iteration
  nbasis<-s2<-lam<-beta.prec<-NA
  cmod<-F # indicator for whether we have changed models since last storing
  model.lookup<-NA # lookup table between models and mcmc iterations

  log.post.cold<-rep(NA,nmcmc) # log posterior for cold chain (the one we care about)

  if(save.yhat){
    yhat.sum<-0#rep(0,data$n) # if we don't want to store all yhat draws, can still get running average
    yhat<-array(dim=c(nmod.max,data$ndes,data$nfunc))#matrix(nrow=nmod.max,ncol=data$n)
  }

  # temperature index
  cold.chain<-1 # to start, the cold chain is curr.list[[1]]
  temp.ind<-1:ntemps # we will change this vector as we swap temperatures
  count.swap<-count.swap.disp<-rep(0,ntemps-1) # number of swaps between each set of neighbors
  swap<-NA # to keep track of swaps
  #require(parallel) # for tempering


  ########################################################################
  ## MCMC

  if(verbose)
    cat('MCMC Start',timestamp(prefix='#--',suffix='--#',quiet=T),'nbasis:',curr.list[[cold.chain]]$nbasis,'\n')
  n.models<-keep.sample<-0 # indexes for storage
  for(i in 2:nmcmc){

    ## update model for each temperature
    #curr.list<-parLapply(cluster,curr.list,updateMCMC)
    ncores<-1
    curr.list<-parallel::mclapply(curr.list,updateMCMC,prior=prior,data=data,funcs=funcs,mc.preschedule=T,mc.cores=ncores)
    # TODO: DO SOMETHING LIKE THIS BUT KEEP EVERYTHING SEPARATE ON THE CLUSTER, all we need is lpost, cmod
    
    #if(i%%1000==0)
    #plot(c(data$xxt.des),data$y,col=as.numeric(unlist(data$xx.cat)))
    #points(c(data$xxt.des),curr.list[[1]]$dc.basis%*%curr.list[[1]]$beta,col=as.numeric(unlist(data$xx.cat)),cex=.5)
    #browser()
    

    ## parallel tempering swap
    if(i>start.temper){# & (i%%20==0)){ #only start after a certain point, and only try every 20
      # sample temp.ind.swap from 1:(ntemps-1), then swap with temp.ind.swap+1
      temp.ind.swap1<-sample(1:(ntemps-1),size=1) # corresponds to temperature temp.ladder[temp.ind.swap1]
      temp.ind.swap2<-temp.ind.swap1+1 # always use the neighboring chain on the right
      chain.ind1<-which(temp.ind==temp.ind.swap1) # which chain has temperature temp.ladder[temp.ind.swap1]
      chain.ind2<-which(temp.ind==temp.ind.swap2)
      alpha.swap<-(data$temp.ladder[temp.ind.swap1]-data$temp.ladder[temp.ind.swap2])*(curr.list[[chain.ind2]]$lpost-curr.list[[chain.ind1]]$lpost)
      if(log(runif(1)) < alpha.swap){
        # swap temperatures
        temp.ind[chain.ind1]<-temp.ind.swap2
        temp.ind[chain.ind2]<-temp.ind.swap1
        curr.list[[chain.ind1]]$temp.ind<-temp.ind.swap2
        curr.list[[chain.ind2]]$temp.ind<-temp.ind.swap1

        count.swap[temp.ind.swap1]<-count.swap[temp.ind.swap1]+1
        count.swap.disp[temp.ind.swap1]<-count.swap.disp[temp.ind.swap1]+1
        swap[i]<-temp.ind.swap1
        if(temp.ind.swap1==1){
          cmod<-T # we changed models
          cold.chain<-chain.ind2 #which(temp.ind==1)
        }
      }
    }

    log.post.cold[i]<-curr.list[[cold.chain]]$lpost
    temp.val[i,]<-temp.ind



    ## write current model if past burnin and model is unique
    # if(i==nburn & test.data){
    #   # make test basis functions
    # }
    if((i>nburn) & (((i-nburn)%%thin)==0)){
      # these things are updated every time
      keep.sample<-keep.sample+1 # indexes samples
      nb<-curr.list[[cold.chain]]$nbasis
      nbasis[keep.sample]<-nb
      beta[keep.sample,1:(nb+1)]<-curr.list[[cold.chain]]$beta
      s2[keep.sample]<-curr.list[[cold.chain]]$s2
      lam[keep.sample]<-curr.list[[cold.chain]]$lam
      beta.prec[keep.sample]<-curr.list[[cold.chain]]$beta.prec
      if(save.yhat){
        yhat.current<-funcs$getYhat(curr.list[[cold.chain]],nb)
        if(func){
          yhat[keep.sample,,]<-yhat.current
        } else{
          yhat[keep.sample,]<-yhat.current
        }
        yhat.sum<-yhat.sum+yhat.current

      }
      # save cold chain basis parms if they are different from previous (cmod=T)
      if(cmod || curr.list[[cold.chain]]$cmod){ # can I actually get curr.list[[cold.chain]]$cmod easily from the core it is on?
        n.models<-n.models+1 # indexes models
        if(nb>0){
          if(des){
            vars.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.des)
            signs.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$signs.des)
            knotInd.des[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$knotInd.des)
            n.int.des[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.des)
          }
          if(cat){
            vars.cat[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.cat)
            sub.size[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$sub.size)
            sub.list[[n.models]]<-curr.list[[cold.chain]]$sub.list
            n.int.cat[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.cat)
          }
          if(func){
            vars.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$vars.func)
            signs.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$signs.func)
            knotInd.func[n.models,1:nb,]<-as.integer(curr.list[[cold.chain]]$knotInd.func)
            n.int.func[n.models,1:nb]<-as.integer(curr.list[[cold.chain]]$n.int.func)
          }
        }
        cmod<-F # reset change model indicator after writing current model
        curr.list[[cold.chain]]$cmod<-F
      }
      model.lookup[keep.sample]<-n.models # update lookup table
    }

    if(verbose & i%%1000==0){
      pr<-c('MCMC iteration',i,timestamp(prefix='#--',suffix='--#',quiet=T),'nbasis:',curr.list[[cold.chain]]$nbasis)
      if(i>start.temper)
        pr<-c(pr,'tempering acc',count.swap.disp/(1000/(ntemps-1)))
      cat(pr,'\n')
      count.swap.disp<-rep(0,ntemps-1) # reset after displaying
    }

  }

  ########################################################################
  ## return

  out.yhat<-list()
  if(save.yhat){
    out.yhat<-list(yhat.mean=yhat.sum/nmod.max,yhat=yhat,y=y)
  }

  out<-list(
       beta=beta,
       s2=s2,
       lam=lam,
       nbasis=nbasis,
       degree=degree,
       nmcmc=nmcmc,
       nburn=nburn,
       thin=thin,
       p=data$p,
       beta.prec=beta.prec,
       step=step,
       y=y,
       log.post.cold=log.post.cold,
       curr.list=curr.list, # for restarting
       swap=swap,
       count.swap=count.swap,
       temp.val=temp.val,
       n.models=n.models,
       model.lookup=model.lookup,
       des=des,func=func,cat=cat,type=type,cx=cx
  )

  mb<-max(nbasis)
  
  out.des<-list()
  if(des){
    out.des<-list(
      knotInd.des=knotInd.des[1:n.models,1:mb,,drop=F],
      signs.des=signs.des[1:n.models,1:mb,,drop=F],
      vars.des=vars.des[1:n.models,1:mb,,drop=F],
      n.int.des=n.int.des[1:n.models,1:mb,drop=F],
      maxInt.des=maxInt.des,
      des.basis=curr.list[[cold.chain]]$des.basis,
      pdes=pdes,
      xx.des=xx.des,range.des=range.des,
      unique.ind.des=data$unique.ind.des
    )
  }

  out.cat<-list()
  if(cat){
    out.cat<-list(
      vars.cat=vars.cat[1:n.models,1:mb,,drop=F],
      sub.size=sub.size[1:n.models,1:mb,,drop=F],
      sub.list=sub.list,
      n.int.cat=n.int.cat[1:n.models,1:mb,drop=F],
      maxInt.cat=maxInt.cat,
      cat.basis=curr.list[[cold.chain]]$cat.basis,
      pcat=pcat,
      xx.cat=xx.cat,
      nlevels=data$nlevels
    )
  }

  out.func<-list()
  if(func){
    out.func<-list(
      knotInd.func=knotInd.func[1:n.models,1:mb,,drop=F],
      signs.func=signs.func[1:n.models,1:mb,,drop=F],
      vars.func=vars.func[1:n.models,1:mb,,drop=F],
      n.int.func=n.int.func[1:n.models,1:mb,drop=F],
      maxInt.func=maxInt.func,
      func.basis=curr.list[[cold.chain]]$func.basis,
      pfunc=pfunc,
      xx.func=xx.func,range.func=range.func,
      unique.ind.func=data$unique.ind.func
    )
  }

  #stopCluster(cluster)
  ret<-c(out.yhat,out,out.des,out.cat,out.func)
  class(ret)<-'bass'
  return(ret)
}