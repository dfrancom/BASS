########################################################################
## functions for prediction
########################################################################

## make basis functions for model i
makeBasisMatrix<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train){
  n<-ncol(xxt)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(nbasis>0){
    for(m in 1:nbasis){
      if(n.int[i,m]==0){
        tbasis.mat[m+1,]<-1 # could do this at beginning
      } else{
        use<-1:n.int[i,m]
        knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])] # get knots from knots.ind
        tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],vars[i,m,use],knots,xxt,q)
      }
    }
  }
  return(tbasis.mat)
} # I think prediction on a test set should be kept separate from the BMARS function, simplifies tempering (only would want to predict on the cold chain).

makeBasisMatrixCat<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,Xt){
  return(1)
}

#' @title BMARS Prediction
#'
#' @description Predict function for BMARS.  Outputs the posterior predictive samples based on the specified MCMC iterations.
#' @param bmars a fitted model, output from the \code{BMARS} function.
#' @param xx.pred a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{BMARS} function.
#' @param xx.pred.func a matrix of new values of the functional variable.  If none, the same values will be used as in the training data.
#' @param mcmc.use a vector indexing which MCMC iterations to use for prediction.
#' @param verbose logical; should progress be displayed?
#' @details Efficiently predicts when two MCMC iterations have the same basis functions (but different weights).
#' @return If model output is a scalar, this returns a matrix with the same number of rows as \code{xx.pred} and columns corresponding to the the MCMC iterations \code{mcmc.use}.  These are samples from the posterior predictive distribution.  If model output is functional, this returns an array with first dimension corresponding to MCMC iteration, second dimension corresponding to the rows of \code{xx.pred}, and third dimension corresponding to the rows of \code{xx.pred.func}.
#' @keywords BMARS
#' @seealso \link{BMARS} for model fitting and \link{sobolBMARS} for sensitivity analysis.
#' @export
#' @examples
#' # See examples in BMARS documentation.
#'
predictBMARS<-function(bmars,xx.pred,xx.pred.func=NULL,mcmc.use=NULL,verbose=FALSE){
  if(is.null(mcmc.use)){ # if null, use all
    mcmc.use<-1:((bmars$nmcmc-bmars$nburn)/bmars$thin)
  }
  if(bmars$func){
    if(is.null(xx.pred.func))
      xx.pred.func<-bmars$xx.func
    else{
      dxf<-dim(xx.pred.func)
      if(is.null(dxf))
        xx.pred.func<-matrix(xx.pred.func)
      for(i in 1:ncol(xx.pred.func)){
        xx.pred.func[,i]<-scale.range(xx.pred.func[,i],bmars$range.func[,i])
      }
    }
  } else{
    xx.pred.func<-t(1) # placeholder
  }

  xx.pred<-as.data.frame(xx.pred)

  for(i in 1:ncol(xx.pred)){
   xx.pred[,i]<-scale.range(xx.pred[,i],bmars$range.des[,i])
  }
  out<-array(dim=c(length(mcmc.use),nrow(xx.pred),nrow(xx.pred.func)))
  k<-0
  models<-bmars$model.lookup[mcmc.use]
  if(verbose)
    cat('Predict Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')
  mod.ind<-0
  for(j in unique(models)){ # loop though models, could be parallel?
    mod.ind<-mod.ind+1
    mcmc.use.j<-mcmc.use[models==j]
    ind<-k+(1:length(mcmc.use.j)) # index for storage
    k<-k+length(ind) # used for start of index
    out[ind,,]<-eval(parse(text=paste('mult',bmars$type,sep='')))(model=j,mcmc.use.mod=mcmc.use.j,bmars=bmars,xx.pred=xx.pred,xx.pred.func=xx.pred.func)
    if(verbose & mod.ind%%100==0)
      cat('Predict',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.ind,'\n')
  }
  return(drop(out))
}

mult_des<-function(model,mcmc.use.mod,bmars,xx.pred,xx.pred.func){
  M<-bmars$nbasis[mcmc.use.mod[1]]
  tmat<-makeBasisMatrix(model,M,bmars$vars,bmars$signs,bmars$knotInd,bmars$degree,t(xx.pred),bmars$n.int,bmars$xx.des) ## need to get rid of these transposes
  out<-bmars$beta[mcmc.use.mod,1:(M+1),drop=F]%*%tmat
  return(out)
}

mult_des_func<-function(model,mcmc.use.mod,bmars,xx.pred,xx.pred.func){
  M<-bmars$nbasis[mcmc.use.mod[1]]
  tmat.des<-makeBasisMatrix(model,M,bmars$vars.des,bmars$signs.des,bmars$knotInd.des,bmars$degree,t(xx.pred),bmars$n.int.des,bmars$xx.des)
  tmat.func<-makeBasisMatrix(model,M,bmars$vars.func,bmars$signs.func,bmars$knotInd.func,bmars$degree,t(xx.pred.func),bmars$n.int.func,bmars$xx.func)
  out<-array(dim=c(length(mcmc.use.mod),nrow(xx.pred),nrow(xx.pred.func)))
  for(i in 1:length(mcmc.use.mod)){
    out[i,,]<-crossprod(diag(c(bmars$beta[mcmc.use.mod[i],1:(M+1)]),M+1)%*%tmat.des,tmat.func)
  }

  return(out)
}
