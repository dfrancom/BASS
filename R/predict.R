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
} # I think prediction on a test set should be kept separate from the bass function, simplifies tempering (only would want to predict on the cold chain).

makeBasisMatrixCat<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,Xt){
  return(1)
}

#' @title BASS Prediction
#'
#' @description Predict function for BASS.  Outputs the posterior predictive samples based on the specified MCMC iterations.
#' @param object a fitted model, output from the \code{bass} function.
#' @param newdata a matrix of new input values at which to predict.  The columns should correspond to the same variables used in the \code{bass} function.
#' @param newdata.func a matrix of new values of the functional variable.  If none, the same values will be used as in the training data.
#' @param mcmc.use a vector indexing which MCMC iterations to use for prediction.
#' @param verbose logical; should progress be displayed?
#' @details Efficiently predicts when two MCMC iterations have the same basis functions (but different weights).
#' @return If model output is a scalar, this returns a matrix with the same number of rows as \code{newdata} and columns corresponding to the the MCMC iterations \code{mcmc.use}.  These are samples from the posterior predictive distribution.  If model output is functional, this returns an array with first dimension corresponding to MCMC iteration, second dimension corresponding to the rows of \code{newdata}, and third dimension corresponding to the rows of \code{newdata.func}.
#' @keywords BMARS
#' @seealso \link{bass} for model fitting and \link{sobol} for sensitivity analysis.
#' @export
#' @examples
#' # See examples in bass documentation.
#'
predict.bass<-function(object,newdata,newdata.func=NULL,mcmc.use=NULL,verbose=FALSE){
  if(is.null(mcmc.use)){ # if null, use all
    mcmc.use<-1:((object$nmcmc-object$nburn)/object$thin)
  }
  if(object$func){
    if(is.null(newdata.func))
      newdata.func<-object$xx.func
    else{
      dxf<-dim(newdata.func)
      if(is.null(dxf))
        newdata.func<-matrix(newdata.func)
      for(i in 1:ncol(newdata.func)){
        newdata.func[,i]<-scale.range(newdata.func[,i],object$range.func[,i])
      }
    }
  } else{
    newdata.func<-t(1) # placeholder
  }

  newdata<-as.data.frame(newdata)

  for(i in 1:ncol(newdata)){
   newdata[,i]<-scale.range(newdata[,i],object$range.des[,i])
  }
  out<-array(dim=c(length(mcmc.use),nrow(newdata),nrow(newdata.func)))
  k<-0
  models<-object$model.lookup[mcmc.use]
  if(verbose)
    cat('Predict Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')
  mod.ind<-0
  for(j in unique(models)){ # loop though models, could be parallel?
    mod.ind<-mod.ind+1
    mcmc.use.j<-mcmc.use[models==j]
    ind<-k+(1:length(mcmc.use.j)) # index for storage
    k<-k+length(ind) # used for start of index
    out[ind,,]<-eval(parse(text=paste('mult',object$type,sep='')))(model=j,mcmc.use.mod=mcmc.use.j,object=object,newdata=newdata,newdata.func=newdata.func)
    if(verbose & mod.ind%%100==0)
      cat('Predict',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.ind,'\n')
  }
  return(drop(out))
}

mult_des<-function(model,mcmc.use.mod,object,newdata,newdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat<-makeBasisMatrix(model,M,object$vars,object$signs,object$knotInd,object$degree,t(newdata),object$n.int,object$xx.des) ## need to get rid of these transposes
  out<-object$beta[mcmc.use.mod,1:(M+1),drop=F]%*%tmat
  return(out)
}

mult_des_func<-function(model,mcmc.use.mod,object,newdata,newdata.func){
  M<-object$nbasis[mcmc.use.mod[1]]
  tmat.des<-makeBasisMatrix(model,M,object$vars.des,object$signs.des,object$knotInd.des,object$degree,t(newdata),object$n.int.des,object$xx.des)
  tmat.func<-makeBasisMatrix(model,M,object$vars.func,object$signs.func,object$knotInd.func,object$degree,t(newdata.func),object$n.int.func,object$xx.func)
  out<-array(dim=c(length(mcmc.use.mod),nrow(newdata),nrow(newdata.func)))
  for(i in 1:length(mcmc.use.mod)){
    out[i,,]<-crossprod(diag(c(object$beta[mcmc.use.mod[i],1:(M+1)]),M+1)%*%tmat.des,tmat.func)
  }

  return(out)
}
