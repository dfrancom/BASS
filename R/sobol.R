########################################################################
## functions for Sobol decomposition - these all use scaling from const function
########################################################################
pCoef<-function(i,q){ # refer to paper
  factorial(q)^2*(-1)^i/(factorial(q-i)*factorial(q+1+i))
}
intabq<-function(a,b,t1,t2,q){
  #integral from a to b of [(x-t1)(x-t2)]^q when q positive integer
  sum(pCoef(0:q,q)*(b-t1)^(q-0:q)*(b-t2)^(q+1+0:q)) - sum(pCoef(0:q,q)*(a-t1)^(q-0:q)*(a-t2)^(q+1+0:q))
}
C2<-function(k,m,n,tl){ # integral of two pieces of tensor that have same variable - deals with sign, truncation
  q<-tl$q
  t1<-tl$t[n,k]
  s1<-tl$s[n,k]
  t2<-tl$t[m,k]
  s2<-tl$s[m,k]
  cc<-const(signs=c(s1,s2),knots=c(t1,t2),degree=q)
  if((s1*s2)==0){
    return(0)
  }
  if(t2<t1){
    t1<-tl$t[m,k]
    s1<-tl$s[m,k]
    t2<-tl$t[n,k]
    s2<-tl$s[n,k]
  }
  if(m==n){ #t1=t2, s1=s2
    return(1/(2*q+1)*((s1+1)/2-s1*t1)^(2*q+1)/cc)
  } else{
    if(s1==1){
      if(s2==1){
        return(intabq(t2,1,t1,t2,q)/cc)
      } else{
        return(intabq(t1,t2,t1,t2,q)*(-1)^q/cc)
      }
    } else{
      if(s2==1){
        return(0)
      } else{
        return(intabq(0,t1,t1,t2,q)/cc)
      }
    }
  }
}

# Vu<-function(u,tl){ # sobol main effect variances - where most of the time is spent
#   CCu<-apply(tl$C1.all.prod3[,,u,drop=F],1:2,prod) # TODO: utilize symmetry
#   C2.temp<-apply(tl$C2.all2[,,u,drop=F],1:2,prod)
#   mat<-tl$CC*(C2.temp/CCu-1)
#   return(apply(tl$a,1,function(x) t(x)%*%mat%*%x))
# }

VuMat<-function(u,tl){ # sobol main effect variances - where most of the time is spent
  CCu<-apply(tl$C1.all.prod3[,,u,drop=F],1:2,prod) # TODO: utilize symmetry
  C2.temp<-apply(tl$C2.all2[,,u,drop=F],1:2,prod)
  mat<-tl$CC*(C2.temp/CCu-1)
  #browser()
  return(mat)
}

Vu<-function(u,tl){ # sobol main effect variances - where most of the time is spent
  mat<-VuMat(u,tl)
  out<-apply(tl$a,1,function(x) t(x)%*%mat%*%x)
  # if(any(out<0))
  #   browser()
  return(out)
}

Vu_des_func<-function(u,tl){ # sobol main effect variances - where most of the time is spent
  #browser() # something wrong if length of u is more than 1
  mat<-VuMat(u,tl)
  nx<-length(tl$xx)
  nmodels<-length(tl$a[,1])
  out<-matrix(nrow=nmodels,ncol=nx)
  for(i in 1:nmodels){
    for(j in 1:nx){
      tt<-tl$a[i,]*tl$tfunc.basis[,j]
      out[i,j]<-t(tt)%*%mat%*%tt
    }
  }
  return(out)#-tl$f0.sq)
  #apply(func.basis,1,function(r) t(r)%*%mat%*%(r)) # if func.basis is nx by nbasis
  #return(apply(tl$a,1,function(x) t(x)%*%mat%*%x)) # tl$a is nbasis by nmodels
}

VuInt<-function(u,tl){ # sobol interaction variances
  add<-0
  len<-length(u)
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # sum(cs.num.ind[l-1]) makes it 0 when it should be # this gets index for which combs are subsets of u
    add<-add+(-1)^(len-l)*rowSums(tl$temp[,ind,drop=F])
  }
  add[abs(add)<1e-15]<-0
  if(any(add<0))
    browser()
  return(add)
}

VuInt_des_func<-function(u,tl){ # sobol interaction variances
  add<-0
  len<-length(u)
  #browser()
  for(l in 1:len){
    ind<-((sum(tl$cs.num.ind[l-1])+1):tl$cs.num.ind[l])[apply(tl$combs[[l]],2,function(x) all(x%in%u))] # sum(cs.num.ind[l-1]) makes it 0 when it should be # this gets index for which combs are subsets of u
    add<-add+(-1)^(len-l)*apply(tl$temp[,ind,,drop=F],c(1,3),sum)#colSums(tl$temp[,ind,,drop=F],dims=2) # something like this
  }
  return(add)
}



#' @title BASS Sensitivity Analysis
#'
#' @description Decomposes the variance of the BASS model into variance due to main effects, two way interactions, and so on, similar to the ANOVA decomposition for linear models.  Uses the Sobol' decomposition, which can be done analytically for MARS models.
#' @param mod a fitted model output from the \code{bass} function.
#' @param mcmc.use an integer vector indexing which MCMC iterations to use for sensitivity analysis.
#' @param func.var an integer indicating which functional variable to make sensitivity indices a function of.  Disregard if \code{mod} is non-functional or if scalar sensitivity indices are desired.
#' @param xx.func.var grid for functional variable specified by \code{func.var}.  Disregard if \code{func.var} is not specified.  If \code{func.var} is specified and \code{xx.func.var} not specified, the grid used to fit \code{mod} will be used.
#' @param verbose logical; should progress be displayed?
#' @details Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
#' @return If non-functional (\code{func.var = NULL}), a list with two elements:
#'  \item{S}{a data frame of sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular main effect or interaction.  The values are the proportion of variance in the model that is due to each main effect or interaction.}
#'  \item{T}{a data frame of total sensitivity indices with number of rows matching the length of \code{mcmc.use}.  The columns are named with a particular variable.}
#'  Otherwise, a list with four elements:
#'  \item{S}{an array with first dimension corresponding to MCMC samples (same length as \code{mcmc.use}), second dimension corresponding to different main effects and interactions (labeled in \code{names.ind}), and third dimension corresponding to the grid used for the functional variable.  The elements of the array are sensitivity indices.}
#'  \item{S.var}{same as \code{S}, but scaled in terms of total variance rather than percent of variance.}
#'  \item{names.ind}{a vector of names of the main effects and interactions used.}
#'  \item{xx}{the grid used for the functional variable.}
#'
#' @keywords BMARS
#' @seealso \link{bass} for model fitting and \link{predict.bass} for prediction.
#' @export
#' @examples
#' # See examples in bass documentation.
#'
sobol<-function(mod,mcmc.use=NULL,func.var=NULL,xx.func.var=NULL,verbose=TRUE){ # note: requires inputs to be scaled to [0,1]
  if(mod$p==1 & !mod$func)
    stop('Sobol only used for multiple input models')
  mcmc.use.poss<-1:((mod$nmcmc-mod$nburn)/mod$thin)
  if(any(!(mcmc.use%in%mcmc.use.poss))){
    mcmc.use<-mcmc.use.poss
    warning('disregarding mcmc.use because of bad values')
  }
  if(any(is.null(mcmc.use))){
    mcmc.use<-mcmc.use.poss
  }

  if(is.null(func.var)){
    func<-F
  } else{
    func<-T
    if(!mod$func){
      func<-F
      warning('disregarding func.var because mod parameter is not functional')
    }
  }



  if(func){
    #if(is.null(func.var))
    #  func.var<-1
    if(!(func.var%in%(1:ncol(mod$xx.func))))
      stop('func.var in wrong range of values')
    if(is.null(xx.func.var)){
      xx.func.var<-mod$xx.func[,func.var,drop=F]
    } else{
      #if(dim(xx.func.var)[2])
      rr<-range(xx.func.var)
      if(rr[1]<mod$range.func[1,func.var] | rr[2]>mod$range.func[2,func.var])
        warning(paste('range of func.var in bass function (',mod$range.func[1,func.var],',',mod$range.func[2,func.var],') is smaller than range of xx.func.var (',rr[1],',',rr[2],'), indicating some extrapolation',sep=''))
      xx.func.var<-scale.range(xx.func.var,mod$range.func[,func.var])
    }

    return(sobol_des_func(mod=mod,mcmc.use=mcmc.use,verbose=verbose,func.var=func.var,xx.func.var=xx.func.var))
  } else{
    return(sobol_des(mod=mod,mcmc.use=mcmc.use,verbose=verbose)) # applies to both des & func as long as functional sobol indices are not desired
  }
}











getCombs<-function(mod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var=NULL){
  vf<-mod$vars.func[uniq.models,,]
  if(!is.null(func.var))
    vf[vf==func.var]<-NA
  n.un<-array(c(as.integer(mod$vars.des[uniq.models,,]),as.integer(vf+mod$pdes)),dim=c(nmodels,maxBasis,maxInt.tot))
  n.un<-apply(n.un,1:2,sort)
  n.un<-unique(c(n.un))
  n.un[sapply(n.un,length)==0]<-NULL
  temp<-list()
  for(ii in 1:length(n.un)){
    pp<-length(n.un[[ii]])
    if(pp==1){
      temp<-c(temp,n.un[[ii]])
    } else{
      temp<-c(temp,do.call(c,sapply(1:pp,function(x) combn(n.un[[ii]],x,simplify=F))))
    }
  }
  temp<-lapply(temp,as.integer)
  n.un<-unique(c(n.un,unique(temp)))
  ord<-order(sapply(n.un,length))
  n.un<-n.un[ord]
  a<-NA
  for(ii in 1:maxInt.tot)
    a[ii]<-which(sapply(n.un,length)==ii)[1]
  a[maxInt.tot+1]<-length(n.un)+1

  combs<-names.ind<-ll<-mat<-list()
  aa<-a
  a<-na.omit(a)
  k<-0
  for(ii in (1:maxInt.tot)[!is.na(aa[-(maxInt.tot+1)])]){
    k<-k+1
    if(!is.na(a[ii])){
      mat<-do.call(rbind,n.un[a[k]:(a[k+1]-1)])
      mat<-mat[do.call(order, as.data.frame(mat)),,drop=F]
      combs[[ii]]<-t(mat)
      names.ind[[ii]]<-apply(combs[[ii]],2,paste,collapse='x') # labels for output
      ll[[ii]]<-split(combs[[ii]],rep(1:ncol(combs[[ii]]),each=nrow(combs[[ii]])))
    }
  }
  #browser()
  num.ind<-sapply(ll,length)
  cs.num.ind<-cumsum(num.ind) # used for indexing
  return(list(combs=combs,names.ind=names.ind,ll=ll,num.ind=num.ind,cs.num.ind=cs.num.ind,aa=aa))
}





get_tl<-function(mod,mcmc.use.m,M,m,p,q,cs.num.ind,combs,func.var=NULL,xx.func.var=NULL){
  a<-mod$beta[mcmc.use.m,2:(M+1),drop=F] # basis coefficients excluding intercept
  vf<-mod$vars.func[m,1:M,]
  #browser()
  if(!is.null(func.var)){
    vf[vf==func.var]<-NA
    #vf[which(vf>func.var,arr.ind = T)]<-vf[which(vf>func.var,arr.ind = T)]-1
  }
  Kind<-cbind(mod$vars.des[m,1:M,],vf+mod$pdes)
  #browser()
  if(M==1){
    Kind<-t(Kind)
  }
  t<-s<-matrix(0,nrow=M,ncol=p)
  for(k in 1:M){ # these matrices mimic the output of earth
    n.int.des<-mod$n.int.des[m,k]
    knotInd.des<-mod$knotInd.des[m,k,1:n.int.des]
    vars.des<-mod$vars.des[m,k,1:n.int.des]
    t[k,vars.des]<-mod$xx.des[cbind(knotInd.des,vars.des)]
    s[k,vars.des]<-mod$signs.des[m,k,1:n.int.des]
    if(mod$func){
      n.int.func<-mod$n.int.func[m,k]
      knotInd.func<-mod$knotInd.func[m,k,1:n.int.func]
      vars.func<-mod$vars.func[m,k,1:n.int.func]
      t[k,vars.func+mod$pdes]<-mod$xx.func[cbind(knotInd.func,vars.func)]
      s[k,vars.func+mod$pdes]<-mod$signs.func[m,k,1:n.int.func]
    }
  }
  #ind<-1:p
  if(!is.null(func.var)){
    #ind<-ind[-(mod$pdes+func.var)]
    s[,mod$pdes+func.var]<-t[,mod$pdes+func.var]<-0
  }
    
  tl<-list(s=s,t=t,q=q,a=a,M=M,Kind=Kind,cs.num.ind=cs.num.ind,combs=combs,xx=xx.func.var) #temporary list
  #tl<-list(s=s[,ind,drop=F],t=t[,ind,drop=F],q=q,a=a,M=M,Kind=Kind,cs.num.ind=cs.num.ind,combs=combs,xx=xx.func.var)
  return(tl)
}

add_tl<-function(tl,p){
  # TODO: these need better names
  C1.all<-(1/(tl$q+1)*((tl$s+1)/2-tl$s*tl$t))*tl$s^2 # so I don't need C function anymore
  C1.all2<-replace(C1.all,which(C1.all==0,arr.ind=T),1) # for products, make 0's 1's
  C1.all.prod<-apply(C1.all2,1,prod)
  tl$CC<-tcrossprod(C1.all.prod)
  C2.all<-C1.all.prod2<-both.ind<-array(0,dim=c(tl$M,tl$M,p))
  # TODO: probably more efficient way of storing since symmetric
  for(ii in 1:tl$M){
    for(jj in ii:tl$M){
      bb<-intersect(na.omit(tl$Kind[ii,]),na.omit(tl$Kind[jj,])) # variables that basis functions ii and jj have in common
      if(length(bb)>0){
        both.ind[ii,jj,bb]<-both.ind[jj,ii,bb]<-1
        C2.all[ii,jj,]<-C2.all[jj,ii,]<-apply(t(1:p),2,C2,m=ii,n=jj,tl=tl)
        C1.all.prod2[ii,jj,]<-C1.all.prod2[jj,ii,]<-C1.all[ii,]*C1.all[jj,] # pairwise products of C1.all
      }
    }
  }
  #browser()
  tl$C1.all.prod3<-C1.all.prod2
  tl$C1.all.prod3[C1.all.prod2==0]<-1
  tl$C2.all2<-C2.all
  tl$C2.all2[!as.logical(both.ind)]<-1
  return(tl)
}





getTot<-function(ll,sob,names.ind,p,maxInt.tot,aa){
  #browser()
  vars.use<-unique(unlist(ll))
  puse<-length(vars.use)
  ll[[1]]<-numeric(0) # for proper lengths below
  tot<-sob[,1:length(names.ind[[1]])]#matrix(0,nrow=nrow(sob),ncol=p)
  #tot[,as.numeric(names.ind[[1]])]<-sob[,as.numeric(names.ind[[1]])] # get total sensitivity indices
  if(maxInt.tot>1){
    for(pp in 1:puse){
      for(l in (2:maxInt.tot)[!is.na(aa[-c(1,maxInt.tot+1)])]){
        tot[,pp]<-tot[,pp]+rowSums(
          sob[,puse+length(ll[[l-1]])+which(apply(do.call(rbind,ll[[l]]),1,function(r){vars.use[pp]%in%r})),drop=F]
        )
      }
    }
  }
  return(tot)
}



sobol_des<-function(mod,mcmc.use,verbose){
  models<-mod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-mod$maxInt.des
  maxBasis<-dim(mod$vars.des)[2]
  q<-mod$degree
  i<-1
  p<-mod$pdes

  if(mod$func){
    p<-p+mod$pfunc
    maxInt.tot<-maxInt.tot+mod$maxInt.func
  }

  ################################################
  # get combs & ll including functional variables
  ################################################
  tt<-getCombs(mod,uniq.models,nmodels,maxBasis,maxInt.tot)
  combs<-tt$combs
  names.ind<-tt$names.ind
  ll<-tt$ll
  num.ind<-tt$num.ind
  cs.num.ind<-tt$cs.num.ind
  ################################################
  sob<-array(0,dim=c(length(mcmc.use),sum(num.ind)))#,length(xx.func.var)))

  if(verbose)
    cat('Sobol Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')

  mod.ind<-0
  for(m in uniq.models){ #do this in parallel?
    mod.ind<-mod.ind+1
    mcmc.use.m<-mcmc.use[models==m] # which part of mcmc.use does this correspond to?
    mod.m.ind<-i:(i+length(mcmc.use.m)-1)
    M<-mod$nbasis[mcmc.use.m][1] # number of basis functions in this model
    if(M>0){
      lens<-mod$n.int.des[m,1:M]
      if(mod$func)
        lens<-lens+mod$n.int.func[m,1:M]
      tl<-get_tl(mod,mcmc.use.m,M,m,p,q,cs.num.ind,combs)
      tl<-add_tl(tl,p)


      var.tot<-Vu(1:p,tl) # total variance
      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used?
      vars.used<-sort(vars.used)

      tl$temp<-matrix(0,nrow=length(mcmc.use.m),ncol=max(cs.num.ind)) # where we store all the integrals (not normalized by subtraction) - matches dim of sob
      tl$temp[,which(combs[[1]]%in%vars.used)]<-apply(t(vars.used),2,Vu,tl=tl)
      sob[mod.m.ind,1:cs.num.ind[1]]<-tl$temp[,1:cs.num.ind[1]]

      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$temp is made sequentially)
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l]
          mm<-matrix(nrow=M,ncol=length(ll[[l]]))
          for(o in 1:M){
            mm[o,]<-unlist(lapply(ll[[l]],function(el){prod(el%in%tl$Kind[o,])})) #all the variables in question must be in the basis
          }
          use<-colSums(mm)>0
          tl$temp[,int.l.ind[use]]<-apply(combs[[l]][,use,drop=F],2,Vu,tl=tl) # perform the necessary integration
          cat<-matrix(0,nrow=length(mcmc.use.m),ncol=length(ll[[l]]))
          for(l2 in (1:length(ll[[l]]))[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)
            cat[,l2]<-VuInt(ll[[l]][[l2]],tl) # do the normalizing
          }
          sob[mod.m.ind,int.l.ind]<-cat # dividing matrix by a vector
        }
      }
      sob[mod.m.ind,]<-sob[mod.m.ind,]/var.tot # sobol indices
      i<-i+length(mcmc.use.m) # for index


      if(verbose & mod.ind%%10==0)
        cat('Sobol',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.ind,'\n')


    }
  }

  sob<-as.data.frame(sob)
  names(sob)<-unlist(names.ind) # give labels

  if(verbose)
    cat('Total Sensitivity',timestamp(prefix='#--',suffix='--#',quiet=T),'\n')

  tot<-getTot(ll,sob,names.ind,p,maxInt.tot,tt$aa)
  ret<-list(S=sob,T=tot,func=F)#,var.tot=var.tot))
  class(ret)<-'bassSob'
  return(ret)
}





makeBasisMatrixVar<-function(i,nbasis,vars,signs,knots.ind,q,xxt,n.int,xx.train,var){ # make for only one variable
  n<-ncol(xxt)
  tbasis.mat<-matrix(nrow=nbasis+1,ncol=n)
  tbasis.mat[1,]<-1
  if(nbasis>0){
    for(m in 1:nbasis){
      if(all(na.omit(vars[i,m,])!=var)){
        tbasis.mat[m+1,]<-1 # could do this at beginning
      } else{
        use<-which(vars[i,m,]==var)#1:n.int[i,m]
        knots<-xx.train[cbind(knots.ind[i,m,use],vars[i,m,use])] # get knots from knots.ind
        tbasis.mat[m+1,]<-makeBasis(signs[i,m,use],1,knots,xxt,q)
      }
    }
  }
  return(tbasis.mat)
}


sobol_des_func<-function(mod,mcmc.use,verbose,func.var,xx.func.var){
  models<-mod$model.lookup[mcmc.use] # only do the heavy lifting once for each model
  uniq.models<-unique(models)
  nmodels<-length(uniq.models)
  maxInt.tot<-mod$maxInt.des
  maxBasis<-dim(mod$vars.des)[2]
  q<-mod$degree
  i<-1
  p<-mod$pdes

  if(mod$func){
    p<-p+mod$pfunc
    maxInt.tot<-maxInt.tot+mod$maxInt.func
  }

  ################################################
  # get combs & ll including functional variables
  ################################################
  tt<-getCombs(mod,uniq.models,nmodels,maxBasis,maxInt.tot,func.var)
  combs<-tt$combs
  names.ind<-tt$names.ind
  ll<-tt$ll
  num.ind<-tt$num.ind
  cs.num.ind<-tt$cs.num.ind
  #browser()
  ################################################
  sob<-sob2<-array(0,dim=c(length(mcmc.use),sum(num.ind),length(xx.func.var)))

  if(verbose)
    cat('Sobol Start',timestamp(prefix='#--',suffix='--#',quiet=T),'Models:',length(unique(models)),'\n')

  mod.ind<-0
  for(m in uniq.models){ #do this in parallel?
    mod.ind<-mod.ind+1
    mcmc.use.m<-mcmc.use[models==m] # which part of mcmc.use does this correspond to?
    mod.m.ind<-i:(i+length(mcmc.use.m)-1)
    M<-mod$nbasis[mcmc.use.m][1] # number of basis functions in this model
    if(M>0){
      # lens<-mod$n.int.des[m,1:M]
      # if(mod$func){
      #   nint<-mod$n.int.func[m,1:M]
      #   ii<-which(mod$vars.func[m,1:M,,drop=F]==func.var,arr.ind = T)[1,]
      #   nint[ii]<-nint[ii]-1
      #   lens<-lens+nint
      # }

      tl<-get_tl(mod,mcmc.use.m,M,m,p,q,cs.num.ind,combs,func.var,xx.func.var)
      tl<-add_tl(tl,p)
      lens<-apply(tl$Kind,1,function(x) length(na.omit(x)))

      tl$tfunc.basis<-makeBasisMatrixVar(m,M,vars=mod$vars.func,signs=mod$signs.func,knots.ind=mod$knotInd.func,q=mod$degree,xxt=t(tl$xx),n.int=mod$n.int.func,xx.train=mod$xx.func,var=func.var)[-1,]

      #browser()
      var.tot<-Vu_des_func(1:p,tl) # total variance
      vars.used<-unique(unlist(na.omit(c(tl$Kind)))) # which variables are used?
      vars.used<-sort(vars.used)
      
      tl$temp<-array(0,dim=c(length(mcmc.use.m),max(cs.num.ind),length(xx.func.var))) # where we store all the integrals (not normalized by subtraction)
      jj=0
      for(pp in vars.used){
        jj=jj+1
        tl$temp[,jj,]<-Vu_des_func(pp,tl)#t(apply(t(vars.used),2,Vu_des_func,tl=tl))
      }
      #browser()
      sob[mod.m.ind,1:cs.num.ind[1],]<-tl$temp[,1:cs.num.ind[1],]
      
      #tl$temp<-matrix(0,nrow=length(mcmc.use.m),ncol=max(cs.num.ind)) # where we store all the integrals (not normalized by subtraction) - matches dim of sob
      #tl$temp[,which(combs[[1]]%in%vars.used)]<-apply(t(vars.used),2,Vu,tl=tl)
      #sob[mod.m.ind,1:cs.num.ind[1]]<-tl$temp[,1:cs.num.ind[1]]

#browser()
      if(max(lens)>1){ # if there are any interactions
        for(l in 2:max(lens)){ # must go in order for it to work (tl$temp is made sequentially)
          #browser()
          int.l.ind<-(cs.num.ind[l-1]+1):cs.num.ind[l]
          mm<-matrix(nrow=M,ncol=length(ll[[l]]))
          for(o in 1:M){
            mm[o,]<-unlist(lapply(ll[[l]],function(el){prod(el%in%tl$Kind[o,])})) #all the variables in question must be in the basis
          }
          use<-colSums(mm)>0
          for(pp in which(use)){
            tl$temp[,int.l.ind[pp],]<-Vu_des_func(combs[[l]][,pp,drop=F],tl)#apply(combs[[l]][,pp,drop=F],2,Vu_des_func,tl=tl) # perform the necessary integration
          }
          cat<-array(0,dim=c(length(mcmc.use.m),length(ll[[l]]),length(xx.func.var)))

          for(l2 in (1:length(ll[[l]]))[use]){ # only go through the interactions that are actually in Kind (but still allow for 2-way when actual is 3-way, etc.)

            cat[,l2,]<-VuInt_des_func(ll[[l]][[l2]],tl) # do the normalizing
          }
          sob[mod.m.ind,int.l.ind,]<-cat # dividing matrix by a vector
        }
      }
      #browser()
      kk<-0
      for(ii in mod.m.ind){
        kk=kk+1
        sob2[ii,,]<-t(t(sob[ii,,])/var.tot[kk,]) # sobol indices
      }
      i<-i+length(mcmc.use.m) # for index


      if(verbose & mod.ind%%10==0)
        cat('Sobol',timestamp(prefix='#--',suffix='--#',quiet=T),'Model:',mod.ind,'\n')


    }
  }

  #sob<-as.data.frame(sob)
  #names(sob)<-unlist(names.ind) # give labels

  ret<-list(S=sob2,S.var=sob,names.ind=unlist(names.ind),xx=tl$xx,func=T)#,var.tot=var.tot))
  class(ret)<-'bassSob'
  return(ret)
}











