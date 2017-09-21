## Author: Carl F. Falk
## First version assumed no access to internal lavaan functions.
## A re-write could make this more efficient.
##
## General purpose likelihood-based confidence interval function using lavaan
## Requires that the parameter or quantity of interest is clearly labeled in the parameter table 
lci<-function(object, label, level=.95, bound=c("lower","upper"),
              optimizer="Rsolnp",
              ci.method="NealeMiller1997",
              start=NULL,diff.method="default",Dtol=.05,
              iterlim=50,reoptimize=FALSE){

  ## input checking
  if(!is.null(start)){
    stop("Custom starting values not yet supported")
  }
  if(reoptimize){
    stop("Re-optimization is not yet supported")
  }
  if(class(object)!="lavaan"){
    stop("object must be a fitted lavaan model")
  }
  if(level <=0 | level >=1){
    stop("level must be between 0 and 1")
  }
  if(ci.method=="bisect" & object@Options$se=="none"){
    stop("Bisection method currently relies in part on standard errors to aid in determining where to start the algorithm. Change se='none' to something else.")
  }
  
  ## extract parameter table
  ptable<-parTable(object)

  ## Look for parameter label
  pindx<-which(ptable$label==label)
  
  if(length(pindx)<1){
    stop("Parameter label not found in lavaan parameter table.")
  } else if(length(pindx)>1){
    warning("FIXME: More than one parameter label match")
  }
  
  ## point estimate based on fitted model
  est<-ptable$est[pindx[1]]
  #ptable$ustart<-ptable$start<-ptable$est # not sure if necessary
  if(is.null(start)){
    start<-est
  }
  
  ## Start list of output before embarking on optimization
  result<-list()
  result$est<-est
    
  crit<-qchisq(level,1) ## hardcoded 1 df for now
  result$crit<-crit
  
  ## Check whether test is satorra.2000
  ## Code borrowed directly from lavTestLRT function
  estimator<-object@Options$estimator
  mods.scaled <- unlist(lapply(list(object), function(x) {
    any(c("satorra.bentler", "yuan.bentler", "mean.var.adjusted", 
          "scaled.shifted") %in% unlist(sapply(slot(x, "test"), 
                                               "[", "test")))
  }))
  if(all(mods.scaled)){
    scaled <- TRUE
    TEST <- object@test[[2]]$test
  } else {
    scaled<-FALSE
  }
  if(scaled){
    if (diff.method[1] == "default") {
      if (estimator == "PML") {
        diff.method[1] <- "mean.var.adjusted.PLRT"
      }
      else if (TEST %in% c("satorra.bentler", "yuan.bentler")) {
        diff.method[1] <- "satorra.bentler.2001"
      }
      else {
        diff.method[1] <- "satorra.2000"
      }
    }
    else if (diff.method[1] == "meanvaradjustedplrt") {
      diff.method[1] <- "mean.var.adjusted.PLRT"
      stopifnot(estimator == "PML")
    }
    else if (diff.method[1] == "satorra.2000") {
      diff.method[1]<- "satorra.2000"
    }
    else if (diff.method[1] == "satorra.bentler.2001") {
      diff.method[1]<- "satorra.bentler.2001"
    }
    else if (diff.method[1]== "satorra.bentler.2010") {
      diff.method[1]<- "satorra.bentler.2010"
    }
    else {
      stop("lavaan ERROR: unknown method for scaled difference test: ", 
           diff.method)
    }
  }
  
  ## If satorra.2000, obtain scaling constant before optimizing and adjust crit accordingly
  ## Then, estimate using normal theory ML
  if(diff.method[1]=="satorra.2000"){
    const<-list()
    const[[label]]<-est
    M0<-lci_refit(object,const)
    # try again if not converged; apparently just rounding est can work
    if(!M0@Fit@converged){
      const[[label]]<-round(est,2)
      M0<-lci_refit(object,const)
    }

    chat<-lav_test_diff_Satorra2000(object, M0, A.method="exact")$scaling.factor
    crit<-chat*crit

    object<-lci_refit(object,estimator="ML")

  } else {
    object<-lci_refit(object)
  }
  
  ## Begin optimization
  for(b in bound){
    D<-NULL
    cont<-NULL
    est.bound<-NULL
    if(ci.method[1]=="NealeMiller1997"){
      fitfunc<-lci_nealemiller1997

      if(ptable$op[pindx]==":=" & (diff.method[1]=="satorra.2000"|estimator=="ML")){
        label<-all.vars(parse(file="", text=ptable$rhs[pindx]))
        if(length(start)!=length(label)){
          start<-ptable$est[match(label,ptable$label)]
        }
        diff.method<-"default"
      }
      
      if(optimizer[1]=="Rsolnp"){
        LCI<-try(Rsolnp::solnp(start,fitfunc,fitmodel=object,
                      label=label, pindx=pindx, crit=crit, bound=b,
                      diff.method=diff.method,control=list(trace=0)))

        if(class(LCI)!="try-error"){
          D<-lci_diff_test(LCI$pars,object,label,diff.method=diff.method)
          est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
          attr(D,"mod")<-NULL
          conv<-LCI$convergence
        } else{
          est.bound<-NA
        }
      } else if (optimizer[1]=="optim"){
        LCI<-try(optim(start, fitfunc,fitmodel=object,
                      label=label, pindx=pindx, crit=crit, bound=b,
                      diff.method=diff.method,method="BFGS"))
        if(class(LCI)!="try-error"){
          D<-lci_diff_test(LCI$par,object,label,diff.method=diff.method)
          est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
          attr(D,"mod")<-NULL          
          conv<-LCI$convergence
        } else{
          est.bound<-NA
        }
      } else if (optimizer[1]=="nlminb"){
        LCI<-try(nlminb(start, fitfunc,fitmodel=object,
                      label=label, pindx=pindx, crit=crit, bound=b,
                      diff.method=diff.method))
        if(class(LCI)!="try-error"){
          D<-lci_diff_test(LCI$par,object,label,diff.method=diff.method)
          est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
          attr(D,"mod")<-NULL
          conv<-LCI$convergence
        } else{
          p.est.bound<-NA
        }
      } else {
        stop("Unsupported optimizer")
      }
    } else if (ci.method[1]=="WuNeale2012" & optimizer[1]=="Rsolnp"){
      stop("Wu and Neale is not yet supported")
      eqfun<-lci_diff_test

      if(class(LCI)!="try-error"){
        D<-lci_diff_test(LCI$pars,object,label,diff.method=diff.method)
        est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
        attr(D,"mod")<-NULL
        conv<-LCI$convergence
      } else{
        est.bound<-NA
      }
      
    } else if (ci.method[1]=="bisect"){
      set.seed(1234)
      if(diff.method[1]=="satorra.2000"){
        diff.method<-"default"
      }
      LCI<-lci_bisect(object,label,crit, bound=b,diff.method=diff.method[1],iterlim=iterlim)
      if(class(LCI)!="try-error"){
        D<-lci_diff_test(LCI$est,object,label,diff.method=diff.method)
        attr(D,"mod")<-NULL        
        est.bound<-LCI$est
        conv<-LCI$iter
      } else{
        est.bound<-NA
      }
    } else {
      stop("Unsupported ci.type and/or optimizer")
    }
    
    if(b=="upper"){
      result$upper<-est.bound
      result$convupper<-conv
      result$Dupper<-(result$crit/crit)*D
    }
    if(b=="lower"){
      result$lower<-est.bound
      result$convlower<-conv
      result$Dlower<-(result$crit/crit)*D
    }
  }
  
  return(result)
}

# Just computes constrains quantity in label to p and computes difference test for use with lci functions
lci_diff_test<-function(p,fitmodel,label,diff.method="default"){

  const<-list()
  for(i in 1:length(p)){
    const[[label[i]]]<-p[i]
  }
  M0<-lci_refit(fitmodel,const)

  if(class(M0)=="try-error"){
    return(NA)
  } else {
    if(M0@Fit@converged){
      D<-try(lavTestLRT(fitmodel,M0,method=diff.method)$`Chisq diff`[2])
      if(class(D)!="try-error"){
        attr(D,"mod")<-M0
        return(D)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  }
}

## Optimization for function of model parameters using a modified version of Neale and Miller (1997)
## Differs from original in that the function only takes input for the parmeter (or function of paramters) of interest
## and allows lavaan to fully optimize with each function evaluation
lci_nealemiller1997<-function(p, fitmodel, label, pindx, crit, bound=c("lower","upper"),
                              diff.method="default"){
  
  D<-lci_diff_test(p,fitmodel,label,diff.method)
  if(!is.na(D)){
    fit<-(D-crit)^2
  
    # just in case this is a function of parameters
    p<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
  
    if(bound=="lower"){
      fit<-fit+p
    } else if (bound=="upper"){
      fit<-fit-p
    } else {
      stop("No valid bound specified")
    }
    attr(fit,"D")<-D
  } else {
    fit<-NA
  }
  return(fit)
}

lci_bisect<-function(fitmodel,label,crit,tol=1e-5,iterlim=25,init=2,bound=c("lower","upper"),
                     diff.method="default",lb=-Inf,ub=Inf){
  ## FIXME: negative init would break things
  
  ## extract parameter table
  ptable<-parTable(fitmodel)
  
  ## Look for parameter label
  pindx<-which(ptable$label==label)
  
  if(length(pindx)<1){
    stop("Parameter label not found in lavaan parameter table.")
  }
  
  ## point estimate based on fitted model
  est<-ptable$est[pindx[1]]
  if(is.null(start)){
    start<-est
  }
  
  ## standard error
  se<-ptable$se[pindx[1]]
  
  ## First part - find out two points - between which is the boundary
  # starting points are function of point estimate and multiplier w/ se  
  p0<-est
  if(bound=="lower"){
    sign<- -1
  } else {
    sign<-1
  }
  p2<-est+sign*init*se
  
  ## enforce hard boundaries
  if(!is.null(lb)){
    p2[p2<lb]<-lb
  }
  if(!is.null(ub)){
    p2[p2>ub]<-ub
  }

  ## Try to find some value of the quantity of interest beyond which the
  ## difference test is significant
  fariter<-0
  inc1<-0
  inc2<-1
  inc3<-1.5
  flag<-FALSE
  while(fariter<iterlim){
    D2<-lci_diff_test(p2,fitmodel,label,diff.method=diff.method)
    # FIXME: a bit of troubleshooting if we end up with an NA value
    # Assume this means the boundary is too far from the MLE?
    # If so, actually try doing bisection here to troubleshoot
    if(is.na(D2)){
      flag<-TRUE # to indicate that an NA was encountered already
      inc3<-inc2
      inc2<-mean(c(inc1,inc3))
      p2<-est+sign*init*se*inc2
    } else if (D2>crit) {
      break
    } else if (p2<=lb | p2>=ub){
      warning("LCI reached upper or lower boundary")
      # go ahead and return endpoint
      return(list(est=p2,D=D2,iter=0))
    } else {
      # if an NA was encountered before, proceed w/ bisection
      if(flag){
        inc1<-inc2
        inc2<-mean(c(inc1,inc3))
        p2<-est+sign*init*se*inc2
      } else {
        init<-init*1.1 # ad-hoc way of increasing upper boundary
      }
      p2<-est+sign*init*se*inc2
    }
    
    if(!is.null(lb)){
      p2[p2<lb]<-lb
    }
    if(!is.null(ub)){
      p2[p2>ub]<-ub
    }
    fariter<-fariter+1
  }
  
  ## Start bisection
  ## Set up initial values
  p0<-est
  p1<-(p2+p0)/2
  ## Fit three models
  ## starting point (MLE or closest to it)
  ## middle point
  ## farthest point from MLE
  D0<-lci_diff_test(p0,fitmodel,label,diff.method=diff.method)
  D1<-lci_diff_test(p1,fitmodel,label,diff.method=diff.method)
  D2<-lci_diff_test(p2,fitmodel,label,diff.method=diff.method)  
  count<-1
  flag<-FALSE
  while(!flag){
    
    ## enforce hard boundaries
    if(!is.null(lb)){
      p0[p0<lb]<-lb
      p1[p1<lb]<-lb
      p2[p2<lb]<-lb
    }
    if(!is.null(ub)){
      p0[p0>ub]<-ub
      p1[p1>ub]<-ub
      p2[p2>ub]<-ub
    }

    if(any(is.na(c(D0,D1,D2)))){
      return(NA)
    } else {
      
      ## Check which section is most likely to contain the estimate we want
      ## between closest and middle
      if(D0<crit & D1>crit){
        p2<-p1 ## farthest away moves to the middle
        D2<-D1
        p1<-(p2+p0)/2
      } else if (D1<crit & D2>crit){
        p0<-p1 ## closest to MLE moves to the middle
        D0<-D1
        p1<-(p2+p0)/2
        firstrun<-FALSE
      }      
      
      ## check for convergence
      D1<-lci_diff_test(p1,fitmodel,label,diff.method=diff.method)
      if(abs(D1-crit)<tol){
        flag<-TRUE
      }
      
      count<-count+1
      if(count>iterlim){break}
    }
  }
  return(list(est=p1,D=D1,iter=count-1,fariter=fariter))
  
}

profile_lci<-function(object,label,diff.method,grid){
  
  object<-lci_refit(object)

  D<-vector("numeric")
  for(j in 1:length(grid)){
    Dtmp<-try(lci_diff_test(grid[j],object,label,diff.method))
    if(!is.na(Dtmp)&class(Dtmp)!="try-error"){
      D<-c(D,Dtmp)
    } else {
      D<-c(D,NA)
    }
  }
  return(D)
}

# utility function for re-fitting models w/ constraints or different estimator, etc.
lci_refit<-function(object,const=NULL,estimator=NULL){
  
  # extract some info from fitted model
  prevmodel<-as.list(object@call)
  
  # extract function used to fit model (e.g., cfa,lavaan,etc.)
  f<-strsplit(as.character(prevmodel[1]),"::")[[1]]
  if(length(f)==1){
    f<-f[1]
  } else {
    f<-get(f[2],asNamespace(f[1]))
  }
  
  # parameter table
  ptab<-parTable(object)
  
  # add constraints, if desired
  if(!is.null(const)){
    if(is.list(const)){
      for(i in 1:length(const)){
        ptab.row<-lavaanify(paste(names(const)[i],"==",const[i],sep=""))[1,]
        ptab.row$plabel<-""
        ptab<-lav_partable_merge(ptab, ptab.row)
      }
    } else {
      stop("const should be a paired list of labels and numerical values")
    }
  }
  
  # change estimator, if desired
  if(!is.null(estimator)){
    prevmodel$estimator<-estimator
  }
  
  prevmodel$model<-ptab
  M<-try(do.call(f,prevmodel[-1]),silent=TRUE)
  
}
