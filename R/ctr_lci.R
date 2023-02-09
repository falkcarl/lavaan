## Author: Carl F. Falk
##
## General purpose likelihood-based confidence interval function using lavaan
## Requires that the parameter or quantity of interest is clearly labeled in the parameter table 
lci<-function(object, label, level=.95, bound=c("lower","upper"),
              optimizer="Rsolnp",
              ci.method="NealeMiller1997",
              start=NULL,diff.method="default",Dtol=.05,
              reoptimize=FALSE,
              iterlim=50,control=list(),...){

  ## input checking
  if(class(object)!="lavaan"){
    stop("Object must be a fitted lavaan model")
  }
  if(!is.null(start)&length(bound)>1){
    stop("Custom starting values only supported when obtaining one boundary at a time. Use 'lower' or 'upper' for the bound argument, not both.")
  }
  if(level <=0 | level >=1){
    stop("level must be between 0 and 1")
  }
  if(is.null(start)){
    if((ci.method=="bisect"|ci.method=="uniroot") & object@Options$se=="none"){
      warning("Bisection requires either standard errors to aid in determining where to start the algorithm (change se='none' to something else) or a custom start value.")
      object@Options$se<-"standard"
      object<-lci_refit(object)
    }
  }
  if(!ci.method %in% c("NealeMiller1997","bisect","uniroot")){
    stop("Unsupported CI estimation method.")
  }
  
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
    #crit<-chat*crit
    
    object@Options$estimator<-"ML"
    object@Options$se<-"standard"
    object@Options$test<-"standard"
    fit<-lci_refit(object)
    
  } else {
    chat<-1
    fit<-lci_refit(object)
  }
  
  ## Start list of output before embarking on optimization
  result<-list()
  result$est<-est
  
  ## Prepare critical value of chi-square
  crit<-qchisq(level,1) ## hardcoded 1 df for now
  result$crit<-crit
  
  # Begin optimization for each bound
  for(b in bound){
    
    LCI<-lci_internal(fit, label, level, est, ptable, pindx,
                      crit, b, chat, optimizer, ci.method,
                      start,diff.method, Dtol, iterlim, control,...)
    
    # Troubleshoot any problems
    if(!is.na(LCI$bound)){
      if(LCI$bound<=result$est & b=="upper" | LCI$bound>=result$est & b=="lower"){
        warning("CI boundary on wrong side of estimate ", b)
        if(reoptimize & ci.method=="NealeMiller1997"){
          warning("Attempting to re-optimize by moving starting values")

          # Obtain SE
          #pest<-parameterEstimates(object)
          se<-ptable$se[pindx[1]]
        
          # Obtain new starting model with function of interest moved slightly towards boundary
          newstart<-ifelse(b=="upper", est+.25*se,est-.25*se)
        
          # Make necessary adjustments to start values in case of function of model parameters
          if(ptable$op[pindx]==":=" & (diff.method[1]=="satorra.2000"|object@Options$estimator=="ML"& !diff.method[1] %in% c("satorra.bentler.2010","satorra.bentler.2001"))){
            const<-list()
            const[[label]]<-newstart
            tmp<-lci_refit(fit,const)
            ptabtmp<-parTable(tmp)
            tmplabel<-all.vars(parse(file="", text=ptabtmp$rhs[pindx]))
            newstart<-ptabtmp$est[match(tmplabel,ptabtmp$label)]
          }
        
          # Call lci_internal again with start replaced by newstart
          LCI<-lci_internal(fit, label, level, est, ptable, pindx,
                            crit, b, chat, optimizer, ci.method,
                            newstart,diff.method, Dtol, iterlim, control,...)
        }
      } else if (abs(LCI$D-crit)>Dtol) {
        warning("Tolerance level for chi-square difference test not met for CI boundary ", b,".
                  An obtained difference test at the CI boundary might not be close to the desired critical value")
        if(reoptimize & ci.method=="NealeMiller1997"){
          warning("Attempting to re-optimize by bisection")
        
          tmpstart<-LCI$bound
        
          if (LCI$D<crit){
            # If there is bias, usually we miss with a D that is too high
            # Could be a sign that no boundary can be found due to insufficient information
          
            # Obtain standard errors
            #pest<-parameterEstimates(object)
            se<-ptable$se[pindx[1]]
          
            # Adjust starting values
            tmpstart<-LCI$bound
            tmpstart<-ifelse(b=="upper", tmpstart+.5*se,tmpstart-.5*se)
          }
        
          # Call lci_internal via bisection
                            crit, b, chat, optimizer, "bisect",
          LCI<-lci_internal(fit, label, level, est, ptable, pindx,
                            tmpstart,diff.method, Dtol, iterlim, control,...)
        }
      }
    } else {
      warning(b, " boundary not found")
    }
    
    # Save results
    if(b=="upper"){
      result$upper<-LCI$bound
      result$convupper<-LCI$conv
      result$Dupper<-LCI$D
    }
    if(b=="lower"){
      result$lower<-LCI$bound
      result$convlower<-LCI$conv
      result$Dlower<-LCI$D
    }
  }
  return(result)
}

lci_internal<-function(object, label, level, est, ptable,
                       pindx, crit, bound, chat,
                       optimizer="Rsolnp",
                       ci.method="NealeMiller1997",
                       start=NULL,diff.method="default",Dtol=.05,
                       iterlim=50,control=list(),...){
  
  ret<-list()
  D<-NULL
  cont<-NULL
  est.bound<-NULL
  if(ci.method[1]=="NealeMiller1997"){
    fitfunc<-lci_nealemiller1997
    
    if(ptable$op[pindx]==":=" & (diff.method[1]=="satorra.2000"|object@Options$estimator=="ML"& !diff.method[1] %in% c("satorra.bentler.2010","satorra.bentler.2001"))){
      label<-all.vars(parse(file="", text=ptable$rhs[pindx]))
      if(length(start)!=length(label)&!is.null(start)){
        warning("Number of starting values did not match number of parameters in quantity of interest. See documentation for \"start\" argument.")
      }
      if(length(start)!=length(label)){
        start<-ptable$est[match(label,ptable$label)]
      }
      diff.method<-"default"
    } else {
      if(is.null(start)){
        start<-est
      }
    }
    
    if(optimizer[1]=="Rsolnp"){
      if(is.null(control$trace)){
        control$trace<-0
      }
      LCI<-suppressWarnings(try(Rsolnp::solnp(start,fitfunc,fitmodel=object,
                                              label=label, pindx=pindx, crit=crit, bound=bound,
                                              diff.method=diff.method,chat=chat,control=control,...)))
      
      if(class(LCI)!="try-error"){
        D<-lci_diff_test(LCI$pars,object,label,diff.method=diff.method,chat=chat)
        est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
        attr(D,"mod")<-NULL
        conv<-LCI$convergence
      } else{
        est.bound<-NA
      }
    } else if (optimizer[1]=="optim"){
      LCI<-suppressWarnings(try(optim(start, fitfunc,fitmodel=object,
                                      label=label, pindx=pindx, crit=crit, bound=bound,
                                      diff.method=diff.method,chat=chat,control=control,...)))
      if(class(LCI)!="try-error"){
        D<-lci_diff_test(LCI$par,object,label,diff.method=diff.method,chat=chat)
        est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
        attr(D,"mod")<-NULL          
        conv<-LCI$convergence
      } else{
        est.bound<-NA
      }
    } else if (optimizer[1]=="nlminb"){
      LCI<-suppressWarnings(try(nlminb(start, fitfunc,fitmodel=object,
                                       label=label, pindx=pindx, crit=crit, bound=bound,
                                       diff.method=diff.method,chat=chat,control=control,...)))
      if(class(LCI)!="try-error"){
        D<-lci_diff_test(LCI$par,object,label,diff.method=diff.method,chat=chat)
        est.bound<-parameterEstimates(attr(D,"mod"))$est[pindx[1]]
        attr(D,"mod")<-NULL
        conv<-LCI$convergence
      } else{
        est.bound<-NA
      }
    } else {
      stop("Unsupported optimizer")
    }
  } else if (ci.method[1]=="bisect"){
    if(diff.method[1]=="satorra.2000"){
      diff.method<-"default"
    }
    LCI<-suppressWarnings(lci_bisect(object,label,level,crit,tol=Dtol*.002,bound=bound,
                                     diff.method=diff.method,iterlim=iterlim,
    if(class(LCI)!="try-error"){
                                     chat=chat,start=start,...))
      D<-lci_diff_test(LCI$est,object,label,diff.method=diff.method,chat=chat)
      attr(D,"mod")<-NULL
      est.bound<-LCI$est
      conv<-LCI$iter
    } else{
      est.bound<-NA
    }
  } else if (ci.method[1]=="uniroot"){
    if(diff.method[1]=="satorra.2000"){
      diff.method<-"default"
    }
    LCI<-suppressWarnings(lci_uniroot(object,label,level,crit,tol=Dtol*.002,bound=bound,
                                     diff.method=diff.method,iterlim=iterlim,
                                     chat=chat,start=start,...))
    if(class(LCI)!="try-error" & any(!is.na(LCI))){
      D<-lci_diff_test(LCI$est,object,label,diff.method=diff.method,chat=chat)
      attr(D,"mod")<-NULL
      est.bound<-LCI$est
      conv<-LCI$iter
    } else{
      est.bound<-NA
      conv<-NA
      D<-NA
    }
  } else {
    stop("Unsupported ci.type and/or optimizer")
  }
  
  ret$bound<-est.bound
  ret$conv<-conv
  ret$D<-D
  
  return(ret)
}

# Just computes constrains quantity in label to p and computes difference test for use with lci functions
lci_diff_test<-function(p,fitmodel,label,diff.method="default",chat=1){

  const<-list()
  for(i in 1:length(p)){
    const[[label[i]]]<-p[i]
  }
  M0<-lci_refit(fitmodel,const)

  if(class(M0)=="try-error"){
    return(NA)
  } else {
    if(M0@Fit@converged){
      D<-try(lavTestLRT(fitmodel,M0,method=diff.method,A.method="exact")$`Chisq diff`[2])
      if(class(D)!="try-error"){
        D<-D/chat
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
                              diff.method="default",chat=1){
  
  D<-lci_diff_test(p,fitmodel,label,diff.method,chat=chat)
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

lci_bisect<-function(fitmodel,label,level,crit,tol=1e-5,iterlim=50,init=qnorm((1-level)/2,lower.tail=FALSE),
                     bound=c("lower","upper"),diff.method="default",lb=-Inf,ub=Inf,chat=1,start=NULL,...){
  ## FIXME: init<1 might break things
  ## FIXME: ub and lb aren't available to user
  ## FIXME: tol and init also not available to user
  
  ## extract parameter table
  ptable<-parTable(fitmodel)
  
  ## Look for parameter label
  pindx<-which(ptable$label==label)
  
  if(length(pindx)<1){
    stop("Parameter label not found in lavaan parameter table.")
  }
  
  ## point estimate based on fitted model
  est<-ptable$est[pindx[1]]
  ## standard error
  se<-ptable$se[pindx[1]]
  
  ## p0 is at the MLE
  ## We need to find p2, which is a point at which the difference test is significant
  p0<-est
  if(bound=="lower"){
    sign<- -1
  } else {
    sign<-1
  }
  
  if(!is.null(start)){
    ## Custom p2
    if((start<est & bound=="upper")|(start>est & bound=="lower")){
      stop("Custom starting value for bisection on wrong side of estimate.") 
    }
    
    p2<-start
    
  } else {
    ## Guess a starting point for p2 based on point estimate and multiplier w/ se
    p2<-est+sign*init*se
  }
  
  ## enforce hard boundaries
  if(!is.null(lb)){
    p2[p2<lb]<-lb
  }
  if(!is.null(ub)){
    p2[p2>ub]<-ub
  }

  ## Try p2 and troubleshoot as necessary
  fariter<-0
  inc1<-0
  inc2<-ifelse(is.null(start),1,(start-est)/(sign*init*se))
  inc3<-ifelse(is.null(start),1.5,1.5*(start-est)/(sign*init*se))
  flag<-FALSE
  while(fariter<iterlim){
    D2<-lci_diff_test(p2,fitmodel,label,diff.method=diff.method,chat=chat)
    # A bit of troubleshooting if we end up with an NA value
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
  D0<-lci_diff_test(p0,fitmodel,label,diff.method=diff.method,chat=chat)
  D1<-lci_diff_test(p1,fitmodel,label,diff.method=diff.method,chat=chat)
  D2<-lci_diff_test(p2,fitmodel,label,diff.method=diff.method,chat=chat)  
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
      D1<-lci_diff_test(p1,fitmodel,label,diff.method=diff.method,chat=chat)
      if(abs(D1-crit)<tol){
        flag<-TRUE
      }
      
      count<-count+1
      if(count>iterlim){break}
    }
  }
  return(list(est=p1,D=D1,iter=count-1,fariter=fariter))
  
}

lci_uniroot<-function(fitmodel,label,level,crit,tol=1e-5,iterlim=50,init=qnorm((1-level)/2,lower.tail=FALSE),
                     bound=c("lower","upper"),diff.method="default",lb=-Inf,ub=Inf,chat=1,start=NULL,...){

  ## extract parameter table
  ptable<-parTable(fitmodel)
  
  ## Look for parameter label
  pindx<-which(ptable$label==label)
  
  if(length(pindx)<1){
    stop("Parameter label not found in lavaan parameter table.")
  }
  
  ## point estimate based on fitted model
  est<-ptable$est[pindx[1]]
  ## standard error
  se<-ptable$se[pindx[1]]
  
  # function to minimize: difference for difference test
  f<-function(x,fitmodel,label,diff.method,chat,crit){
    D<-lci_diff_test(x,fitmodel,label,diff.method=diff.method,chat=chat)
    if(is.na(D)){
      return(.Machine$double.xmax) # ensures that uniroot has something to work with
    } else{
      return(D-crit)
    }
  }
  
  ## p0 is at the MLE
  p0<-est
  
  ## Pick another endpoint for boundary
  sign<-ifelse(bound=="lower",-1,1)
  
  if(!is.null(start)){
    
    ## Custom p2
    if((start<est & bound=="upper")|(start>est & bound=="lower")){
      stop("Custom starting value for bisection on wrong side of estimate.") 
    }
    p2<-start
    
  } else {
    
    ## Guess a starting point for p2 based on point estimate and multiplier w/ se
    p2<-est+sign*init*se*sqrt(crit)
    
    tmp<-f(p2,fitmodel,label,diff.method,chat,crit)
    count<-1
    while(!tmp>0 & count < iterlim){
      init<-init*1.1
      p2<-est+sign*init*se*sqrt(crit)
      tmp<-f(p2,fitmodel,label,diff.method,chat,crit)
      count<-count+1
    }
  }
  
  # obtain result using uniroot
  # interval is between p0 and p2
  result<-try(uniroot(f,c(p0,p2), fitmodel=fitmodel, label=label, diff.method=diff.method,
                  chat=chat, crit=crit, tol=tol, maxiter=iterlim),silent=TRUE)
  
  if(class(result)!="try-error"){
    endpoint<-result$root
    D<-lci_diff_test(est,fitmodel,label,diff.method=diff.method,chat=chat)
    iter<-result$iter
    return(list(est=endpoint,D=D,iter=iter))
  } else {
    return(NA)
  }
  
}

lciProfile<-function(object,label,diff.method="default",grid=NULL){
  
  if(class(object)!="lavaan"){
    stop("object must be a fitted lavaan model")
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
  
  # if grid is empty, try to come up with something reasonable
  if(is.null(grid)){
    est<-ptable$est[pindx][1]
    se<-try(ptable$se[pindx][1])
    if(object@Options$se=="none"|class(se)=="try-error"|is.na(se))
      stop("If grid is not specified, standard errors must be available for the fitted model to determine range of x-axis values.")
    grid<-seq(est-se*2,est+se*2,length.out=101)
  }
  
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
  return(list(grid=grid,D=D))
}

# utility function for re-fitting models w/ constraints or different estimator, etc.
lci_refit<-function(object,const=NULL,estimator=NULL){
  
  # extract some info from fitted model
  #prevmodel<-as.list(object@call)
  
  # extract function used to fit model (e.g., cfa,lavaan,etc.)
  #f<-strsplit(as.character(prevmodel[1]),"::")[[1]]
  #if(length(f)==1){
  #  f<-f[1]
  #} else {
  #  f<-get(f[2],asNamespace(f[1]))
  #}
  #f<-"lavaan"
  
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
    #prevmodel$estimator<-estimator
    object@Options$estimator<-estimator
  }
  
  #prevmodel$model<-ptab
  #object@Model<-ptab
  object@ParTable<-as.list(ptab)

  M<-lavaan(slotOptions = object@Options,
         slotParTable = object@ParTable,
         slotSampleStats = object@SampleStats,
         slotData = object@Data,
         slotCache = object@Cache)
  #M<-try(do.call(f,prevmodel[-1]),silent=TRUE)
  
  return(M)
  
}
