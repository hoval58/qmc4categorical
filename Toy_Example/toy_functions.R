
source("../rsobol.R")

f<-function(x){
  return(exp(-x^2)*cos(x))
}

alloc_pow2 = function(alpha_list,M,tol=1e-9,verbose=FALSE,ansatz=0,rho=3){
  L=length(alpha_list)
  n=2^M
  if(rho=="infinity"){
    b=floor(log2(n/L))#all strata get initalized to this sample size
    m=rep(b,L)
    todo=n-L*2^b
    k=floor(todo/2^b)#number of strata to double
    remainders<-n/L-2^b
    if(k>0){
      top_k_indices <- order(remainders, decreasing = TRUE)[1:k]
      m[top_k_indices]=m[top_k_indices]+1
    }
    return(m)
  }
  else{
    
    if( ansatz==0 )
      beta = alpha_list^(2/(rho+1))
    else
      beta = alpha_list^(2/(rho+2))
    beta = beta/sum(beta) #renormalizes
    reportit = function(){
      if( verbose ){
        cat(2^m)
        cat("     Todo",sum(2^m)-n,"\n")
      }
    }        
    
    if( abs(sum(beta)-1.0) > tol )
      stop("Weights don't sum within tolerance of 1.")
    
    
    if( 2^M<L )
      stop("2^M was below L")
    m = rep(0,L)
    
    # step up the sample sizes
    
    todo = n - sum(2^m)
    reportit()
    while( todo > 0 ){
      if( min(2^m) > todo )
        stop("Samples cannot be allocated properly.")
      l = which.max( (beta/2^m)*(2^m<=todo))
      m[l] = m[l]+1
      todo = n - sum(2^m)
      reportit()
    }
    
    m
  }
}

alloc_int<-function(alpha_list,rho=2,M,ansatz=0){
  
  L=length(alpha_list)
  n=2^M
  if(rho=="infinity"){
    vec<-rep(n/L,L)
  }
  else{
    if( ansatz==0 ){
      beta = alpha_list^(2/(rho+1))}
    else{
      beta = alpha_list^(2/(rho+2))}
    beta = beta/sum(beta) 
    vec<-floor(beta*n)
    k<-n-sum(vec) #todo
    remainders<-n*beta-floor(beta*n)
    if(k>0){
      top_k_indices <- order(remainders, decreasing = TRUE)[1:k]
      vec[top_k_indices]=vec[top_k_indices]+1
    }
  }
  return(vec)
}

#a function that allocates a point in (0,1) to its stratum
which_stratum<-function(v,setting,alpha_list,m,rho){
  if (setting=="mc"|setting=="plain-rqmc"){
    cum_betas<-c(0,cumsum(alpha_list))
  }
  if (setting=="rqmc-adjusted"){
    betas<-alloc_int(alpha_list,rho=rho,M=m)/2^m
    cum_betas<-c(0,cumsum(betas))
  }
  if (setting=="rqmc-pow2"){
    betas<-2^alloc_pow2(alpha_list,m,rho=rho)/2^m
    cum_betas<-c(0,cumsum(betas))
  }
  stratum <- findInterval(v, cum_betas)
  return(stratum)
}

compute_estimator<-function(m, setting, seed, alpha_list, theta_list, rho=3){
  L=length(alpha_list)
  
  if (setting=="mc"){
    set.seed(seed)
    u<-matrix(runif((2^m)*s), nrow = 2^m, ncol = s)
  }
  if (setting=="plain-rqmc"|setting=="rqmc-adjusted"|setting=="rqmc-pow2"){
    u <- rsobol(m = m, s = s,seed=seed)
  }
  if (setting=="rqmc-indep"){
    n_l<-2^alloc_pow2(alpha_list,m,rho=rho)
    results=vector("list",L)
    for(i in 1:L){
      u<-rsobol(m=log2(n_l)[i],s=1,seed=i*seed)
      x2<-qnorm(u,mean=theta_list[i],sd=1)
      results[[i]]<-mean(f(x2))
    }
    return(sum(alpha_list*unlist(results)))
    
  }
  
  
  allocations <- which_stratum(u[, 1],setting,alpha_list,m,rho)
  # Find unique strata and sort them
  sorted_strata <- sort(unique(allocations))
  
  all_strata<-1:L
  stratum_groups <- lapply(all_strata, function(stratum) {
    subset <- u[allocations == stratum, , drop = FALSE]  # Keep matrix structure even if one row
    if (nrow(subset) == 0) {
      NULL  # Explicitly assign NULL for empty strata
    } else {
      as.matrix(subset)
    }
  })
  results=vector("list",L)
  for(i in 1:L){
    if(is.null(stratum_groups[[i]])){
      results[[i]]=0
    }
    else{
      x2<-qnorm(stratum_groups[[i]][,2],mean=theta_list[i],sd = 1)
      if (setting=="mc"|setting=="plain-rqmc"){
        results[[i]]<-sum(f(x2))
      }
      if (setting=="rqmc-pow2"){
        results[[i]]<-mean(f(x2))
      }
      if (setting=="rqmc-adjusted"){
        betas<-alloc_int(alpha_list,rho=rho,M=m)/2^m
        results[[i]]<-(alpha_list[i]/betas[i])*sum(f(x2))
      }
    }
  }
  if (setting=="mc"|setting=="plain-rqmc"|setting=="rqmc-adjusted"){
    return(sum(unlist(results))/2^m)
  }
  if (setting=="rqmc-pow2"){
    return(sum(alpha_list*unlist(results)))
  }
  
  
}

