
source("rsobol.R")

library(VGAM)

H<-function(Q,K_s,Z_m,Z_v){
  slope<-(Z_m-Z_v)/L
  return((Q/(K_s*B*sqrt(slope)))^0.6)
}

alloc_pow2 <- function(alpha_list,M,tol=1e-9,verbose=FALSE,ansatz=0,rho=3){
  # Computes an allocation of 2^M samples across L strata such that each 
  # stratum receives a power-of-two number of samples, chosen to approximate 
  # an ideal proportional allocation.
  L=length(alpha_list)
  n=2^M
  if(rho=="infinity"){
    b=floor(log2(n/L))#all strata get initialized to the same sample size
    m=rep(b,L)
    todo=n-L*2^b
    k=floor(todo/2^b)#number of strata to double
    remainders<-rep(n/L - 2^b, L)
    if(k>0){
      top_k_indices <- order(remainders, decreasing = TRUE)[1:k]
      m[top_k_indices]=m[top_k_indices]+1
    }
    return(m)
  }
  else{
    
    if( ansatz==0 )
      beta = alpha_list^(2/(rho+1)) #ideal proportion of samples under ansatz 0
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
      l = which.max( (beta/2^m)*(2^m<=todo)) #increment the most under represented of the eligible strata
      m[l] = m[l]+1
      todo = n - sum(2^m)
      reportit()
    }
    
    m # a vector containing the exponents for each stratum
  }
}


alloc_int<-function(alpha_list,rho=2,M,ansatz=0){
  # Computes an integer allocation of 2^M samples across L strata by rounding 
  # an optimal proportional allocation (based on Eq. 14 of the paper).
  
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

compute_estimator<-function(m, setting,seed,alpha_list,shape_K,scale_K,shape_Q,scale_Q,rho=3){
  # Computes an estimator of E[H(Q, K_s, Z_m, Z_v)] using either Monte Carlo (MC), 
  # plain randomized quasi–Monte Carlo (RQMC), RQMC with adjusted stratified 
  # weights, or RQMC with power-of-two independent sampling. 
  #
  # The function supports several sampling settings:
  #
  #   • "mc"                – Standard Monte Carlo sampling using IID uniforms.
  #   • "plain-rqmc"        – Sobol-based RQMC without stratification adjustment.
  #   • "rqmc-adjusted"     – RQMC with importance-adjusted stratum weights 
  #                           using integer allocations (alloc_int).
  #   • "rqmc-pow2"         – RQMC with power-of-two allocation per stratum 
  #                           (alloc_pow2), averaging within each stratum.
  #   • "rqmc-indep"        – Independent Sobol sequences for each stratum, each 
  #                           using 2^{m_l} samples where m_l comes from alloc_pow2.
  
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
      u<-rsobol(m=log2(n_l)[i],s=4,seed=i*seed)
      Z_m<-qunif(u[,1],54,56)
      Z_v<-qunif(u[,2],49,51)
      Q<-qfrechet(u[,3],scale=scale_Q[i],shape=shape_Q[i])
      K_s<-qgamma(u[,4],shape=shape_K[i],scale=scale_K[i])
      results[[i]]<-mean(H(Q,K_s,Z_m,Z_v))
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
      Z_m<-qunif(stratum_groups[[i]][,2],54,56)
      Z_v<-qunif(stratum_groups[[i]][,3],49,51)
      Q<-qfrechet(stratum_groups[[i]][,4],scale=scale_Q[i],shape=shape_Q[i])
      K_s<-qgamma(stratum_groups[[i]][,5],shape=shape_K[i],scale=scale_K[i])
      if (setting=="mc"|setting=="plain-rqmc"){
        results[[i]]<-sum(H(Q,K_s,Z_m,Z_v))
      }
      if (setting=="rqmc-pow2"){
        results[[i]]<-mean(H(Q,K_s,Z_m,Z_v))
      }
      if (setting=="rqmc-adjusted"){
        betas<-alloc_int(alpha_list,rho=rho,M=m)/2^m
        results[[i]]<-(alpha_list[i]/betas[i])*sum(H(Q,K_s,Z_m,Z_v))
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