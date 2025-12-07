

source("Toy_Example/toy_functions.R")
source("Toy_Example/plot_results_toy.R")

dir.create("Toy_Example/toy_results", recursive = TRUE, showWarnings = FALSE)

alphas=c(0.5,0.44,0.01,0.01,0.01,0.01,0.01,0.01) #mixture weights
thetas=c(0.7,1,1.5,1.6,1.7,1.8,1.9,2) #mixture components parameters
s=2 #dimension needed for all are random vectors: the first to allocate to one stratum, the second to sample

m_list=seq(3,12) #sample sizes (2^m_list)
n_rep=500 #number of replicates

### --- compute variance of the 5 estimators considered in the paper --- ###

var_mc_list<-list()
var_plain_rqmc_list<-list()
var_rqmc_adjusted_list<-list()
var_rqmc_pow2_list<-list()
var_rqmc_indep_list<-list()

for (j in 1:length(m_list)){
  cat("Running m =", m_list[[j]], "\n")
  mu_list_mc<-list()
  mu_list_plain_rqmc<-list()
  mu_list_rqmc_adjusted<-list()
  mu_list_rqmc_pow2<-list()
  mu_list_rqmc_indep<-list()
  
  for(i in 1:n_rep){
    mu_list_mc[[i]]<-compute_estimator(m_list[j],setting="mc",alpha_list=alphas,theta_list=thetas, seed=i)
    mu_list_plain_rqmc[[i]]<-compute_estimator(m_list[j],setting="plain-rqmc",alpha_list=alphas,theta_list=thetas, seed=i)
    mu_list_rqmc_adjusted[[i]]<-compute_estimator(m_list[j],setting="rqmc-adjusted",alpha_list=alphas,theta_list=thetas, seed=i,rho=2)
    mu_list_rqmc_pow2[[i]]<-compute_estimator(m_list[j],setting="rqmc-pow2",alpha_list=alphas,theta_list=thetas, seed=i,rho=3)
    mu_list_rqmc_indep[[i]]<-compute_estimator(m_list[j],setting="rqmc-indep",alpha_list=alphas,theta_list=thetas, seed=i,rho=3)
  }
  
  
  var_mc_list[[j]]<-var(unlist(mu_list_mc))
  var_plain_rqmc_list[[j]]<-var(unlist(mu_list_plain_rqmc))
  var_rqmc_adjusted_list[[j]]<-var(unlist(mu_list_rqmc_adjusted))
  var_rqmc_pow2_list[[j]]<-var(unlist(mu_list_rqmc_pow2))
  var_rqmc_indep_list[[j]]<-var(unlist(mu_list_rqmc_indep))
  
}
sample_sizes <- 2^seq(3, 12)

plot_variance(var_mc_list,var_plain_rqmc_list,var_rqmc_pow2_list,
              var_rqmc_adjusted_list,var_rqmc_indep_list, sample_sizes, 
              file = "Toy_Example/toy_results/variances_toy.png")
