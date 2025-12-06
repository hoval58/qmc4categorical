
source("Flood_Example/flood_model.R")
source("Flood_Example/plot_results_flood.R")

dir.create("Flood_Example/flood_results", recursive = TRUE, showWarnings = FALSE)

#fixed parameters of the Saint-Venant flood model

B=300 #fixed
L=5000
s=5 #dimension of thw random vectors: 1 for selecting the stratum + 4 for sampling the 4 variables

#parameters of the Frechet and Gamma distributions in the mixture

alphas<-c(0.95,0.02,0.02,0.01) #mixture weights
shape_Q<-c(6,6,6,6)
scale_Q<-c(1300,3900,1300,3900)
shape_K<-c(90,90,15,15)
scale_K<-c(1/3,1/3,1,1)

m_list=seq(3,12) #sample sizes
n_rep=500 #number of replicates

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
    mu_list_mc[[i]]<-compute_estimator(m_list[j],setting="mc",alpha_list=alphas,shape_K,scale_K,shape_Q,scale_Q, seed=i)
    mu_list_plain_rqmc[[i]]<-compute_estimator(m_list[j],setting="plain-rqmc",alpha_list=alphas,shape_K,scale_K,shape_Q,scale_Q, seed=i)
    mu_list_rqmc_adjusted[[i]]<-compute_estimator(m_list[j],setting="rqmc-adjusted",alpha_list=alphas,shape_K,scale_K,shape_Q,scale_Q, seed=i,rho=2)
    mu_list_rqmc_pow2[[i]]<-compute_estimator(m_list[j],setting="rqmc-pow2",alpha_list=alphas,shape_K,scale_K,shape_Q,scale_Q, seed=i,rho=2)
    mu_list_rqmc_indep[[i]]<-compute_estimator(m_list[j],setting="rqmc-indep",alpha_list=alphas,shape_K,scale_K,shape_Q,scale_Q, seed=i,rho=2)
  }
  
  
  var_mc_list[[j]]<-var(unlist(mu_list_mc))
  var_plain_rqmc_list[[j]]<-var(unlist(mu_list_plain_rqmc))
  var_rqmc_adjusted_list[[j]]<-var(unlist(mu_list_rqmc_adjusted))
  var_rqmc_pow2_list[[j]]<-var(unlist(mu_list_rqmc_pow2))
  var_rqmc_indep_list[[j]]<-var(unlist(mu_list_rqmc_indep))
  
}


sample_sizes <- 2^seq(3, 12)

var_lists <- list(
  "MC" = var_mc_list,
  "RQMC" = var_plain_rqmc_list,
  "RQMC powers of 2" = var_rqmc_pow2_list,
  "RQMC L independent" = var_rqmc_indep_list,
  "RQMC adjusted" = var_rqmc_adjusted_list
)

# Plot and save figure
plot_variance(var_lists, sample_sizes, file = "Flood_example/flood_results/variance_plot.png")


### --- compute variance of the importance adjusted estimator with different rhos --- ###


var_infty_list <- list()
var_rho_1_list <- list()
var_rho_2_list <- list()
var_rho_3_list <- list()

for (j in 1:length(m_list)) {
  cat("Running m =", m_list[[j]], "\n")
  
  # temp lists for each replication
  mu_list_infty <- vector("list", n_rep)
  mu_list_rho_1 <- vector("list", n_rep)
  mu_list_rho_2 <- vector("list", n_rep)
  mu_list_rho_3 <- vector("list", n_rep)
  
  # replications
  for (i in 1:n_rep) {
    mu_list_infty[[i]] <- compute_estimator(
      m_list[j], setting = "rqmc-adjusted", alpha_list = alphas,
      shape_K, scale_K, shape_Q, scale_Q, seed = i, rho = "infinity"
    )
    mu_list_rho_1[[i]] <- compute_estimator(
      m_list[j], setting = "rqmc-adjusted", alpha_list = alphas,
      shape_K, scale_K, shape_Q, scale_Q, seed = i, rho = 1
    )
    mu_list_rho_2[[i]] <- compute_estimator(
      m_list[j], setting = "rqmc-adjusted", alpha_list = alphas,
      shape_K, scale_K, shape_Q, scale_Q, seed = i, rho = 2
    )
    mu_list_rho_3[[i]] <- compute_estimator(
      m_list[j], setting = "rqmc-adjusted", alpha_list = alphas,
      shape_K, scale_K, shape_Q, scale_Q, seed = i, rho = 3
    )
  }
  
  # store variances
  var_infty_list[[j]] <- var(unlist(mu_list_infty))
  var_rho_1_list[[j]] <- var(unlist(mu_list_rho_1))
  var_rho_2_list[[j]] <- var(unlist(mu_list_rho_2))
  var_rho_3_list[[j]] <- var(unlist(mu_list_rho_3))
}



var_lists_diff_rhos <- list(
  "rho=infinity" = var_infty_list,
  "rho=1" = var_rho_1_list,
  "rho=2" = var_rho_2_list,
  "rho=3" = var_rho_3_list
)


# plot and save figure

plot_rho_variance(
  var_lists_diff_rhos,
  sample_sizes,
  slope = -1.8,
  file = "Flood_Example/flood_results/diff_rhos.png"
)


