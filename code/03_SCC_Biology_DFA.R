library(dplyr)
library(reshape2)
library(bayesdfa)

#Organize SCC biology data
load("data/ewidata.rda")

max_year = 2017 # hold out the 2018 data as a test set

sub_data<-ewidata[which(ewidata$year<=max_year),] # use to biology dataset to match time period of climate data
dfa_data = sub_data[which(sub_data$system=="SCC"&sub_data$subtype!="climate"),]
dfa_data = dplyr::arrange(dfa_data, code, year)

## Double check that time series are the correct ones
unique(dfa_data$code)

# Run if want to apply DFA to RREAS time series only CalCOFI time series only or
#keep <- grep("RREAS.",dfa_data$code)
#dfa_data=dfa_data[keep,]
#dfa_data=dfa_data[,c(1:3)]

## Log transform those data that are skewed.
melted = melt(dfa_data[, c("code", "year", "value")], id.vars = c("code", "year"))
dat <- dcast(melted, year ~ code)

write.csv(names(dat), file="ts_names_biology.csv")

n1 <- names(dat)[grepl('calcofi.',names(dat))]
n2 <- names(dat)[grepl('RREAS.',names(dat))]
ids <- c(n1,n2,"ZALOPHUS.PUPCT","ZALOPHUS.PUPWT")

for(i in 1:ncol(dat)){
  if(names(dat)[i] %in% ids){
    dat[,i] <- log(dat[,i])
  }
}

remelt = melt(dat,id.vars = "year")
names(remelt)<-c("year","code","value")
Y <- dcast(remelt, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

n_chains = 3
n_iter = 4000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(estimate_trend_ma = FALSE,
  estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(TRUE, FALSE),
  var_index = c("survey"), num_trends = 1:3,
  elpd_kfold = NA, se_elpd_kfold=NA)

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

for(i in 1:nrow(model_df)) {

  # These are a bbunch of switches for turning thigns on / off
  sigma_str = ""
  if(model_df$estimate_process_sigma[i] == TRUE) sigma_str = "_estSig"

  if(model_df$var_index[i]=="equal") varIndx = rep(1,nspp)
  if(model_df$var_index[i]=="unequal") varIndx = seq(1,nspp)
  if(model_df$var_index[i]=="survey") varIndx = c(rep(1,14), rep(2, 13), rep(3, 8), rep(4,3))

  str = "normal"
  if(model_df$est_nu[i]==TRUE) str = "student-t"
  theta_str = ""
  if(model_df$estimate_trend_ma[i]) theta_str = "_theta"
  phi_str = ""
  if(model_df$estimate_trend_ar[i]) phi_str = "_phi"
  sigma_str = ""
  if(model_df$estimate_process_sigma[i]) sigma_str = "_estSig"

  # Note that I have to set par_list = "all" to enable diagnostics. Using k-fold cross validation because of unstable Pareto-k
  n_cv = 10 # insted of folds, we'll make predictions for the last 10 years
  log_lik = matrix(0, nrow=n_chains*n_iter/2, ncol = n_cv)
  for(k in 1:n_cv) {
    ytrain = Y[,1:(ncol(Y)-k)]
    ytest = Y[,(ncol(Y)-k+1)]

    fit.mod = fit_dfa(y = ytrain,
      num_trends = model_df$num_trends[i],
      iter=n_iter,
      varIndx = varIndx,
      chains=n_chains, estimate_nu=model_df$est_nu[i],
      estimate_trend_ma = model_df$estimate_trend_ma[i],
      estimate_trend_ar = model_df$estimate_trend_ar[i],
      estimate_process_sigma = model_df$estimate_process_sigma[i],
      seed=123)
    # extract log_lik for new data
    pars = rstan::extract(fit.mod$model)

    for(j in 1:nrow(log_lik)) {
      # loop over iterations
      pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1)
      log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))

  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))

}

saveRDS(model_df, file="results_biology/biology_model_summary.rds")

# finally run the best model from the lfo-cv above -- largest is best
model_df = dplyr::arrange(model_df,-elpd_kfold)

if(model_df$var_index[1]=="equal") varIndx = rep(1,nspp)
if(model_df$var_index[1]=="unequal") varIndx = seq(1,nspp)
if(model_df$var_index[1]=="survey") varIndx = c(rep(1,14), rep(2, 13), rep(3, 8), rep(4,3))

fit.mod = fit_dfa(y = Y,
  num_trends = model_df$num_trends[1],
  iter=n_iter,
  varIndx = varIndx,
  chains=n_chains, estimate_nu=model_df$est_nu[1],
  estimate_trend_ma = model_df$estimate_trend_ma[1],
  estimate_trend_ar = model_df$estimate_trend_ar[1],
  estimate_process_sigma = model_df$estimate_process_sigma[1],
  seed=123)

# These are a bbunch of switches for turning thigns on / off
sigma_str = ""
if(model_df$estimate_process_sigma[i] == TRUE) sigma_str = "_estSig"

if(model_df$var_index[i]=="equal") varIndx = rep(1,nspp)
if(model_df$var_index[i]=="unequal") varIndx = seq(1,nspp)
if(model_df$var_index[i]=="survey") varIndx = c(rep(1,14), rep(2, 13), rep(3, 8), rep(4,3))

str = "normal"
if(model_df$est_nu[1]==TRUE) str = "student-t"
theta_str = ""
if(model_df$estimate_trend_ma[1]) theta_str = "_theta"
phi_str = ""
if(model_df$estimate_trend_ar[1]) phi_str = "_phi"
sigma_str = ""
if(model_df$estimate_process_sigma[1]) sigma_str = "_estSig"

saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[1],
  "_",model_df$var_index[1],"_",str,theta_str,phi_str,sigma_str,".rds"))
