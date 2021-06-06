library(dplyr)
library(reshape2)
library(bayesdfa)

#Organize SCC biology data
max_year = 2017
min_year = 1981
load("data/ewidata.rda")
sub_data<-ewidata[which(ewidata$year>=min_year&ewidata$year<=max_year),] 
dfa_data = sub_data[which(sub_data$system=="SCC"&sub_data$subtype!="climate"),]

## Double check that time series are the correct ones
unique(dfa_data$code)

# Run if want to apply DFA to RREAS time series only or CalCOFI time series only
#keep <- grep("RREAS.",dfa_data$code)
#dfa_data=dfa_data[keep,]
#dfa_data=dfa_data[,c(1:3)]

## Log transform those data that are skewed.
melted = melt(dfa_data[, c("code", "year", "value")], id.vars = c("code", "year"))
dat <- dcast(melted, year ~ code)

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

##### Local Environmental Variables

#Using data from central California current 
dat<-read.csv("data/roms.july_june.mean.csv") # the south and central data (July-June are also in the ewidata file
dat=dat[,c(1:7)] 
avgROMs=dat
colnames(avgROMs)<-c("year","sst","ssh","ild","bv","cuti","beuti")
future_roms = dplyr::filter(avgROMs, year==max_year+1)
avgROMs = dplyr::filter(avgROMs,year <= max_year)

#Environmental data can also be accessed via the ewidata file
#load("/Users/mary.hunsicker/ewidata/data/ewidata.rda")
#sub_data<-ewidata[which(ewidata$year>= 1981&ewidata$year<2019),]
#dfa_data = sub_data[which(sub_data$system=="SCC"&sub_data$subtype=="climate"),]

#subset avgROMs for years of interest
avgROMs.short = subset(avgROMs, year>= min_year & year<=(max_year+1),select=year:beuti)

nspp=nrow(Y)
nyear=ncol(Y)
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

#Create separate data frames for each env variable
##sst
sst.value<-avgROMs.short$sst
cov.sst<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.sst$timeseries<-rep(1:nspp, each=nyear)
cov.sst$time<-rep(seq(1:nyear), times=nspp)
cov.sst$value<-rep(sst.value, times=nspp)
cov.sst$covariate<-as.integer("1")

##ssh
ssh.value<-avgROMs.short$ssh
cov.ssh<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.ssh$timeseries<-rep(1:nspp, each=nyear)
cov.ssh$time<-rep(seq(1:nyear), times=nspp)
cov.ssh$value<-rep(ssh.value, times=nspp)
cov.ssh$covariate<-as.integer("1")

##ild
ild.value<-avgROMs.short$ild
cov.ild<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.ild$timeseries<-rep(1:nspp, each=nyear)
cov.ild$time<-rep(seq(1:nyear), times=nspp)
cov.ild$value<-rep(ild.value, times=nspp)
cov.ild$covariate<-as.integer("1")

##BV
bv.value<-avgROMs.short$bv
cov.bv<-data.frame(time=numeric(n_row),
                   timeseries=numeric(n_row),
                   covariate=numeric(n_row),
                   value=numeric(n_row))
cov.bv$timeseries<-rep(1:nspp, each=nyear)
cov.bv$time<-rep(seq(1:nyear), times=nspp)
cov.bv$value<-rep(bv.value, times=nspp)
cov.bv$covariate<-as.integer("1")

# scale bv
cov.bv$value = cov.bv$value*100 

##BEUTI
beuti.value<-avgROMs.short$beuti
cov.beuti<-data.frame(time=numeric(n_row),
                      timeseries=numeric(n_row),
                      covariate=numeric(n_row),
                      value=numeric(n_row))
cov.beuti$timeseries<-rep(1:nspp, each=nyear)
cov.beuti$time<-rep(seq(1:nyear), times=nspp)
cov.beuti$value<-rep(beuti.value, times=nspp) 
cov.beuti$covariate<-as.integer("1")

##cuti
cuti.value<-avgROMs.short$cuti
cov.cuti<-data.frame(time=numeric(n_row),
                     timeseries=numeric(n_row),
                     covariate=numeric(n_row),
                     value=numeric(n_row))
cov.cuti$timeseries<-rep(1:nspp, each=nyear)
cov.cuti$time<-rep(seq(1:nyear), times=nspp)
cov.cuti$value<-rep(cuti.value, times=nspp)
cov.cuti$covariate<-as.integer("1")

######## Run models
n_chains = 3
n_iter = 4000 # bump up to 4000?

options(mc.cores = parallel::detectCores())


model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(FALSE),
                       var_index = c("survey"), num_trends = 1, #option 1:3
                       elpd_kfold = NA, se_elpd_kfold=NA,
                       cov1 = c("sst","ssh","ild","bv","beuti", "cuti"),
                       cov2 = c("sst","ssh","ild","bv","beuti", "cuti"),
                       future_roms = NA
)

df_keep = expand.grid(cov1 = c("sst","ssh","ild","bv","beuti", "cuti"),
                      cov2= c("sst","ssh","ild","bv","beuti", "cuti"))
df_keep$keep = 1
for(i in 1:nrow(df_keep)) {
  i_indx = match(df_keep$cov1[i], c("sst","ssh","ild","bv","beuti", "cuti"))
  j_indx = match(df_keep$cov2[i], c("sst","ssh","ild","bv","beuti", "cuti"))
  if(j_indx >= i_indx) df_keep$keep[i] = 0
}
model_df = left_join(model_df,df_keep) %>% dplyr::filter(keep == 1)

# add in predictions for next year
# model_df$future_roms[which(model_df$cov=="sst")] = future_roms$sst
# model_df$future_roms[which(model_df$cov=="ssh")] = future_roms$ssh
# model_df$future_roms[which(model_df$cov=="ild")] = future_roms$ild
# model_df$future_roms[which(model_df$cov=="bv")] = future_roms$bv
# model_df$future_roms[which(model_df$cov=="beuti")] = future_roms$beuti
# model_df$future_roms[which(model_df$cov=="cuti")] = future_roms$cuti

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}
for(i in c(15)) {
  
  # These are a bunch of switches for turning things on / off
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
  
  if(model_df$cov1[i] == "sst") {
    obs_covar1 = cov.sst
  }
  if(model_df$cov1[i] == "ssh") {
    obs_covar1 = cov.ssh
  }
  if(model_df$cov1[i] == "ild") {
    obs_covar1 = cov.ild
  }
  if(model_df$cov1[i] == "bv") {
    obs_covar1 = cov.bv
  }
  if(model_df$cov1[i] == "beuti") {
    obs_covar1 = cov.beuti
  }
  if(model_df$cov1[i] == "cuti") {
    obs_covar1 = cov.cuti
  }
  
  if(model_df$cov2[i] == "sst") {
    obs_covar2 = cov.sst
  }
  if(model_df$cov2[i] == "ssh") {
    obs_covar2 = cov.ssh
  }
  if(model_df$cov2[i] == "ild") {
    obs_covar2 = cov.ild
  }
  if(model_df$cov2[i] == "bv") {
    obs_covar2 = cov.bv
  }
  if(model_df$cov2[i] == "beuti") {
    obs_covar2 = cov.beuti
  }
  if(model_df$cov2[i] == "cuti") {
    obs_covar2 = cov.cuti
  } 
  obs_covar2$covariate = 2
  obs_covar = rbind(obs_covar1, obs_covar2)
  
  # Note that we have to set par_list = "all" to enable diagnostics. 
  # Added k-fold cross validation because of unstable Pareto-k
  n_cv = 10 # instead of folds, we'll make predictions for the last 10 years
  log_lik = matrix(0, nrow=n_chains*n_iter/2, ncol = n_cv)
  for(k in 1:n_cv) {
    year_train_max = (max_year - k)
    year_test = max_year - k + 1
    ytrain = Y[,1:(ncol(Y)-k)]
    ytest = Y[,(ncol(Y)-k+1)]
    covar_train = dplyr::filter(obs_covar, time <= (year_train_max - min(avgROMs.short$year) + 1))
    covar_test = dplyr::filter(obs_covar, time==(year_test - min(avgROMs.short$year) + 1))
    fit.mod = fit_dfa(y = ytrain,
                      num_trends = model_df$num_trends[i],
                      iter=n_iter,
                      varIndx = varIndx,
                      chains=n_chains, 
                      estimate_nu=model_df$est_nu[i],
                      estimate_trend_ma = model_df$estimate_trend_ma[i],
                      estimate_trend_ar = model_df$estimate_trend_ar[i],
                      estimate_process_sigma = model_df$estimate_process_sigma[i],
                      seed=123,
                      obs_covar = covar_train)
    # extract log_lik for new data
    pars = rstan::extract(fit.mod$model)
    
    for(j in 1:nrow(log_lik)) {
      # loop over iterations
      pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1) + 
        pars$b_obs[j,1,] * covar_test$value[which(covar_test$covariate==1)] + 
        pars$b_obs[j,2,] * covar_test$value[which(covar_test$covariate==2)]
      log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))
  
  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  saveRDS(model_df, file="results_covariate/covariate_model_summary_2covariates.rds")
}
# finally run the best model from the lfo-cv above -- largest is best
model_df = readRDS("results_covariate/covariate_model_summary.rds")
model_df = dplyr::arrange(model_df,-elpd_kfold)

i = 1 # ran with 1, 2, 3 for the 3 best models
if(model_df$var_index[i]=="equal") varIndx = rep(1,nspp)
if(model_df$var_index[i]=="unequal") varIndx = seq(1,nspp)
if(model_df$var_index[i]=="survey") varIndx = c(rep(1,14), rep(2, 13), rep(3, 8), rep(4,3))

if(model_df$cov[i] == "sst") {
  obs_covar = cov.sst
}
if(model_df$cov[i] == "ssh") {
  obs_covar = cov.ssh
}
if(model_df$cov[i] == "ild") {
  obs_covar = cov.ild
}
if(model_df$cov[i] == "bv") {
  obs_covar = cov.bv
}
if(model_df$cov[i] == "beuti") {
  obs_covar = cov.beuti
}
if(model_df$cov[i] == "cuti") {
  obs_covar = cov.cuti
}

fit.mod = fit_dfa(y = Y,
                  num_trends = model_df$num_trends[i],
                  iter=n_iter,
                  varIndx = varIndx,
                  chains=n_chains, estimate_nu=model_df$est_nu[i],
                  estimate_trend_ma = model_df$estimate_trend_ma[i],
                  estimate_trend_ar = model_df$estimate_trend_ar[i],
                  estimate_process_sigma = model_df$estimate_process_sigma[i],
                  seed=123,
                  obs_covar = obs_covar)

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

saveRDS(fit.mod, file = paste0("results_covariate/",model_df$cov[i],"_",model_df$num_trends[i],
                               "_",model_df$var_index[i],"_",str,theta_str,phi_str,sigma_str,".rds"))





# This is just an optional loop to do the forecasts with a model without covariates

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE, estimate_process_sigma = c(TRUE, FALSE),
                       var_index = c("survey"), num_trends = 1,
                       elpd_kfold = NA, se_elpd_kfold=NA)

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
  
  # We have to set par_list = "all" to enable diagnostics. 
  # Added k-fold cross validation because of unstable Pareto-k
  n_cv = 10 # insted of folds, we'll make predictions for the last 10 years
  log_lik = matrix(0, nrow=n_chains*n_iter/2, ncol = n_cv)
  for(k in 1:n_cv) {
    year_train_max = (max_year - k)
    year_test = max_year - k + 1
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
      pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1)# + pars$b_obs[j,,] * covar_test$value
      log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))
  
  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  saveRDS(model_df, file="results_covariate/covariate_model_summary_null.rds")
  
}

# And as a final check, the full model with all covariates





##### Local Environmental Variables
##sst
sst.value<-scale(avgROMs.short$sst)
cov.sst<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.sst$timeseries<-rep(1:nspp, each=nyear)
cov.sst$time<-rep(seq(1:nyear), times=nspp)
cov.sst$value<-rep(sst.value, times=nspp)
cov.sst$covariate<-as.integer("1")

##ssh
ssh.value<-scale(avgROMs.short$ssh)
cov.ssh<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.ssh$timeseries<-rep(1:nspp, each=nyear)
cov.ssh$time<-rep(seq(1:nyear), times=nspp)
cov.ssh$value<-rep(ssh.value, times=nspp)
cov.ssh$covariate<-as.integer("1")

##ild
ild.value<-scale(avgROMs.short$ild)
cov.ild<-data.frame(time=numeric(n_row),
                    timeseries=numeric(n_row),
                    covariate=numeric(n_row),
                    value=numeric(n_row))
cov.ild$timeseries<-rep(1:nspp, each=nyear)
cov.ild$time<-rep(seq(1:nyear), times=nspp)
cov.ild$value<-rep(ild.value, times=nspp)
cov.ild$covariate<-as.integer("1")

##BV
bv.value<-scale(avgROMs.short$bv)
cov.bv<-data.frame(time=numeric(n_row),
                   timeseries=numeric(n_row),
                   covariate=numeric(n_row),
                   value=numeric(n_row))
cov.bv$timeseries<-rep(1:nspp, each=nyear)
cov.bv$time<-rep(seq(1:nyear), times=nspp)
cov.bv$value<-rep(bv.value, times=nspp)
cov.bv$covariate<-as.integer("1")

##BEUTI
beuti.value<-scale(avgROMs.short$beuti)
cov.beuti<-data.frame(time=numeric(n_row),
                      timeseries=numeric(n_row),
                      covariate=numeric(n_row),
                      value=numeric(n_row))
cov.beuti$timeseries<-rep(1:nspp, each=nyear)
cov.beuti$time<-rep(seq(1:nyear), times=nspp)
cov.beuti$value<-rep(beuti.value, times=nspp)
cov.beuti$covariate<-as.integer("1")

##cuti
cuti.value<-scale(avgROMs.short$cuti)
cov.cuti<-data.frame(time=numeric(n_row),
                     timeseries=numeric(n_row),
                     covariate=numeric(n_row),
                     value=numeric(n_row))
cov.cuti$timeseries<-rep(1:nspp, each=nyear)
cov.cuti$time<-rep(seq(1:nyear), times=nspp)
cov.cuti$value<-rep(cuti.value, times=nspp)
cov.cuti$covariate<-as.integer("1")



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
  
  cov.sst$covariate=1
  cov.beuti$covariate=2
  cov.bv$covariate=3
  cov.cuti$covariate=4
  cov.ild$covariate=5
  cov.ssh$covariate=6
  
  obs_covar = rbind(cov.sst, cov.beuti, cov.bv, cov.cuti, cov.ild, cov.ssh)
  
  #  Have to set par_list = "all" to enable diagnostics. 
  #  Added k-fold cross validation because of unstable Pareto-k
  n_cv = 10 # insted of folds, we'll make predictions for the last 10 years
  log_lik = matrix(0, nrow=n_chains*n_iter/2, ncol = n_cv)
  for(k in 1:n_cv) {
    year_train_max = (max_year - k)
    year_test = max_year - k + 1
    ytrain = Y[,1:(ncol(Y)-k)]
    ytest = Y[,(ncol(Y)-k+1)]
    
    #for(p in 1:6) {
    #  obs_covar$value[which(obs_covar$covariate==p)] = (obs_covar$value[which(obs_covar$covariate==p)] - mean(obs_covar$value[which(obs_covar$covariate==p)])) / sd(obs_covar$value[which(obs_covar$covariate==p)])
    #}
    #obs_covar = dplyr::group_by(obs_covar, covariate, timeseries) %>%
    #  dplyr::mutate(new_value = scale(value)) %>%
    #  dplyr::select(-value) %>%
    #  dplyr::rename(value=new_value)
    covar_train = dplyr::filter(obs_covar, time <= (year_train_max - min(avgROMs$year) + 1))
    covar_test = dplyr::filter(obs_covar, time==(year_test - min(avgROMs$year) + 1))
    covar_pred = dplyr::filter(obs_covar, time==(year_test - min(avgROMs$year) + 1), timeseries==1)
    
    fit.mod = fit_dfa(y = ytrain,
                      num_trends = model_df$num_trends[i],
                      iter=n_iter,
                      varIndx = varIndx,
                      chains=n_chains, estimate_nu=model_df$est_nu[i],
                      estimate_trend_ma = model_df$estimate_trend_ma[i],
                      estimate_trend_ar = model_df$estimate_trend_ar[i],
                      estimate_process_sigma = model_df$estimate_process_sigma[i],
                      seed=123,
                      obs_covar = covar_train)
    # extract log_lik for new data
    pars = rstan::extract(fit.mod$model)
    
    for(j in 1:nrow(log_lik)) {
      # loop over iterations
      pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1) + t(pars$b_obs[j,,]) %*% matrix(covar_pred$value,ncol=1)
      log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))
  
  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  saveRDS(model_df, file="results_covariate/covariate_model_summary_full.rds")
  
}
