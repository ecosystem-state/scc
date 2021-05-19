library(rstan)
library(dplyr)
library(reshape2)

#Organize SCC biology data
max_year = 2017
load("data/ewidata.rda")

sub_data<-ewidata[which(ewidata$year>1979&ewidata$year<=max_year),] # use to biology dataset to match time period of climate data
dfa_data = sub_data[which(sub_data$system=="SCC"&sub_data$subtype!="climate"),]

## Double check that time series are the correct ones
unique(dfa_data$code)

# Run if want to apply DFA to RREAS time series only CalCOFI time series only or
#keep <- grep("RREAS.",dfa_data$code)
#dfa_data=dfa_data[keep,]
#dfa_data=dfa_data[,c(1:3)]

## Log transform those data that are skewed.
melted = melt(dfa_data[, c("code", "year", "value")], id.vars = c("code", "year"))
dat <- dcast(melted, year ~ code)

write.csv(names(dat), file="ts_names_covariate.csv")
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
nspp=nrow(Y)
nyear=ncol(Y)
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

# load environmental covariate data, and average over spring months to have one value per covar per year

dat<-read.csv("data/roms_data_south.csv") # calcofi

avgROMs<-dat %>%
  group_by(year) %>%
  summarise(spr.sst=mean(sst[month%in%c(3:5)],na.rm=T), # Was 3:5
    spr.ssh=mean(ssh[month%in%c(3:5)],na.rm=T),
    spr.ild=mean(ild[month%in%c(3:5)],na.rm=T),
    spr.BV=mean(BV[month%in%c(3:5)],na.rm=T),
    spr.CUTI=mean(CUTI[month%in%c(3:5)],na.rm=T),
    spr.BEUTI=mean(BEUTI[month%in%c(3:5)],na.rm=T))
future_roms = dplyr::filter(avgROMs, year==max_year+1)
avgROMs = dplyr::filter(avgROMs,year <= max_year)

colnames(avgROMs)<-c("year","spr.sst","spr.ssh","spr.ild","spr.BV","spr.CUTI","spr.BEUTI")
#write.csv(avgROMs,".csv")

#subset avgROMs for years of interest
avgROMs.short = subset(avgROMs, year>"1979" & year<"2019",select=year:spr.BEUTI)

#create separate dataframes for each env variable
##### Local Environmental Variables
##sst
sst.value<-avgROMs.short$spr.sst
cov.sst<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.sst$timeseries<-rep(1:nspp, each=nyear)
cov.sst$time<-rep(seq(1:nyear), times=nspp)
cov.sst$value<-rep(sst.value, times=nspp)
cov.sst$covariate<-as.integer("1")

##ssh
ssh.value<-avgROMs.short$spr.ssh
cov.ssh<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.ssh$timeseries<-rep(1:nspp, each=nyear)
cov.ssh$time<-rep(seq(1:nyear), times=nspp)
cov.ssh$value<-rep(ssh.value, times=nspp)
cov.ssh$covariate<-as.integer("1")

##ild
ild.value<-avgROMs.short$spr.ild
cov.ild<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.ild$timeseries<-rep(1:nspp, each=nyear)
cov.ild$time<-rep(seq(1:nyear), times=nspp)
cov.ild$value<-rep(ild.value, times=nspp)
cov.ild$covariate<-as.integer("1")

##BV
bv.value<-avgROMs.short$spr.BV
cov.bv<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.bv$timeseries<-rep(1:nspp, each=nyear)
cov.bv$time<-rep(seq(1:nyear), times=nspp)
cov.bv$value<-rep(bv.value, times=nspp)
cov.bv$covariate<-as.integer("1")

##BEUTI
beuti.value<-avgROMs.short$spr.BEUTI
cov.beuti<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.beuti$timeseries<-rep(1:nspp, each=nyear)
cov.beuti$time<-rep(seq(1:nyear), times=nspp)
cov.beuti$value<-rep(beuti.value, times=nspp)
cov.beuti$covariate<-as.integer("1")

##cuti
cuti.value<-avgROMs.short$spr.CUTI
cov.cuti<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov.cuti$timeseries<-rep(1:nspp, each=nyear)
cov.cuti$time<-rep(seq(1:nyear), times=nspp)
cov.cuti$value<-rep(cuti.value, times=nspp)
cov.cuti$covariate<-as.integer("1")

######## Run models
# specs for best model "mod1_t" = fit_dfa(y = Y, num_trends = 1, iter=4000, varIndx = equal,
#chains=3, estimate_nu=FALSE, seed=123)
n_chains = 3
n_iter = 4000

options(mc.cores = parallel::detectCores())

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]
nyear=dim(Y)[2]
ntrend=1
n_row=nspp*nyear
n_row_pro=ntrend*nyear

model_df = expand.grid(
  var_index = c("survey"),
  elpd_kfold = NA, se_elpd_kfold=NA,
  cov = c("sst","ssh","ild","bv","beuti", "cuti")
)

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

  # these are just Mary's same calls with blown out arguments. Note that I have to set par_list = "all"
  # to enable diagnostics. Eric added k-fold cross validation because of unstable Pareto-k
  n_cv = 10 # insted of folds, we'll make predictions for the last 10 years
  log_lik = matrix(0, nrow=n_chains*n_iter/2, ncol = n_cv)
  for(k in 1:n_cv) {
    year_train_max = (max_year - k)
    year_test = max_year - k + 1
    ytrain = Y[,1:(ncol(Y)-k)]
    ytest = Y[,(ncol(Y)-k+1)]
    covar_train = dplyr::filter(obs_covar, time <= (year_train_max - min(avgROMs$year) + 1))
    covar_test = dplyr::filter(obs_covar, time==(year_test - min(avgROMs$year) + 1))

    x = covar_train$value
    y = c(t(ytrain))
    indx = which(!is.na(y))
    x = x[indx]
    y = y[indx]
    time_series = covar_train$timeseries[indx]
    P = length(varIndx)
    n_pos = length(indx)

    fit.mod = stan(file = "multregr.stan",
      data = list("x"=x,"y"=y,"time_series"=time_series,"P"=P,"n_pos"=n_pos,"varIndx" = varIndx),
      chains=n_chains, iter=n_iter)
    # extract log_lik for new data
    pars = rstan::extract(fit.mod)

    for(j in 1:nrow(log_lik)) {
      # loop over iterations
      pred = pars$b[j,] * covar_test$value
      log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))

  #saveRDS(fit.mod, file = paste0("results_biology/biology_",model_df$num_trends[i],
  #  "_",var_index,"_",str,theta_str,phi_str,sigma_str,".rds"))
  saveRDS(model_df, file="results_covariate/multregr_model_summary.rds")

}
