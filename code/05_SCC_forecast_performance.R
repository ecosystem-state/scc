library(rstan)
library(dplyr)
library(reshape2)

#Organize SCC biology data
max_year = 2018 # last year forecasted in paper
min_year = 1981
load("data/ewidata.rda")

dfa_data <- dplyr::filter(ewidata, year >= min_year, year <= max_year) %>%
  dplyr::filter(system=="SCC", subtype!="climate")

## Double check that time series are the correct ones
unique(dfa_data$code)

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
dat<-read.csv("data/roms.july_june.mean.csv") # the south and central data (July-June are also in the ewidata file
dat=dat[,c(1:7)] 
avgROMs=dat
colnames(avgROMs)<-c("year","sst","ssh","ild","bv","cuti","beuti")
future_roms = dplyr::filter(avgROMs, year==max_year+1)
avgROMs = dplyr::filter(avgROMs,year <= max_year)
#subset avgROMs for years of interest
avgROMs.short = subset(avgROMs, year>=min_year & year<(max_year+1),select=year:beuti)

##### Local Environmental Variables
nspp=nrow(Y)
nyear=ncol(Y)
n_row=nspp*nyear

# Use BEUTI from the best model
value<-avgROMs.short$beuti

cov<-data.frame(time=numeric(n_row),
  timeseries=numeric(n_row),
  covariate=numeric(n_row),
  value=numeric(n_row))
cov$timeseries<-rep(1:nspp, each=nyear)
cov$time<-rep(seq(1:nyear), times=nspp)
cov$value<-rep(value, times=nspp)
cov$covariate<-as.integer("1")

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

model_df = expand.grid(estimate_trend_ma = FALSE,
                       estimate_trend_ar = TRUE, est_nu = TRUE,
                       estimate_process_sigma = c(FALSE),
                       var_index = c("survey"), num_trends = 1,
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

# fit models to predict the last 10 years
n_cv = 10

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

  # Note that I have to set par_list = "all" to enable diagnostics. 
  # Added k-fold cross validation because of unstable Pareto-k
  for(k in 0:n_cv) {
    # commented out to put this above
    ytrain = Y[,1:(ncol(Y)-k)]
    #ytest = Y[,(ncol(Y)-k+1)]

    # # standardize each
    for(j in 1:nrow(Y)) {
       ytrain[j,] = (ytrain[j,] - mean(ytrain[j,],na.rm=T)) / sd(ytrain[j,],na.rm=T)
    }

    # could be more generalized to filted on time as year rather than dim
    cov_train = dplyr::filter(cov, time <= dim(ytrain)[2])
    #cov_train = cov

    # we need to fit this model 2x. first we'll fit it to the
    # full dataset ytrain
    fit.mod = fit_dfa(y = ytrain,
      num_trends = model_df$num_trends[i],
      iter=n_iter,
      varIndx = varIndx,
      thin=10,
      chains=n_chains, estimate_nu=model_df$est_nu[i],
      estimate_trend_ma = model_df$estimate_trend_ma[i],
      estimate_trend_ar = model_df$estimate_trend_ar[i],
      estimate_process_sigma = model_df$estimate_process_sigma[i],
      seed=123,
      obs_covar = cov_train,
      scale="none")

    saveRDS(fit.mod, paste0("results_forecasting/model_",k,".rds"))

    # and second, we need to drop the last column of NAs in ytrain
    # to keep things interpretable, we don't re-standardize
    #ytrain = ytrain[,-ncol(ytrain)]
    ytrain[,ncol(ytrain)] <- NA
    #cov_train = dplyr::filter(cov, time <= dim(ytrain)[2])

    fit.mod = fit_dfa(y = ytrain,
                      num_trends = model_df$num_trends[i],
                      iter=n_iter,
                      varIndx = varIndx,
                      thin=10,
                      chains=n_chains, estimate_nu=model_df$est_nu[i],
                      estimate_trend_ma = model_df$estimate_trend_ma[i],
                      estimate_trend_ar = model_df$estimate_trend_ar[i],
                      estimate_process_sigma = model_df$estimate_process_sigma[i],
                      seed=123,
                      obs_covar = cov_train,
                      scale="none")
    saveRDS(fit.mod, paste0("results_forecasting/model_",k,"_forecast.rds"))

  }
}

df = data.frame("k"=0:10, "x_star_mean"=NA,
                "x_star_lo"=NA,
                "x_star_hi"=NA,
                "x_mean"=NA,
                "x_lo"=NA,
                "x_hi"=NA)
df$last_year = max_year - df$k
df$pred_year = max_year - df$k + 1


# flipping is a little complicated. for each data point, we have the
# full model and the forecast, and we need to flip them to be in the
# same direction. and to interpret across time points, we need to make sure
# everything is in same direction. use year 0 as reference
fit.mod = readRDS(paste0("results_forecasting/","model_",0,".rds"))
rot_ref = rotate_trends(fit.mod)

for(k in 0:n_cv) {
  #
  # load in DFA model using this held out time step
  fit.mod = readRDS(paste0("results_forecasting/","model_",k,".rds"))
  rot_full = rotate_trends(fit.mod)
  trends_full = rot_full$trends
  n_y = length(rot_full$trends_mean)
  if(sum((-rot_full$trends_mean - rot_ref$trends_mean[1:n_y])^2) < sum((rot_full$trends_mean - rot_ref$trends_mean[1:n_y])^2)) {
    trends_full <- -trends_full # flip
  }
  # summarize mean and quantiles
  df$x_mean[k+1] = mean(trends_full[,,dim(trends_full)[3]])
  df$x_lo[k+1] = quantile(trends_full[,,dim(trends_full)[3]],0.025)
  df$x_hi[k+1] = quantile(trends_full[,,dim(trends_full)[3]], 0.975)

  # do the same with the forecasted value for this
  fit.mod = readRDS(paste0("results_forecasting/","model_",k,"_forecast.rds"))
  rot = rotate_trends(fit.mod)
  trends = rot$trends

  if(sum((-rot$trends_mean - rot_ref$trends_mean[1:n_y])^2) < sum((rot$trends_mean - rot_ref$trends_mean[1:n_y])^2)) {
    trends <- -trends # flip
  }
  # when k = 1, we're making prediction on year 37
  #xtt1 = fit.mod$samples_permuted$xstar
  df$x_star_mean[k+1] = mean(trends[,,dim(trends)[3]])
  df$x_star_lo[k+1] = quantile(trends[,,dim(trends)[3]],0.025)
  df$x_star_hi[k+1] = quantile(trends[,,dim(trends)[3]], 0.975)

  #points(length(x_full), df$x_star_mean[k+1],col="blue")
}

# change to long format
df_true = df[,c("k","x_mean","x_lo","x_hi","last_year","pred_year")]
df_true$state = "true"
df_pred = df[,c("k","x_star_mean","x_star_lo","x_star_hi","last_year","pred_year")]
df_pred = dplyr::rename(df_pred, x_mean = x_star_mean,
                        x_lo = x_star_lo,
                        x_hi = x_star_hi)
df_pred$state = "forecast"
df = rbind(df_true, df_pred)

pdf("results_forecasting/1-step ahead predictions.pdf")
ggplot(df, aes(pred_year, x_mean,group=state,col=state,fill=state)) +
  geom_pointrange(aes(ymin=x_lo,ymax=x_hi),position=position_dodge(0.4)) +
  theme_bw() +
  xlab("Year predicted") +
  ylab("State index") +
  scale_color_viridis(discrete=TRUE,end=0.8)
dev.off()


