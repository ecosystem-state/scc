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
# add in the dataset
Y$dataset = substr(as.character(Y$code),1,4)

dataset_combos1 = data.frame(x1 = unique(Y$dataset), x2=NA,x3=NA,x4=NA)
dataset_combos2 = expand.grid(x1 = unique(Y$dataset), x2 = unique(Y$dataset),x3=NA,x4=NA)
dataset_combos2 = dplyr::filter(dataset_combos2, x1!=x2)
dataset_combos3 = expand.grid(x1 = unique(Y$dataset), x2 = unique(Y$dataset),
                              x3 = unique(Y$dataset), x4=NA)
dataset_combos3 = dplyr::filter(dataset_combos3, x1!=x2, x1!=x3, x2!=x3)
dataset_combos4 = data.frame(x1 = "calc",x2="RREA",x3="SBRD",x4="ZALO")
dataset_combos = rbind(dataset_combos1, dataset_combos2, dataset_combos3, dataset_combos4)
dataset_combos$id = seq(1,nrow(dataset_combos))

model_df = expand.grid(estimate_trend_ma = c(TRUE,FALSE),
                       estimate_trend_ar = c(TRUE,FALSE), 
                       est_nu = c(TRUE,FALSE), 
                       estimate_process_sigma = c(TRUE,FALSE),
                       var_index = c("equal","survey"), 
                       num_trends = 1:3,
                       elpd_kfold = NA, se_elpd_kfold=NA,
                       trend_model = "rw",n_knots = NA,
                       id = dataset_combos$id)
# also join in the b-spline models 
model_df2 = expand.grid(estimate_trend_ma = FALSE,
                        estimate_trend_ar = FALSE,
                        est_nu = FALSE,
                        estimate_process_sigma = FALSE,
                        var_index = c("equal","survey"),
                        num_trends = 1:3,
                        elpd_kfold = NA, se_elpd_kfold=NA,
                        trend_model = "spline",
                        n_knots = c("1/3", "2/3", "nyear"),
                        id = dataset_combos$id)
model_df = rbind(model_df, model_df2)

# load environmental covariate data, and average over spring months to have one value per covar per year
nspp=dim(Y)[1]

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

n_cv = 10 # insted of folds, we'll make predictions for the last 10 years
output_samples = 3000
n_chains = 1

for(i in 1:nrow(model_df)) {
  
  if(model_df$var_index[i]=="equal") varIndx = rep(1,nspp)
  #if(model_df$var_index[i]=="unequal") varIndx = seq(1,nspp)
  if(model_df$var_index[i]=="survey") varIndx = as.numeric(as.factor(Y$dataset))
  
  log_lik = matrix(0, nrow=output_samples, ncol = n_cv)
  for(k in 1:n_cv) {
    ysub = Y[,which(names(Y)%in%c("code","dataset")==FALSE)]
    ysub = ysub[,1:(ncol(ysub)-k+1)]    # get subset of years
    ysub = t(as.matrix(scale(t(ysub)))) #standardized
    # make years we want to predict on NA
    ytrain = ysub
    ytrain[,ncol(ysub)] = NA
    ytest = ysub[,ncol(ysub)]

    n_knots = NULL
    if(model_df$trend_model[i]=="spline") {
      if(model_df$n_knots[i] == "nyear") n_knots = ncol(ytrain) 
      if(model_df$n_knots[i] == "1/3") n_knots = round(3*ncol(ytrain)/3)
      if(model_df$n_knots[i] == "2/3") n_knots = round(2*ncol(ytrain)/3)
    }
    
    fit.mod = try(fit_dfa(y = ytrain,
                          num_trends = model_df$num_trends[i],
                          #iter=n_iter,
                          varIndx = varIndx,
                          scale="none",
                          chains=n_chains, estimate_nu=model_df$est_nu[i],
                          estimate_trend_ma = model_df$estimate_trend_ma[i],
                          estimate_trend_ar = model_df$estimate_trend_ar[i],
                          estimate_process_sigma = model_df$estimate_process_sigma[i],
                          trend_model = as.character(model_df$trend_model[i]),
                          n_knots = n_knots,
                          seed=123,
                          estimation = "vb",
                          #iter=output_samples,
                          tol_rel_obj = 0.005,
                          iter=100000,
                          output_samples = output_samples,
                          #importance_resampling = TRUE,
                          adapt_iter = 100,
                          refresh = 0
    ),silent=TRUE)
    
    if(class(fit.mod) != "try-error") {
      # extract log_lik for new data
      pars = rstan::extract(fit.mod$model)
      
      # loop over iterations
      for(j in 1:nrow(log_lik)) {
        #pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1)
        pred = pars$Z[j,,] %*% matrix(pars$x[j,,ncol(pars$x)],ncol=1)
        log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
      }
    }
  }
  elpds = apply(log_lik,2,log_sum_exp)
  model_df$elpd_kfold[i] = sum(elpds)
  model_df$se_elpd_kfold[i] = sqrt(length(elpds) * var(elpds))
  
  # Also fit the full model
  n_knots = NULL
  if(model_df$trend_model[i]=="spline") {
    if(model_df$n_knots[i] == "nyear") n_knots = ncol(ytrain) 
    if(model_df$n_knots[i] == "1/3") n_knots = round(2*ncol(ytrain)/3)
    if(model_df$n_knots[i] == "2/3") n_knots = round(1*ncol(ytrain)/3)
  } 
  
  ysub = Y[,which(names(Y)%in%c("code","dataset")==FALSE)]
  ysub = t(as.matrix(scale(t(ysub)))) #standardized
  fit.mod = try(fit_dfa(y = ysub,
                        num_trends = model_df$num_trends[i],
                        #iter=n_iter,
                        varIndx = varIndx,
                        scale="none",
                        chains=n_chains, estimate_nu=model_df$est_nu[i],
                        estimate_trend_ma = model_df$estimate_trend_ma[i],
                        estimate_trend_ar = model_df$estimate_trend_ar[i],
                        estimate_process_sigma = model_df$estimate_process_sigma[i],
                        trend_model = as.character(model_df$trend_model[i]),
                        n_knots = n_knots,
                        seed=123,
                        estimation = "vb",
                        tol_rel_obj = 0.005,
                        #iter=output_samples,
                        iter=100000,
                        output_samples = output_samples,
                        #importance_resampling = TRUE,
                        adapt_iter = 100,
                        refresh = 0
  ),silent=TRUE)
  
  if(class(fit.mod) != "try-error") {
    rot = rotate_trends(fit.mod)
    loadings = dfa_loadings(rot)
    fitted = dfa_fitted(fit.mod)
    trends = dfa_trends(rot)
    # don't include model or samples
    fit.mod$model <- NULL
    fit.mod$samples <- NULL
    fit.mod$samples_permuted <- NULL
    
    saveRDS(list(model=fit.mod,fitted=fitted,trends=trends,loadings=loadings), file = paste0("results_biology/biology_",i,".rds"))
  }
  
}

saveRDS(model_df, file="results_biology/biology_model_summary_vb.rds")
