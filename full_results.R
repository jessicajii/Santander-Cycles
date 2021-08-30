# Title:        Utility functions for Hawkes processes
# Author:       Francesco Sanna Passino, Jessica Ji

## Empty list object
weekly_data = list()
## List all the .csv file names (and sort them in case they are not ordered by number)
file_names = sort(list.files('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/', pattern='.csv'))
file_names = file_names[1:which(grepl('261',file_names))]
## Total number of files
n_weeks = length(file_names)

## Training set
weekly_data_train = list()
for(week in (n_weeks-5):(n_weeks-2)){
  weekly_data_train[[week]] = read.table(paste('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/',file_names[week],sep=''), 
                                         sep=',', header=FALSE, col.names=c('start_id','end_id','start_time','duration'))
}
df_train = dplyr::bind_rows(weekly_data_train)
df_train = transform(df_train, end_time = start_time + duration)

## Test set
weekly_data_test = list()
for(week in (n_weeks-1):n_weeks){
  weekly_data_test[[week]] = read.table(paste('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/',file_names[week],sep=''), 
                                        sep=',', header=FALSE, col.names=c('start_id','end_id','start_time','duration'))
}
df_test = dplyr::bind_rows(weekly_data_test)
df_test = transform(df_test, end_time = start_time + duration)

## Import stations
stations = read.table('~/Documents/github/Santander-Cycles-Network/santander_locations.csv', sep=',', header=TRUE)

## Load some utility functions
source('~/Documents/github/Santander-Cycles-Network/hawkes_utility.R')

## Filter **all** times (not separately per station)
set.seed(123)
## Obtain the start of the observation period
t0 = floor(min(c(min(df_train$start_time),min(df_train$end_time))) / 60 / 60 / 24) * (60 * 60 * 24)
## Shift the times for convenience + add random uniform noise (and sort) 
## + transform times to hours (for convenience)
df_train = transform(df_train, 
                     start_time_hours = ((df_train$start_time - t0) + runif(n=dim(df_train)[1])) / 60 / 60, 
                     end_time_hours = ((df_train$end_time - t0) + runif(n=dim(df_train)[1])) / 60 / 60)
T = ceiling(max(c(max(df_train$start_time_hours),max(df_train$end_time_hours))) / 24) * 24
df_test = transform(df_test, 
                     start_time_hours = ((df_test$start_time - t0) + runif(n=dim(df_test)[1])) / 60 / 60, 
                     end_time_hours = ((df_test$end_time - t0) + runif(n=dim(df_test)[1])) / 60 / 60)
T_test = ceiling(max(c(max(df_test$start_time_hours),max(df_test$end_time_hours))) / 24) * 24

## Allocate lists
para = list()
pv_poisson_train = list(); pv_poisson_test = list(); ks_pois_train = list(); ks_pois_test = list()
pv_mod1_train = list(); pv_mod1_test = list(); ks_mod1_train = list(); ks_mod1_test = list()
pv_mod2_train = list(); pv_mod2_test = list(); ks_mod2_train = list(); ks_mod2_test = list()
pv_full_train = list(); pv_full_test = list(); ks_full_train = list(); ks_full_test = list()

## Unique StartIDs
start_ids = sort(unique(df_train$start_id))
start_ids_prime = sort(unique(df_test$start_id))

## Repeat for all source stations
set.seed(1771)
for(id in start_ids){
  message('Evaluating station ', id, ' (', stations$StationName[which(stations$Station.Id == id)],
          ')\t\t\t\t\t\t\t\t\t','\r',appendLF=FALSE)
  flush.console()
  ## Obtain filtered times
  t = sort(df_train$start_time_hours[df_train$start_id == id])
  t_prime = sort(df_train$end_time_hours[df_train$end_id == id])
  t_test = sort(df_test$start_time_hours[df_test$start_id == id])
  t_prime_test = sort(df_test$end_time_hours[df_test$end_id == id])
  ## Obtain B_sequence
  B_sequence_train = B_sequence_calculator(t,t_prime)
  B_sequence_full = B_sequence_train
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    index = ifelse(t_prime_test[1] < t_test[1], max(which(t_prime_test < t_test[1])), 0)
    if(index == 0){
      B_sequence_full[[1+length(t)]] = t_test[1] - t_prime_test[1:index]
    } else {
      B_sequence_full[[1+length(t)]] = numeric(0)
    }
    for(i in 2:length(t_test)){
      B_sequence_full[[i+length(t)]] = t_test[i] - t_prime_test[t_prime_test < t_test[i] & t_prime_test >= t_test[i-1]]
    }
  }
  ### Parameter estimation
  ## Poisson
  res_poisson = optim_poisson(starting_values=0, t=t, t_prime=t_prime, T=T, B_sequence=B_sequence_train)
  p_poisson = transform_parameters(res_poisson$par, to_constrained=TRUE)
  p_poisson[c(2,4)] = 0; p_poisson[c(3,5)] = 1
  ## Model 1
  res_mod1 = optim_logl_constrained(starting_values=rep(0,3), t=t, t_prime=t_prime, T=T, B_sequence=B_sequence_train)
  p_mod1 = transform_parameters(res_mod1$par, to_constrained=TRUE)
  p_mod1[4] = 0; p_mod1[5] = 1
  ## Model 2
  res_mod2 = optim_logl_constrained2(starting_values=rep(0,3), t=t, t_prime=t_prime, T=T, B_sequence=B_sequence_train)
  p_mod2 = transform_parameters(res_mod2$par, to_constrained=TRUE)[c(1,4,5,2,3)]
  p_mod2[2] = 0; p_mod2[3] = 1
  ## Full model
  res_full = optim_logl(starting_values=rep(0,5), t=t, t_prime=t_prime, T=T, B_sequence=B_sequence_train)
  p_full = transform_parameters(res_full$par, to_constrained=TRUE)
  ## p-values - Poisson
  pv_poisson = p_value_calculator(p_poisson, t=c(t, t_test), t_prime=c(t_prime,t_prime_test), 
                                       B_sequence=B_sequence_full)
  pv_poisson_train[[id]] = pv_poisson[1:length(t)]
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    pv_poisson_test[[id]] = pv_poisson[(length(t)+1):(length(t)+length(t_test))]
  }
  ## p-values: model 1
  pv_mod1 = p_value_calculator(p_mod1, t=c(t, t_test), t_prime=c(t_prime,t_prime_test), 
                                  B_sequence=B_sequence_full)
  pv_mod1_train[[id]] = pv_mod1[1:length(t)]
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    pv_mod1_test[[id]] = pv_mod1[(length(t)+1):(length(t)+length(t_test))]
  }
  ## p-values: model 1
  pv_mod2 = p_value_calculator(p_mod2, t=c(t, t_test), t_prime=c(t_prime,t_prime_test), 
                               B_sequence=B_sequence_full)
  pv_mod2_train[[id]] = pv_mod2[1:length(t)]
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    pv_mod2_test[[id]] = pv_mod2[(length(t)+1):(length(t)+length(t_test))]
  }
  ## p-values: full model
  pv_full = p_value_calculator(p_full, t=c(t, t_test), t_prime=c(t_prime,t_prime_test), 
                               B_sequence=B_sequence_full)
  pv_full_train[[id]] = pv_full[1:length(t)]
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    pv_full_test[[id]] = pv_full[(length(t)+1):(length(t)+length(t_test))]
  }
  ## Store parameter values
  para[[id]] = cbind(p_poisson, p_mod1, p_mod2, p_full)
  ## KS tests
  ks_pois_train[id] = ks.test(pv_poisson_train[[id]],'punif')$statistic
  ks_mod1_train[id] = ks.test(pv_mod1_train[[id]],'punif')$statistic
  ks_mod2_train[id] = ks.test(pv_mod2_train[[id]],'punif')$statistic
  ks_full_train[id] = ks.test(pv_full_train[[id]],'punif')$statistic
  if(length(t_test) > 0 & length(t_prime_test) > 0){
    ks_pois_test[id] = ks.test(pv_poisson_test[[id]],'punif')$statistic
    ks_mod1_test[id] = ks.test(pv_mod1_test[[id]],'punif')$statistic
    ks_mod2_test[id] = ks.test(pv_mod2_test[[id]],'punif')$statistic
    ks_full_test[id] = ks.test(pv_full_test[[id]],'punif')$statistic
  }
}

## Save output
save(para, pv_poisson_train, pv_poisson_test, ks_pois_train, ks_pois_test,
  pv_mod1_train, pv_mod1_test, ks_mod1_train, ks_mod1_test,
  pv_mod2_train, pv_mod2_test, ks_mod2_train, ks_mod2_test,
  pv_full_train, pv_full_test, ks_full_train, ks_full_test,
  file = "full_results.RData")

## QQ plots for a given station
qq_plots_station = function(id){
  ## Station name
  id_name = unlist(strsplit(stations$StationName[stations$Station.Id == id], ","))[1]
  ## Thresholds
  prob_thresh = seq(0,1,length.out=250)
  ## Training set
  plot(quantile(pv_poisson_train[[id]], probs=prob_thresh), prob_thresh,type='l', lwd=2, col='black',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main=paste('Q-Q plot -', id_name,'- Training set'))
  lines(quantile(pv_mod1_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
        col='blue',ylab='',xlab='', main='')
  lines(quantile(pv_mod2_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
        col='red',ylab='',xlab='', main='')
  lines(quantile(pv_full_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
        col='green',ylab='',xlab='', main='')
  legend('topleft', legend=c('Poisson','Model 1', 'Model 2', 'Full model'), col=c('black','blue','red','green'), lwd=2)
  lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)
  ## Test set
  plot(quantile(pv_poisson_test[[id]], probs=prob_thresh), prob_thresh,type='l', lwd=2, col='black',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main=paste('Q-Q plot -', id_name,'- Test set'))
  lines(quantile(pv_mod1_test[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
         col='blue',ylab='',xlab='', main='')
  lines(quantile(pv_mod2_test[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
         col='red',ylab='',xlab='', main='')
  lines(quantile(pv_full_test[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=2,
         col='green',ylab='',xlab='', main='')
  legend('topleft', legend=c('Poisson','Model 1', 'Model 2', 'Full model'), col=c('black','blue','red','green'), lwd=2)
  lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)
}

## Hyde Park Corner
id = stations$Station.Id[which(stations$StationName == 'Hyde Park Corner, Hyde Park')]
qq_plots_station(id)

## Queen's Gate (North)
id = stations$Station.Id[which(grepl("Queen's Gate \\(North\\)",stations$StationName))]
qq_plots_station(id)

## IC
id = stations$Station.Id[which(grepl("Imperial College, Knightsbridge",stations$StationName))]
qq_plots_station(id)

## British Museum, Bloomsbury
id = stations$Station.Id[which(grepl("British Museum, Bloomsbury",stations$StationName))]
qq_plots_station(id)

## Top-10 stations (** INCREASE PLOT SIZE **)
par(mfrow=c(3,2), mar=c(4,5,4,1))
ids = as.integer(names(sort(table(df_train$start_id), decreasing=TRUE)[2:10]))
for(id in ids){
  qq_plots_station(id)
}
par(mfrow=c(1,1), mar=c(4,4,4,1))
## Comment: it appears that the *full* model does a better job at modelling popular stations

## Multiple stations in the same plot
####### Poisson
## Training set
plot(quantile(pv_poisson_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Poisson - Training set')
for(id in start_ids_prime){
  if(length(pv_poisson_train[[id]]) > 10){
    lines(quantile(pv_poisson_train[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

## Test set
plot(quantile(pv_poisson_test[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Poisson - Test set')
for(id in start_ids_prime){
  if(length(pv_poisson_test[[id]]) > 10){
    lines(quantile(pv_poisson_test[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)
####### Model 1
## Training set
id = 1
prob_thresh = seq(0,1,length.out=250)
plot(quantile(pv_mod1_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Model 1 - Training set')
for(id in start_ids_prime){
  if(length(pv_mod1_train[[id]]) > 10){
    lines(quantile(pv_mod1_train[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

## Test set
id = 1
prob_thresh = seq(0,1,length.out=250)
plot(quantile(pv_mod1_test[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Model 1 - Test set')
for(id in start_ids_prime){
  if(length(pv_mod1_test[[id]]) > 10){
    lines(quantile(pv_mod1_test[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

####### Model 2
## Training set
plot(quantile(pv_mod2_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Model 2 - Training set')
for(id in start_ids_prime){
  if(length(pv_mod2_train[[id]]) > 10){
    lines(quantile(pv_mod2_train[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

## Test set
plot(quantile(pv_mod2_test[[id]], probs=prob_thresh,na.rm = TRUE), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Model 2 - Test set')
for(id in start_ids_prime){
  if(length(pv_mod2_test[[id]]) > 10){
    lines(quantile(pv_mod2_test[[id]], probs=prob_thresh, na.rm = TRUE), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

####### Full Model 
## Training set
plot(quantile(pv_full_train[[id]], probs=prob_thresh), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Full Model - Training set')
for(id in start_ids_prime){
  if(length(pv_full_train[[id]]) > 10){
    lines(quantile(pv_full_train[[id]], probs=prob_thresh), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

## Test set
plot(quantile(pv_full_test[[id]], probs=prob_thresh, na.rm = TRUE), prob_thresh,type='l',lwd=1, xlim=c(0,1), ylim=c(0,1),
     col='gray',ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot - Full Model - Test set')
for(id in start_ids_prime){
  if(length(pv_full_test[[id]]) > 10){
    lines(quantile(pv_full_test[[id]], probs=prob_thresh, na.rm = TRUE), prob_thresh,lwd=1,type='l',
          col='gray',ylab='',xlab='', main='')
  }
}
lines(c(0,1),c(0,1),lty='dotted',col='black',lwd=2)

## Boxplots
n = length(unlist(ks_pois_train)); n_test = length(unlist(ks_pois_test)) 
ks_factors = as.factor(c(rep('Poisson',n), rep('Model 1',n), rep('Model 2',n), rep('Full model', n)))
ks_factors_test = as.factor(c(rep('Poisson',n_test), rep('Model 1',n_test), rep('Model 2',n_test), rep('Full model', n_test)))
ks_scores_all = c(unlist(ks_pois_train), unlist(ks_mod1_train), unlist(ks_mod2_train), unlist(ks_full_train))
ks_scores_all_test = c(unlist(ks_pois_test), unlist(ks_mod1_test), unlist(ks_mod2_test), unlist(ks_full_test))
boxplot(ks_scores_all ~ ks_factors, col=c('lightgreen','lightblue','pink','gray'), 
        xlab='Model', ylab='Kolmogorov-Smirnov scores', main='Boxplots of KS scores - Training set', ylim=c(0,0.5))
boxplot(ks_scores_all_test ~ ks_factors_test, col=c('lightgreen','lightblue','pink','gray'), 
        xlab='Model', ylab='Kolmogorov-Smirnov scores', main='Boxplots of KS scores - Test set', ylim=c(0,0.5))

## Number of observations vs. KS-score 
## Shows that the full model is much better than the others
len_train = unlist(lapply(pv_full_train, length))
len_train = len_train[len_train > 0]
par(mfrow=c(2,2), mar=c(4,4,4,1))
plot(len_train, unlist(ks_pois_train), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores', main='Poisson (Training set)')
abline(lm(unlist(ks_pois_train) ~ len_train)$coef)
plot(len_train, unlist(ks_mod1_train), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Firebrick3', main='Model 1 (Training set)')
abline(lm(unlist(ks_mod1_train) ~ len_train)$coef, col='Firebrick3')
plot(len_train, unlist(ks_mod2_train), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Blue3', main='Model 2 (Training set)')
abline(lm(unlist(ks_mod2_train) ~ len_train)$coef, col='Blue3')
plot(len_train, unlist(ks_full_train), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Green4', main='Full model (Training set)')
abline(lm(unlist(ks_full_train) ~ len_train)$coef, col='Green4')

### test set
len_test = unlist(lapply(pv_full_test, length))
len_test = len_test[len_test > 0]
plot(len_test, unlist(ks_pois_test), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores', main='Poisson (Test set)')
abline(lm(unlist(ks_pois_test) ~ len_test)$coef)
plot(len_test, unlist(ks_mod1_test), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Firebrick3', main='Model 1 (Test set)')
abline(lm(unlist(ks_mod1_test) ~ len_test)$coef, col='Firebrick3')
plot(len_test, unlist(ks_mod2_test), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Blue3', main='Model 2 (Test set)')
abline(lm(unlist(ks_mod2_test) ~ len_test)$coef, col='Blue3')
plot(len_test, unlist(ks_full_test), ylim=c(0,0.5), pch=19, cex=.25, 
     xlab='Number of observations', ylab='Kolmogorov-Smirnov scores',col='Green4', main='Full model (Test set)')
abline(lm(unlist(ks_full_test) ~ len_test)$coef, col='Green4')

## Intensity function
set.seed(110)
t_k = runif(n=5, min=0, max=10)
intensity = function(t, t_k, params){
  lambda=params[1];  beta=params[2]; theta=params[3]
  return(lambda + sum(beta * exp(-theta * (t-t_k[t_k < t]))))
} 
hawkes = sapply(seq(0,18,length.out=1000), FUN=function(t) intensity(t, t_k=t_k, params=c(0.1,0.15,0.7)))
plot(seq(0,18,length.out=1000), hawkes, type='l', xlab='Time', ylab='Intensity', ylim=c(0,0.7))
points(t_k, rep(0,5), pch=19, cex=0.5)
abline(h=0, lty='dotted')
abline(v=t_k, lty='dotted')
