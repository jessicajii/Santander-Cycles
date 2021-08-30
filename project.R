# Title:        Hawkes process for starting times of bike journeys
# Author:       Francesco Sanna Passino, Jessica Ji

## Empty list object
weekly_data = list()
weekly_data_test = list()
## List all the .csv file names (and sort them in case they are not ordered by number)
file_names = sort(list.files('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/', pattern='.csv'))
file_names = file_names[1:which(grepl('261',file_names))]
## Total number of files
n_weeks = length(file_names)
## Import 6 weeks of data as training set
for(week in (n_weeks-5):n_weeks){
  weekly_data[[week]] = read.table(paste('~/Documents/github/Santander-Cycles-Network/santander_bikes_data/',file_names[week],sep=''), 
                                   sep=',', header=FALSE, 
                                   col.names=c('start_id','end_id','start_time','duration'))
}
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

## Find StationID for Hyde Park Corner
id = stations[grepl('Hyde Park Corner', stations$StationName),]$Station.Id
## 
t1 = df$start_time[df$start_id == id]
t1_prime = df$end_time[df$end_id == id]
## Obtain the start of the observation period
t0 = floor(min(c(min(t1),min(t1_prime))) / 60 / 60 / 24) * (60 * 60 * 24)
## Shift the times for convenience + add random uniform noise (and sort) 
## + transform times to hours (for convenience)
t = sort((t1 - t0) + runif(n=length(t1))) / 60 / 60
t_prime = sort(t1_prime - t0 + runif(n=length(t1_prime))) / 60 / 60
## Total time of observation -- 28 days * 24 hours = 672 hours (as expected -- 4 weeks of data)
T = ceiling(max(c(max(t),max(t_prime))) / 24) * 24
## Train
set.seed(123)
t1_train = df_train$start_time[df_train$start_id == id]
t1_prime_train = df_train$end_time[df_train$end_id == id]
## Obtain the start of the observation period
t0_train = floor(min(c(min(t1_train),min(t1_prime_train))) / 60 / 60 / 24) * (60 * 60 * 24)
## Shift the times for convenience + add random uniform noise (and sort) 
## + transform times to hours (for convenience)
t_train = sort((t1_train - t0_train) + runif(n=length(t1_train))) / 60 / 60
t_prime_train = sort(t1_prime_train - t0_train + runif(n=length(t1_prime_train))) / 60 / 60
## Total time of observation -- 28 days * 24 hours = 672 hours (as expected -- 4 weeks of data)
T_train = ceiling(max(c(max(t_train),max(t_prime_train))) / 24) * 24
##
## Test
set.seed(123)
t1_test = df_test$start_time[df_test$start_id == id]
t1_prime_test = df_test$end_time[df_test$end_id == id]

## Shift the times for convenience + add random uniform noise (and sort) 
## + transform times to hours (for convenience)
t_test = sort((t1_test - t0_train) + runif(n=length(t1_test))) / 60 / 60
t_prime_test = sort(t1_prime_test - t0_train + runif(n=length(t1_prime_test))) / 60 / 60
## Total time of observation -- 28 days * 24 hours = 672 hours (as expected -- 4 weeks of data)
T_test = ceiling(max(c(max(t_test),max(t_prime_test))) / 24) * 24
##

## *Negative* log-lihelihood
log_likelihood = function(params, t, t_prime, T, B_sequence){
  ## Reparametrise
  trans_params = transform_parameters(params, inverse=FALSE)
  lambda = trans_params[1]
  beta = trans_params[2]
  theta = trans_params[3]
  alpha = trans_params[4]
  delta = trans_params[5]
  ## B_function is a list with len(t) elements 
  ## B_sequence[[i]] contains the difference between t[i] and the t_primes such that:
  ##    t_prime < t[i] AND t_prime >= t[i-1]
  ## If B_sequence is not passed to the function, calculate it. 
  ## For optimisation, a B_sequence ** must ** be passed to the function for efficient optimisation. 
  if(missing(B_sequence) == TRUE){
    B_sequence = list()
    index = ifelse(t_prime[1] < t[1], max(which(t_prime < t[1])), 0)
    if(index == 0){
      B_sequence[[1]] = t[1] - t_prime[1:index]
    } else {
      B_sequence[[1]] = numeric(0)
    }
    for(i in 2:length(t)){
      B_sequence[[i]] = t[i] - t_prime[t_prime < t[i] & t_prime >= t[i-1]]
    }
  }
  ## Calculate likelihood 
  logl = -lambda*T + beta/theta * sum(exp(-theta * (T-t)) - 1)
  logl = logl + alpha/delta * sum(exp(-delta * (T-t_prime)) - 1)
  A = 0
  B = sum(exp(-delta*B_sequence[[1]])) ## IMPORTANT: initial value of B is *not* 0
  for(i in 1:length(t)){
    logl = logl + log(lambda + beta * A + alpha * B)
    if(i < length(t)){
      A = exp(-theta*(t[i+1]-t[i])) * (1 + A)
      B = exp(-delta*(t[i+1]-t[i])) * B + sum(exp(-delta*B_sequence[[i+1]]))
    }
  }
  return(-logl)
}

## Function to calculate B_sequence
B_sequence_calculator = function(t, t_prime){
  B_sequence = list()
  index = ifelse(t_prime[1] < t[1], max(which(t_prime < t[1])), 0)
  if(index == 0){
    B_sequence[[1]] = t[1] - t_prime[1:index]
  } else {
    B_sequence[[1]] = numeric(0)
  }
  for(i in 2:length(t)){
    B_sequence[[i]] = t[i] - t_prime[t_prime < t[i] & t_prime >= t[i-1]]
  }
  return(B_sequence)
}

## Optimise the log-likelihood
optim_logl = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=x, t=t, t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}

## Load some utility functions
source('~/Desktop/hawkes_double.R')
## If in hawkes_double.R you change test = FALSE to test = TRUE, you will get some 
## simulated data with parameters specified in the file hawkes_double.R 
## The event times are in t and t_prime ** which would replace the Hyde Park event times! 

## Calculate B_sequence
params = c(0.1, 0.25, 1.0, 0.1, 0.25)
B_seq = B_sequence_calculator(t, t_prime)
B_seq_train = B_sequence_calculator(t_train, t_prime_train)
B_seq_full = B_seq_train
index = ifelse(t_prime_test[1] < t_test[1], max(which(t_prime_test < t_test[1])), 0)
if(index == 0){
  B_seq_full[[1+length(t_train)]] = t_test[1] - t_prime_test[1:index]
} else {
  B_seq_full[[1+length(t_train)]] = numeric(0)
}
for(i in 2:length(t_test)){
  B_seq_full[[i+length(t_train)]] = t_test[i] - t_prime_test[t_prime_test < t_test[i] & t_prime_test >= t_test[i-1]]
}

## MLE
res = optim_logl(starting_values=transform_parameters(params, inverse=FALSE),
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p = transform_parameters(res$par, inverse=FALSE)

## Constrained optimisation
optim_logl_constrained = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x,-Inf,0), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}
res_model1 = optim_logl_constrained(starting_values=transform_parameters(params, inverse=FALSE)[1:3],
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p2 = transform_parameters(res_model1$par, inverse=FALSE)
p2[4] = 0
p2[5] = 1

## Important to check the dependence on starting values, better results here!!!
res = optim_logl(starting_values=rep(0,5),
                 t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p = transform_parameters(res$par, inverse=FALSE)
#
res_train = optim_logl(starting_values = rep(0,5),
                       t=t_train, t_prime = t_prime_train, T = T_train, B_sequence = B_seq_train)
p_train = transform_parameters(res_train$par, inverse = FALSE)

#####
## If simulated data are used, something very close to the 'truth' is obtained here 
## params = c(0.1, 0.25, 1.0, 0.1, 0.25) was used in the simulation 

res_model1 = optim_logl_constrained(starting_values=c(0,0,0),
                                    t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p2 = transform_parameters(res_model1$par, inverse=FALSE)

## Other model: only consider the arrival times of the journeys ending at station i
optim_logl_constrained2 = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x[1],-Inf,0,x[2],x[3]), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}
res_model2 = optim_logl_constrained2(starting_values=c(0,0,0),
                                    t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
p3 = transform_parameters(res_model2$par, inverse=FALSE)[c(1,4,5,2,3)]
p3[2] = 0
p3[3] = 1

#
optim_poisson = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x,-Inf,0,-Inf,0), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl, method = "BFGS"))
}
p_poisson = optim_poisson(starting_values = 0, t=t, t_prime = t_prime, T=T, B_sequence = B_seq) 
p4 = transform_parameters(p_poisson$par, inverse = FALSE)
p4 = c(p4[1], 0, 1, 0, 1)

#####
## p-values
p_ik2 = function(params, t, t_prime, T, B_sequence){
  lambda = params[1]
  beta =  params[2]
  theta =  params[3]
  alpha =  params[4]
  delta =  params[5]
  if(missing(B_sequence) == TRUE){
    B_sequence = list()
    index = ifelse(t_prime[1] < t[1], max(which(t_prime < t[1])), 0)
    if(index == 0){
      B_sequence[[1]] = t[1] - t_prime[1:index]
    } else {
      B_sequence[[1]] = numeric(0)
    }
    for(i in 2:length(t)){
      B_sequence[[i]] = t[i] - t_prime[t_prime < t[i] & t_prime >= t[i-1]]
    }
  }
  A = 0
  A2 = 0
  B = sum(exp(-delta*B_sequence[[1]])) 
  B2 = sum(exp(-delta*B_sequence[[1]])) 
  p = numeric(length(t))
  for(k in 2:length(t)){
    if(k < length(t)){
      A = exp(-theta*(t[k]-t[k-1])) * (1 + A)
      B = exp(-delta*(t[k]-t[k-1])) * B + sum(exp(-delta*B_sequence[[k]]))
    }
    p[k] = exp(-lambda*(t[k] - t[k-1]) 
               + (beta/theta)*(A - A2 - 1) 
               + (alpha/delta)*(B - B2 + length(which(t_prime < t[k-1])) - length(which(t_prime < t[k]))))
    A2 = A
    B2 = B
  }
  return(p)
}
p_value1 <- p_ik2(p, t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
hist(p_value1,freq=FALSE,main='Histogram of p-values')
qqplot(p_value1, seq(0,length(t)-1)/(length(t)-1),type='l',lwd=2,col='red',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot')
abline(0,1,lty='dotted')
ks_result <-ks.test(p_value1,'punif') ## Test for uniformity: fail to reject null H0 --> good model fit
### Test
p_value_full <- p_ik2(p_train, t=c(t_train, t_test), t_prime=c(t_prime_train,t_prime_test), T=T_test, B_sequence=B_seq_full)
p_value_train <- p_value_full[1:length(t_train)]
p_value_test <- p_value_full[(length(t_train)+1):(length(t_train)+length(t_test))]
hist(p_value_full,freq=FALSE,main='Histogram of p-values')
hist(p_value_train,freq=FALSE,main='Histogram of p-values')
qqplot(p_value_train, seq(0,length(t)-1)/(length(t)-1),type='l',lwd=2,col='red',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot')
abline(0,1,lty='dotted')
ks_result_train <-ks.test(p_value_train,'punif')
hist(p_value_test,freq=FALSE,main='Histogram of p-values')
qqplot(p_value_test, seq(0,length(t)-1)/(length(t)-1),type='l',lwd=2,col='red',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot')
abline(0,1,lty='dotted')
ks_result_test <-ks.test(p_value_test,'punif')
#
p_value2 <- p_ik2(p3, t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
hist(p_value2,freq=FALSE,main='Histogram of p-values')
ks_result2 <- ks.test(p_value2,'punif')
#
p_value3 <- p_ik2(p2, t=t, t_prime=t_prime, T=T, B_sequence=B_seq)
hist(p_value3,freq=FALSE,main='Histogram of p-values')
ks_result3 <- ks.test(p_value3,'punif')
##
start_ID <- sort(unique(df$start_id))
D_values <- function(i){
  set.seed(123)
  t=filter_all_station(id=start_ID[i],df=df)$t
  t_prime=filter_all_station(id=start_ID[i],df=df)$t_prime
  T=filter_all_station(id=start_ID[i],df=df)$T
  B_sequence=B_sequence_calculator(filter_all_station(id=start_ID[i],df=df)$t, 
                                   filter_all_station(id=start_ID[i],df=df)$t_prime)
  res = optim_logl(starting_values=rep(0,5),
      t=t, t_prime=t_prime, T=T, B_sequence=B_sequence)
  p = transform_parameters(res$par, inverse=FALSE)
  p_values = p_ik2(p, t=t, t_prime=t_prime, T=T, B_sequence=B_sequence)
  ks_results = ks.test(p_values, 'punif')
  D = ks_results$statistic
  result = c(start_ID[i], D, ks_results$p.value, p, length(t), length(t_prime))
  return(result)
}
D_val1 <- sapply(1:300, D_values)
D_val1 <- t(D_val1)
D_val2 <- sapply(301:500, D_values)
D_val2 <- t(D_val2)
D_val3 <- sapply(501:length(start_ID), D_values)
D_val3 <- t(D_val3)
df2 <- rbind(D_val1, D_val2)
df2 <- rbind(df2, D_val3)
colnames(df2) <- c("Station.Id", "D.value", "KS.pvalue","lambda","beta", "theta", "alpha","delta", "n","n_prime")
write.csv(df2,"~/Desktop/Result2.csv", row.names = TRUE)
