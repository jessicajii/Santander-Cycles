# Title:   Hakwes process model with spatial component
# Author:  Francesco Sanna Passino

## Empty list object
weekly_data = list()
weekly_data_test = list()
## List all the .csv file names (and sort them in case they are not ordered by number)
file_names = sort(list.files('~/Documents/github/Santander-Cycles/santander_bikes_data/', pattern='.csv'))
## Total number of files
n_weeks = length(file_names)
## Import 4 weeks of data as training set
for(week in (n_weeks-5):(n_weeks-2)){
  weekly_data[[week]] = read.table(paste('~/Documents/github/Santander-Cycles/santander_bikes_data/',file_names[week],sep=''), 
                          sep=',', header=FALSE, 
                          col.names=c('start_id','end_id','start_time','duration'))
}
## Last 2 weeks as test set
for(week in (n_weeks-1):n_weeks){
  weekly_data_test[[week]] = read.table(paste('~/Documents/github/Santander-Cycles/santander_bikes_data/',file_names[week],sep=''), 
                                   sep=',', header=FALSE, 
                                   col.names=c('start_id','end_id','start_time','duration'))
}

## Create a unique dataframe
library(dplyr)
df = dplyr::bind_rows(weekly_data)
df = transform(df, end_time = start_time + duration)
df_test = dplyr::bind_rows(weekly_data)
df_test = transform(df, end_time = start_time + duration)
## Import stations
stations = read.table('~/Documents/github/Santander-Cycles/santander_locations.csv', sep=',', header=TRUE)

## Function to calculate C_sequence
##    - t is the vector of reference times (Hyde Park starting times)
##    - starting_times is a list containing all the starting times for the other stations
C_sequence_calculator = function(t, starting_times){
  C_sequence = list()
  for(j in 1:length(starting_times)){
    cat('Progress:', 100*round(j/length(starting_times),4), "%\r") 
    flush.console()
    C_sequence[[j]] = list() 
    t_prime = starting_times[[j]]
    index = ifelse(t_prime[1] < t[1], max(which(t_prime < t[1])), 0)
    if(index != 0){
      C_sequence[[j]][[1]] = t[1] - t_prime[1:index]
    } else {
      C_sequence[[j]][[1]] = numeric(0)
    }
    for(i in 2:length(t)){
      C_sequence[[j]][[i]] = t[i] - t_prime[t_prime < t[i] & t_prime >= t[i-1]]
    }
  }
  return(C_sequence)
}

## IDs of stations
IDs = sort(unique(c(df$start_id,df$end_id)))
subset_stations = which(sapply(stations$Station.Id, function(x) any(x == IDs)))
library(geodist)
distances = geodist(stations[subset_stations,c('longitude','latitude')], 
                    measure='vincenty') / 1000
## Hyde Park distances
id = stations[grepl('Hyde Park Corner', stations$StationName),]$Station.Id
hp_dists = distances[which(IDs == id),]

## Weights must be precomputed
kappa = sapply(hp_dists, function(x,q=1){exp(-q*x)})
## kappa = as.numeric(kappa == 1)
## Hyde Park starting times
t_hyde = df$start_time[df$start_id == id]
tmin = min(t_hyde)
tmax = max(t_hyde)
                               
## Preprocessing 
starting_times = list()
for (j in 1:length(IDs)){
  starting_times[[j]] = df$start_time[df$start_id == IDs[j]]
  m1 = min(starting_times[[j]])
  m2 = max(starting_times[[j]])
  if(m1 < tmin){
    tmin = m1 
  }
  if(m2 > tmax){
    tmax = m2
  }
}
t0 = floor(tmin / 60 / 60 / 24) * (60 * 60 * 24)
T = (ceiling(tmax / 60 / 60 / 24) * (60 * 60 * 24) - t0) / 60 / 60 ## T = 672 hours as expected 

## Preprocess the data
for(j in 1:length(IDs)){
  starting_times[[j]] = sort(starting_times[[j]] - t0 + runif(n=length(starting_times[[j]]))) / 60 / 60
}

## Hyde Park starting times
t_hyde = starting_times[[which(IDs == id)]]
C_seq = C_sequence_calculator(t_hyde, starting_times)                                                            

## Define function for tansforming parameters
transform_parameters_spatial = function(params, to_constrained=TRUE){
  n_stations = (length(params) - 1) / 2
  tp = numeric(1 + n_stations * 2)
  if(to_constrained==FALSE){
    tp[1] = log(params[1])
    tp[seq(2,2*n_stations,2)] = log(params[seq(2,2*n_stations,2)])
    tp[seq(3,2*n_stations+1,2)] = log(params[seq(3,2*n_stations+1,2)] - params[seq(2,2*n_stations,2)])
  } else {
    tp[1] = exp(params[1])
    tp[seq(2,2*n_stations,2)] = exp(params[seq(2,2*n_stations,2)])
    tp[seq(3,2*n_stations+1,2)] = exp(params[seq(3,2*n_stations+1,2)]) + exp(params[seq(2,2*n_stations,2)])
  }
  return(tp)
}

# Log-likelihood for the model with spatial component
spatial_likelihood = function(params, t, T, C_sequence, kappa){
  ## Parameters
  trans_params = transform_parameters_spatial(params, to_constrained=TRUE)
  lambda = trans_params[1]
  beta = trans_params[seq(2,2*length(C_sequence),2)]
  theta = trans_params[seq(3,2*length(C_sequence)+1,2)]
  ## Log-likelihood
  logl = -lambda*T
  for(j in 1:length(C_sequence)){
    logl = logl + kappa[j] * (beta[j]/theta[j] * sum(exp(-theta[j] * (T-starting_times[[j]])) - 1))
  }
  C = numeric(length(C_sequence))
  for(j in 1:length(C_sequence)){
    C[j] = sum(exp(-theta[j]*C_sequence[[j]][[1]])) ## IMPORTANT: initial value of C[j] is *not* 0
  }
  for(i in 1:length(t)){
    logl = logl + log(lambda + sum(kappa * beta * C))
    if(i < length(t)){
      for(j in 1:length(C_sequence)){
        C[j] = exp(-theta[j]*(t[i+1]-t[i])) * C[j] + sum(exp(-theta[j]*C_sequence[[j]][[i+1]]))
      }
    }
  }
  return(logl)
}

## Gradients
loglikelihood_gradients = function(params, t, T, C_sequence, kappa){
  ## Parameters
  trans_params = transform_parameters_spatial(params, to_constrained=TRUE)
  lambda = trans_params[1]
  beta = trans_params[seq(2,2*length(C_sequence),2)]
  theta = trans_params[seq(3,2*length(C_sequence)+1,2)]
  ## Obtain process C
  C = matrix(nrow=length(C_sequence),ncol=length(t))
  C_add = matrix(nrow=length(C_sequence),ncol=length(t))
  for(j in 1:length(C_sequence)){
    C[j,1] = sum(exp(-theta[j]*C_sequence[[j]][[1]]))
  }
  for(i in 1:length(t)){
    if(i < length(t)){
      for(j in 1:length(C_sequence)){
        C_add[j,i+1] = sum(exp(-theta[j]*C_sequence[[j]][[i+1]]))
        C[j,i+1] = exp(-theta[j]*(t[i+1]-t[i])) * C[j,i] + C_add[j,i+1]
      }
    }
  }
  ## Initialise the vector of gradients
  gradients = numeric(length(trans_params))
  ## Lambda
  gradients[1] = lambda * (- T + sum(1 / (lambda + rowSums(beta*kappa*C))))
  ## Theta
  exp_theta_star = exp(params[seq(3,2*length(C_sequence)+1,2)])
  exp_diffs = numeric(length(C_sequence))
  n_events = numeric(length(C_sequence))
  for(j in 1:length(C_sequence)){
    exp_diffs[j] = sum(exp(-theta[j] * (T-starting_times[[j]])))
    n_events[j] = length(starting_times[[j]])
  }
  zeta = lambda + rowSums(beta*kappa*C)
  beta_indices = seq(2,2*length(C_sequence),2)
  theta_indices = seq(3,2*length(C_sequence)+1,2)
  grad_C_beta = C[,1] * beta 
  grad_C_theta = C[,1] * exp_theta_star 
  gradients[beta_indices] = kappa * beta / theta * ((1 - beta / theta) * (exp_diffs - n_events) - beta * exp_diffs)
  gradients[theta_indices] = - kappa * beta / theta * exp_theta_star * ((exp_diffs - n_events) / theta + exp_diffs)
  for(h in 1:length(t)){
    gradients[beta_indices] = gradients[beta_indices] - 1 / zeta * sum(kappa * beta * (C[,h] + grad_C_beta))
    gradients[theta_indices] = gradients[theta_indices] - 1 / zeta * sum(kappa * beta * grad_C_theta)
    if(h < length(t)){
      grad_C_beta = exp(-theta*(t[h+1]-t[h])) * (grad_C_beta - C[,h] * beta) - C_add[,h+1] * beta
      grad_C_theta = exp(-theta*(t[h+1]-t[h])) * (grad_C_theta - C[,h] * exp_theta_star) - C_add[,h+1] * exp_theta_star
    }
  }
  return(gradients)
}

## Adam optimiser
adam = function(starting_values, t, T, C_sequence, kappa, eta=0.1, 
                max_iter=100, tol=1e-6, b1=0.9, b2=0.99, eps=1e-2){
  m = numeric(length(starting_values))
  v = numeric(length(starting_values))
  iter = 1
  convergence = FALSE
  params = starting_values
  liks = numeric()
  while(convergence == FALSE & iter <= max_iter){
    g = loglikelihood_gradients(params=params, t=t, T=T, C_sequence=C_sequence, kappa=kappa)
    m = b1 * m + (1-b1) * g
    v = b2 * v + (1-b2) * (g^2)
    m_hat = m / (1-(b1^iter))
    v_hat = v / (1-(b2^iter))
    params = params + eta * m_hat / (sqrt(v_hat) + eps)
    if(iter > 1){
      lik_old = lik
    }
    lik = spatial_likelihood(params=params, t=t, T=T, C_sequence=C_sequence, kappa=kappa)
    cat('Iteration:', iter, '\tLog-likelihood:', lik, '\r')
    flush.console()
    liks = c(liks,lik)
    if(iter > 1){
      convergence = (abs(lik - lik_old) / abs(lik_old) < tol)
    }
    iter = iter + 1
  }
  return(list(par=params,f=liks))
}

## Important: input consists of *constrained* parameters
p_values = function(params, t, T, C_sequence, kappa, starting_times){
  ## Parameters
  lambda = params[1]
  beta = params[seq(2,2*length(C_sequence),2)]
  theta = params[seq(3,2*length(C_sequence)+1,2)]
  C = numeric(length(C_sequence))
  C2 = numeric(length(C_sequence))
  N = numeric(length(C_sequence))
  N2 = numeric(length(C_sequence))
  for(j in 1:length(C_sequence)){
    C[j] = sum(exp(-theta[j]*C_sequence[[j]][[1]])) 
    C2[j] = C[j]
    N[j] = length(which(starting_times[[j]] < t[1]))
    N2[j] = N[j]
  }
  pvals = numeric(length(t)-1)
  for(i in 1:length(t)){
    if(i < length(t)){
      for(j in 1:length(C_sequence)){
        C[j] = exp(-theta[j]*(t[i+1]-t[i])) * C[j] + sum(exp(-theta[j]*C_sequence[[j]][[i+1]]))
        N[j] = length(which(starting_times[[j]] < t[i+1]))
      }
      pvals[i] = -lambda*(t[i+1] - t[i])
      for(j in 1:length(C_sequence)){
        pvals[i] = pvals[i] + beta[j] / theta[j] * kappa[j] * (C[j] - C2[j] - N[j] + N2[j])
      }
    }
    C2 = C
    N2 = N
  }
  return(exp(pvals))
}

## Optimisation
params = rep(0,2*length(C_seq)+1)
para = adam(starting_values=params, t=t_hyde, T=T, C_sequence=C_seq, kappa=kappa, max_iter=50, eta=0.1)

## p-values
pp = transform_parameters_spatial(para$par, to_constrained=TRUE)
pvals = p_values(params=pp, t=t_hyde, T=T, C_sequence=C_seq, kappa=kappa, starting_times=starting_times)
hist(pvals,freq=FALSE,main='Histogram of p-values')
qqplot(pvals, seq(0,length(t_hyde)-1)/(length(t_hyde)-1),type='l',lwd=2,col='red',
       ylab='Theoretical quantiles',xlab='Observed quantiles', main='Q-Q plot')
abline(0,1,lty='dotted')
ks.test(pvals,'punif')
