# Utility functions for Hawkes processes
# Author:       Francesco Sanna Passino

## Calculate intensity function
#
intensity = function(t,h,t_prime,params,check_t_prime=TRUE){
  ## Parameters
  lambda=params[1];  beta=params[2]; theta=params[3]; alpha=params[4]; delta=params[5]
  ## 
  if(check_t_prime==TRUE){
    t_prime = t_prime[t_prime < t]
  }
  return(lambda + sum(beta * exp(-theta * (t-h))) + sum(alpha * exp(-delta * (t-t_prime))))
}

## Simulate Hawkes process with scaled exponential self- and mutual-excitation function
hawkes_simulator = function(params, t_prime, T, lambda_star){
  ## Parameters
  lambda=params[1];  beta=params[2]; theta=params[3]; alpha=params[4]; delta=params[5]
  ## Initial points
  t_star = 0
  h = numeric(0)
  t_auxiliary = numeric(0)
  index = 0
  ## Simulate Hawkes process
  while(t_star < T){
    ## Propose new arrival time
    t_star = t_star - log(runif(1)) / lambda_star ## Important to pick a good upper bound lambda*
    if(t_star > T){
      break
    }
    if(index < length(t_prime)){
      if(t_star > t_prime[index+1]){
        for(j in (index+1):length(t_prime)){
          condition = (t_star < t_prime[j])
          if(condition == TRUE){
            index = j
            break
          } else {
            t_auxiliary = c(t_auxiliary, t_prime[j])
          }
        }
      }
    }
    lambda_prob = intensity(t_star,h,t_prime=t_auxiliary,
                            params=params,check_t_prime=FALSE) / lambda_star
    accept = sample(c(1,0), size=1, prob=c(lambda_prob, 1-lambda_prob))
    if(accept == 1){
      h = c(h, t_star)
    }
  }
  return(h)
}

## Transform parameters
transform_parameters = function(params, to_constrained=TRUE){
  tp = numeric(5)
  if(to_constrained==FALSE){
    tp[1] = log(params[1])
    tp[2] = log(params[2])
    tp[3] = log(params[3]-params[2])
    tp[4] = log(params[4])
    tp[5] = log(params[5]-params[4])
  } else {
    tp[1] = exp(params[1])
    tp[2] = exp(params[2])
    tp[3] = exp(params[2]) + exp(params[3])
    tp[4] = exp(params[4])
    tp[5] = exp(params[4]) + exp(params[5])
  }
  return(tp)
}

test = FALSE
if(test==TRUE){
  ## Parameters
  params = c(0.1, 0.25, 1.0, 0.1, 0.25)
  T = 10000
  ## Simulate t_prime as a Poisson process
  rho = 0.25
  set.seed(123)
  t_prime = rexp(n=1,rate=rho)
  while(TRUE){
    proposal = tail(t_prime,1)+rexp(n=1,rate=rho)
    if(proposal < T){
      t_prime = c(t_prime, proposal)
    } else {
      break 
    }
  }
  ## Simulate Hawkes process
  t = hawkes_simulator(params, t_prime, T, lambda_star=10)
}

## Negative log-lihelihood
log_likelihood = function(params, t, t_prime, T, B_sequence){
  ## Reparametrise
  trans_params = transform_parameters(params, to_constrained=TRUE)
  lambda = trans_params[1] #exp(params[1])
  beta = trans_params[2]   #exp(params[2])
  theta = trans_params[3]  #exp(params[2]) + exp(params[3])
  alpha = trans_params[4]  #exp(params[4])
  delta = trans_params[5]  #exp(params[4]) + exp(params[5])
  ## B_function is a list with length(t) elements 
  ## B_sequence[[i]] contains the difference between t[i] and the t_primes such that:
  ##    t_prime < t[i] AND t_prime >= t[i-1]
  ## If B_sequence is not passed to the function, calculate it. 
  ## For efficient optimisation, a B_sequence must be passed to the function. 
  if(missing(B_sequence) == TRUE){
    B_sequence = B_sequence_calculator(t, t_prime)
  }
  ## Calculate likelihood 
  logl = -lambda*T + beta/theta * sum(exp(-theta * (T-t)) - 1)
  logl = logl + alpha/delta * sum(exp(-delta * (T-t_prime)) - 1)
  A = 0
  B = sum(exp(-delta*B_sequence[[1]])) 
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

## Likelihood optimisation
## Simple Poisson model
optim_poisson = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x,-Inf,0,-Inf,0), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl, method = "BFGS"))
}

## Model 1
optim_logl_constrained = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x,-Inf,0), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}

## Model 2
optim_logl_constrained2 = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=c(x[1],-Inf,0,x[2],x[3]), t=t, 
                                        t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}

## Unconstrained (full model)
optim_logl = function(starting_values,t,t_prime,T,B_sequence){
  neg_logl = function(x) log_likelihood(params=x, t=t, t_prime=t_prime, T=T, B_sequence=B_sequence)
  return(optim(par=starting_values, fn=neg_logl))
}

## Calculate p-values
p_value_calculator = function(params, t, t_prime, B_sequence){
  lambda = params[1]
  beta =  params[2]
  theta =  params[3]
  alpha =  params[4]
  delta =  params[5]
  if(missing(B_sequence) == TRUE){
    B_sequence = B_sequence_calculator(t, t_prime)
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
