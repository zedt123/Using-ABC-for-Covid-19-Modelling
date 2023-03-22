library(deSolve)
SEIR <- function(t,x, parameters){
  ##Inputs:
  #t: time
  #x: state variables
  #parameters: parameter vector
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  N <- S+E+I+R
  
  #The with statement means we can use the parameter names as in parameters
  with(as.list(parameters),{
    dS <- - S/N * (beta_1 * E + beta_2 * I)
    dE <- S/N * (beta_1 * E + beta_2 * I) - delta * E
    dI <- delta * E - gamma * I
    dR <- gamma * I
    res <- c(dS, dE, dI, dR)
    list(res)
  })
}

run_model <- function(S0, E0, I0, R0, beta_1, beta_2, delta, gamma, t_total){
  out <- lsoda(func = SEIR,
               y = c(S = S0, E = E0, I = I0, R = R0),
               parms = c(beta_1 = beta_1, beta_2 = beta_2, delta = delta, gamma = gamma),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','E','I','R')]))
}

calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}