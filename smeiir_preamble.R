library(deSolve)
SMEIIR <- function(t,x, parameters){
  ##Inputs:
  #t: time
  #x: state variables
  #parameters: parameter vector
  S <- x[1]
  M <- x[2]
  E <- x[3]
  Is <- x[4]
  Ia <- x[5]
  R <- x[6]
  N <- S+M+E+Is+Ia+R
  
  #The with statement means we can use the parameter names as in parameters
  with(as.list(parameters),{
    dS <- - S/N * (beta_e * E + beta_s * Is + beta_a * Ia)
    dM <- - S/N * (beta_Me * E + beta_Ms * Is + beta_Ma * Ia)
    dE <- S/N * (beta_e * E + beta_s * Is + beta_a * Ia + beta_Me * E + beta_Ms * Is + beta_Ma * Ia) - delta * E
    dIs <- f * delta * E - gamma_s * Is
    dIa <- (1 - f) * delta * E - gamma_a * Ia
    dR <- gamma_s * Is + gamma_a * Ia
    res <- c(dS, dM, dE, dIs, dIa, dR)
    list(res)
  })
}

run_model <- function(S0, M0, E0, Is0, Ia0, R0, beta_e, beta_s, beta_a, beta_Me, beta_Ms, beta_Ma, delta, gamma_s, gamma_a, f, t_total){
  out <- lsoda(func = SMEIIR,
               y = c(S = S0, M = M0, E = E0, Is = Is0, Ia = Ia0, R = R0),
               parms = c(beta_e = beta_e, beta_s = beta_s, beta_a = beta_a, beta_Me = beta_Me, beta_Ms = beta_Ms, beta_Ma = beta_Ma, delta = delta, gamma_s = gamma_s, gamma_a = gamma_a, f = f),
               times = seq(1, t_total, by = 1))
  return(data.frame(out[,c('S','M','E','Is','Ia','R')]))
}

calc_distance <- function(D, D_star){
  return(sqrt(sum((D - D_star)**2)))
}