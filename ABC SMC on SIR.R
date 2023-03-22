source('case_1/case_1_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N_t <- 5 # Number of intermediate distributions
N <- 100 # Number of accepted particles
epsilon <- c(50000,45000,40000,39000,38000) # Epsilon values
n_par <- 2 # How many parameters will be estimated
res <- matrix(ncol = n_par + 2, nrow = N) # empty matrix to store values
weights <- rep(0,N)
colnames(res) <- c('i','beta', 'gamma', 'distance')

#### ABC Algorithm ####
t <- 1 # initialize t
i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# For t=1
while(i <= N){
  # Sample from prior
  beta_star <- runif(1, 4.5, 7)
  gamma_star <- runif(1, 4.5, 7)
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], beta_star, gamma_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance <= epsilon[t]){
    weights[i] <- 1
    res[i,] <- c(i,beta_star, gamma_star, distance)
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}
print(j)
weights <- weights / N # Normalise the weights
t <- 2
new_res <- matrix(ncol = n_par + 2, nrow = N) # keep in memory the current samples
colnames(new_res) <- c('i','beta', 'gamma', 'distance')
new_weights <- rep(0,N) # keep in memory the current weights

sigma_1 <- 0.1 # standard deviations for the matrices
sigma_2 <- 0.1

# for t>1
while(t<=N_t){
  i <- 1
  j <- 1
  while(i<=N){
    star <- sample(res[,1],1,prob=weights) # choose a pair according to the weights
    beta_star <- res[[star,2]]
    gamma_star <- res[[star,3]]
    
    k_1 <- runif(1,-1,1) # perturbation matrices
    k_2 <- runif(1,-1,1)
    
    beta_star <- beta_star + sigma_1 * k_1
    gamma_star <- gamma_star + sigma_2 * k_2
    
    if(beta_star < 4.5 | beta_star > 7.0){ # if after perturbation it lies outside the support go to the next iteration
      next
    }
    
    if(gamma_star < 4.5 | gamma_star > 7.0){
      next
    }     
    
    # Simulate dataset from the model
    D_star <- run_model(data$S[1], data$I[1], data$R[1], beta_star, gamma_star, t_final)
    
    # Calculate distance
    distance <- calc_distance(data$I, D_star$I)  
    # Accept or reject the parameters
    if(distance <= epsilon[t]){
      weight_denominator <- 0
      for(k in 1:N){
        if(beta_star<res[[k,2]]-sigma_1 | beta_star>res[[k,2]]+sigma_1){ # beta
          next
        }
        if(gamma_star<res[[k,3]]-sigma_2 | gamma_star>res[[k,3]]+sigma_2){ # gamma
          next
        }
        weight_denominator <- weight_denominator + weights[k]
      }
      new_weights[i] <- 0.4/weight_denominator
      new_res[i,] <- c(i,beta_star, gamma_star, distance)
      i <- i + 1
    }  
    j <- j + 1
    acc_rate <- i / j
    if(j %% 10000 == 0){
      cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
    }
  }
  weights <- new_weights # keep in memory the weights because we will need them for the next iteration
  weights <- weights / sum(weights) # normalise the weights
  res <- new_res # keep in memory the samples because we will need them for the next iteration
  t <- t+1 # move to the next intermediate distribution
  print(j)
}

View(res)
beta <- 5.784274
gamma <- 5.622316
D_star <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,30000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_basic_smc_results.csv')

res <- read.csv('England_basic_smc_results.csv')

# Graph with all the curves
library(ggplot2)

# Rejection ABC
beta_1 <- 5.775527
gamma_1 <- 5.614064
D_1 <- run_model(data$S[1], data$I[1], data$R[1], beta_1, gamma_1, t_final)

# ABC MCMC
beta_2 <- 5.787519
gamma_2 <- 5.624737
D_2 <- run_model(data$S[1], data$I[1], data$R[1], beta_2, gamma_2, t_final)

# ABC SMC
beta_3 <- 5.784274
gamma_3 <- 5.622316
D_3 <- run_model(data$S[1], data$I[1], data$R[1], beta_3, gamma_3, t_final)

points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     rej = D_1$I,
                     mcmc = D_2$I,
                     smc = D_3$I)

ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = rej, color = "Rejection", shape = "Rejection"), size = 3) +
  geom_point(aes(y = mcmc, color = "MCMC", shape = "MCMC"), size = 3) +
  geom_point(aes(y = smc, color = "SMC", shape = "SMC"), size = 3) +
  scale_color_manual(values=c("blue","red","green","purple")) +
  scale_shape_manual(values=c(16, 17, 15, 18)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 


(37404+37411+37420+37429+37485+37511+37516+37546+37582+37597)/10

(37393+37394+37414+37432+37445+37507+37542+37546+37558+37574)/10

(37393+37395+37396+37401+37420+37420+37438+37452+37455+37456)/10

