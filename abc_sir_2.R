source('case_1/case_1_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 100 # Number of accepted particles
epsilon <- 573911.7 # Epsilon value
n_par <- 2 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('beta', 'gamma', 'distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  beta_0 <- runif(1, 0, 1)
  gamma_0 <- runif(1, 0, 1)
  D_star <- run_model(data$S[1], data$I[1], data$R[1], beta_0, gamma_0, t_final)
  distance <- calc_distance(data$I, D_star$I)
  if (distance<epsilon){
    break
  }
}

res[1,] <- c(0.5080517, 0.4675646, 573911.7)
res[1,] <- c(beta_0, gamma_0, distance)

while(i < N){
  # Generate a proposal
  beta_star <- rnorm(1,res[i,1],0.01)
  gamma_star <- rnorm(1,res[i,2],0.01)
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], beta_star, gamma_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)
  
  # Accept or reject the parameters
  if(distance <= epsilon){
    res[i+1,] <- c(beta_star, gamma_star, distance)
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)
beta <- 0.5086747
gamma <- 0.4683127
D_star <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_basic_mcmc_results_analysis_2.csv')