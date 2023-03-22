source('sird_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIRD_data.csv')

# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','R','D')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 500 # Number of accepted particles
epsilon <- 61803.80 # Epsilon value
n_par <- 3 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('beta','gamma','zeta','distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  beta_0 <- runif(1, 0, 10)
  gamma_0 <- runif(1, 0, 10)
  zeta_0 <- runif(1, 0, 10)
  
  D_star <- run_model(data$S[1], data$I[1], data$R[1], data$D[1], beta_0, gamma_0, zeta_0, t_final)
  distance <- calc_distance(data$I, D_star$I) + calc_distance(data$D,D_star$D)
  if (distance<epsilon){
    break
  }
}

res[1,] <- c(5.635351, 5.424854, 0.06097443, 61803.80)
#res[1,] <- c(5.78, 5.51, 0.1138467, 274026.6)
#res[1,] <- c(beta_0, gamma_0, zeta_0, distance)

while(i < N){
  # Generate a proposal
  beta_star <- rnorm(1,res[i,1],0.01)
  gamma_star <- rnorm(1,res[i,2],0.01)
  zeta_star <- rnorm(1,res[i,3],0.01)
  
  if(beta_star<0){
    next
  }
  
  if(gamma_star<0){
    next
  }
  
  if(zeta_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], data$D[1], beta_star, gamma_star, zeta_star, t_final)
 
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I) + calc_distance(data$D,D_star$D)
  
  # Accept or reject the parameters
  if(distance < epsilon){
    res[i+1,] <- c(beta_star, gamma_star, zeta_star, distance)
    epsilon <- distance
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)

beta <- 5.594663
gamma <- 5.383964
zeta <- 0.05875361

D_star <- run_model(data$S[1], data$I[1], data$R[1], data$D[1], beta, gamma, zeta, t_final)

calc_distance(data$I, D_star$I)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,30000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_sird_mcmc_results_4.csv')

# SIR, SEIR Models to compare results

source('case_1/case_1_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]
# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

beta <- 5.784274
gamma <- 5.622316
D_star_2 <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)

source('seir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIR_data.csv')

# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','E','I','R')]

beta_1 <- 4.744693
beta_2 <- 0.7367509
delta <- 5.313627
gamma <- 5.356273

D_star_3 <- run_model(data$S[1], data$E[1], data$I[1], data$R[1], beta_1, beta_2, delta, gamma, t_final)


points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     seiir = D_star$Is + D_star$Ia,
                     seir = D_star_3$I,
                     sir = D_star_2$I
)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = seiir, color = "SEIIR Model", shape = "SEIIR Model"), size = 4) +
  geom_point(aes(y = seir, color = "SEIR Model", shape = "SEIR Model"), size = 3) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  scale_color_manual(values=c("blue",'purple',"red","green")) +
  scale_shape_manual(values=c(16, 18, 15, 17)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

