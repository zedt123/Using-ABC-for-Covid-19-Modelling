source('sird_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_BSIRD_data.csv')

# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','I','R','Dw','D')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 500 # Number of accepted particles
epsilon <- 9745.028 # Epsilon value
n_par <- 5 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('alpha_1','alpha_2','gamma','zeta_1','zeta_2','distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  alpha_1_0 <- runif(1, 2, 7)
  alpha_2_0 <- runif(1, 900, 2000)
  gamma_0 <- runif(1, 2, 7)
  zeta_1_0 <- runif(1, 0, 5)
  zeta_2_0 <- runif(1, 0, 10)
  
  D_star <- run_model(data$S[1], data$I[1], data$R[1], data$Dw[1], data$D[1], alpha_1_0, alpha_2_0, gamma_0, zeta_1_0, zeta_2_0, t_final)
  distance <- calc_distance(data$I, D_star$I) + calc_distance(data$D,D_star$D)
  if (distance<epsilon){
    break
  }
}
# 27274.91
# 37966.26
res[1,] <- c(4.029458, 1822.153, 3.788962, 0.04084259, 0.7572414, 29650.30)
#res[1,] <- c(alpha_1_0, alpha_2_0, gamma_0, zeta_1_0, zeta_2_0, distance)

while(i < N){
  # Generate a proposal
  alpha_1_star <- rnorm(1,res[i,1],0.01)
  alpha_2_star <- rnorm(1,res[i,2],0.01)
  gamma_star <- rnorm(1,res[i,3],0.01)
  zeta_1_star <- rnorm(1,res[i,4],0.001)
  zeta_2_star <- rnorm(1,res[i,5],0.01)
  
  #alpha_1_star <- res[[i,1]]
  #alpha_2_star <- res[[i,2]]
  #gamma_star <- res[[i,3]]
  #zeta_1_star <- res[[i,4]]
  #zeta_2_star <- res[[i,5]]
  
  if(alpha_1_star<0){
    next
  }
  
  if(gamma_star<0){
    next
  }
  
  if(zeta_1_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star <- run_model(data$S[1], data$I[1], data$R[1], data$Dw[1], data$D[1], alpha_1_star, alpha_2_star,gamma_star, zeta_1_star, zeta_2_star, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$I, D_star$I)# + calc_distance(data$D,D_star$D)
  
  # Accept or reject the parameters
  if(distance < epsilon){
    res[i+1,] <- c(alpha_1_star, alpha_2_star,gamma_star, zeta_1_star, zeta_2_star, distance)
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
alpha_1 <- 4.029458
alpha_2 <- 1822.153
gamma <- 3.788962
zeta_1 <- 0.04084259
zeta_2 <- 0.7572414

alpha_1 <- 4.270075
alpha_2 <- 1822.395
gamma <- 3.970365
zeta_1 <- 0.09070381
zeta_2 <- 1.8499032

alpha_1 <- 4.260713
alpha_2 <- 967.5806
gamma <- 3.468241
zeta_1 <- 0.5795331
zeta_2 <- 6.022266

D_star <- run_model(data$S[1], data$I[1], data$R[1], data$Dw[1], data$D[1], alpha_1, alpha_2, gamma, zeta_1, zeta_2, t_final)

calc_distance(data$I, D_star$I)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,30000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$I, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

ts <- seq(1,t_final,by=1)
plot(ts,alpha_1 * exp(-alpha_2 * D_star$Dw/(D_star$S[1]+D_star$I[1]+D_star$R[1]+D_star$Dw[1]+D_star$D[1])),type = 'o', xlab = 'Day',ylab = 'Transmission Rate')


write.csv(res, 'England_seiir_mcmc_results_4.csv')

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

