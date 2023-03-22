source('smeiir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIIR_data.csv')

# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','E','Ia','Is','R')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 100 # Number of accepted particles
epsilon <- 25714.51 # Epsilon value
n_par <- 9 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('beta_e','beta_s','beta_a','beta_Me','beta_Ms','beta_Ma','delta','gamma_s','gamma_a','distance')
f <- 0.4 # https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2787098

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  beta_e_0 <- runif(1, 0, 10)
  beta_s_0 <- runif(1, 0, 10)
  beta_a_0 <- runif(1, 0, 10)
  delta_0 <- runif(1, 0, 10)
  gamma_s_0 <- runif(1, 0, 10)
  gamma_a_0 <- runif(1, 0, 10)
  
  beta_Me_0 <- beta_e_0/2
  beta_Ms_0 <- beta_s_0/2
  beta_Ma_0 <- beta_a_0/2
  
  D_star <- run_model(0.85*data$S[1], 0.15*data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e_0, beta_s_0, beta_a_0, beta_Me_0, beta_Ms_0, beta_Ma_0, delta_0, gamma_s_0, gamma_a_0, f, t_final)
  distance <- calc_distance(data$Is, D_star$Is) + calc_distance(data$Ia,D_star$Ia) # I is the sum of Is and Ia
  if (distance<epsilon){
    break
  }
}


res[1,] <- c(5.750479, 0.008770471, 3.172277, 2.875239, 0.004385235, 1.586138, 11.64997, 2.493874, 6.268006, 25714.51)
#res[1,] <- c(beta_e_0, beta_s_0, beta_a_0, beta_Me_0, beta_Ms_0, beta_Ma_0, delta_0, gamma_s_0, gamma_a_0, distance)

while(i < N){
  # Generate a proposal
  beta_e_star <- rnorm(1,res[i,1],0.1)
  beta_s_star <- rnorm(1,res[i,2],0.1)
  beta_a_star <- rnorm(1,res[i,3],0.1)
  delta_star <- rnorm(1,res[i,7],0.1)
  gamma_s_star <- rnorm(1,res[i,8],0.1)
  gamma_a_star <- rnorm(1,res[i,9],0.1)
  
  beta_Me_star <- beta_e_star/2
  beta_Ms_star <- beta_s_star/2
  beta_Ma_star <- beta_a_star/2
  
  if(beta_e_star<0){
    next
  }
  
  if(beta_s_star<0){
    next
  }
  
  if(beta_a_star<0){
    next
  }
  
  # Simulate dataset from the model
  D_star <- run_model(0.85*data$S[1], 0.15*data$S[1], data$E[1], data$Is[1], data$Ia[1],  data$R[1], beta_e_star, beta_s_star, beta_a_star, beta_Me_star, beta_Ms_star, beta_Ma_star, delta_star, gamma_s_star, gamma_a_star, f, t_final)
  
  # Calculate distance
  distance <- calc_distance(data$Is, D_star$Is) + calc_distance(data$Ia,D_star$Ia)
  
  # Accept or reject the parameters
  if(distance <= epsilon){
    res[i+1,] <- c(beta_e_star, beta_s_star, beta_a_star, beta_Me_star, beta_Ms_star, beta_Ma_star, delta_star, gamma_s_star, gamma_a_star, distance)
    i <- i + 1
  }  
  j <- j + 1
  acc_rate <- i / j
  if(j %% 10000 == 0){
    cat("Iterations:" , j , " " , 'Accepted:' , i, "\n")
  }
}

View(res)

beta_e <- 6.080017
beta_s <- 0.1992957997
beta_a <- 2.945709
delta <- 12.55609
gamma_s <- 2.646483
gamma_a <- 6.209167

beta_Me <- beta_e/2
beta_Ms <- beta_s/2
beta_Ma <- beta_a/2

D_star <- run_model(0.85*data$S[1], 0.15*data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e, beta_s, beta_a, beta_Me, beta_Ms, beta_Ma, delta, gamma_s, gamma_a, f, t_final)

calc_distance(data$Is+data$Ia, D_star$Is+D_star$Ia)

plot(seq(1,t_final,by=1), data$Is+data$Ia, pch=16, ylim = c(0,30000), xlab = 'Day', ylab = 'Number of Cases')
points(seq(1,t_final,by=1), D_star$Is+D_star$Ia, pch=15, col = 'blue')
legend("topright", legend=c("Real Data", "Simulated Data"),
       col=c("black", "blue"), pch=c(16,15), cex=0.8)

write.csv(res, 'England_smeiir_mcmc_results.csv')

# SIR, SEIR, SEIIR Models to compare results

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


source('seiir_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SEIIR_data.csv')

# We will try to estimate the first wave of COVID in England: 2020-03-01 - 2020-06-30
data <- data[data$date >= '2020-03-06' & data$date <= '2020-06-30',]
rownames(data) <- NULL
data <- data[c('date','S','E','Ia','Is','R')]

beta_e <- 7.203261
beta_s <- 0.4934008
beta_a <- 4.523558
delta <- 15.12028
gamma_s <- 2.457595
gamma_a <- 5.832009

D_star_4 <- run_model(data$S[1], data$E[1], data$Is[1], data$Ia[1], data$R[1], beta_e, beta_s, beta_a, delta, gamma_s, gamma_a, f, t_final)

points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$Is + data$Ia,
                     smeiir = D_star$Is + D_star$Ia,
                     seiir = D_star_4$Is + D_star_4$Ia,
                     seir = D_star_3$I,
                     sir = D_star_2$I)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +
  geom_point(aes(y = seiir, color = "SEIIR Model", shape = "SEIIR Model"), size = 4) +
  geom_point(aes(y = seir, color = "SEIR Model", shape = "SEIR Model"), size = 3) +
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3) +
  geom_point(aes(y = smeiir, color = "SMEIIR Model", shape = "SMEIIR Model"), size = 3) +
  scale_color_manual(values=c("blue",'purple',"red","green", "orange")) +
  scale_shape_manual(values=c(16, 18, 15, 17, 19)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

