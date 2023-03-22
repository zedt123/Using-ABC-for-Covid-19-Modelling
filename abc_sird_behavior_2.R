source('sird_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_BSIRD_data.csv')

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','R','Dw','D')]

t_final <- length(data$date) #last day of the epidemic

#### ABC Setup ####
N <- 500 # Number of accepted particles
epsilon <- 585895 # Epsilon value
n_par <- 5 # How many parameters will be estimated
res <- matrix(ncol = n_par + 1, nrow = N) # empty matrix to store values
colnames(res) <- c('alpha_1','alpha_2','gamma','zeta_1','zeta_2','distance')

#### ABC Algorithm ####

i <- 1 # initialize counter of accepted particles
j <- 1 # initialize counter of  proposed particles

# Initialise using the prior

while (TRUE){
  alpha_1_0 <- runif(1, 0.4, 0.6)
  alpha_2_0 <- runif(1, 1000, 20000)
  gamma_0 <- runif(1, 0.4, 0.6)
  zeta_1_0 <- runif(1, 0, 0.01)
  zeta_2_0 <- runif(1, 0, 4)
  
  D_star <- run_model(data$S[1], data$I[1], data$R[1], data$Dw[1], data$D[1], alpha_1_0, alpha_2_0, gamma_0, zeta_1_0, zeta_2_0, t_final)
  distance <- calc_distance(data$I, D_star$I) + calc_distance(data$D,D_star$D)
  if (distance<epsilon){
    break
  }
}

#res[1,] <- c(0.5088866, 732.0650, 0.4682686, 1.639989e-04, 7.108274, 622417.6)
res[1,] <- c(0.4926692, 1040.646, 0.4484920, 0.004006870, 3.058122, 663598.7)
res[1,] <- c(0.4182294, 13327.57, 0.3608298, 0.003771990, 1.601466, 828283.8)
res[1,] <- c(0.4947707, 13327.26, 0.4461707, 0.0013732647, 1.281142, 669800.5)
res[1,] <- c(alpha_1_0, alpha_2_0, gamma_0, zeta_1_0, zeta_2_0, distance)

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
  
  if(zeta_2_star<0){
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
alpha_1 <- 0.5088866
alpha_2 <- 732.0650
gamma <- 0.4682686
zeta_1 <- 1.639989e-04
zeta_2 <- 7.108274

alpha_1 <- 0.4926692
alpha_2 <- 1040.646
gamma <- 0.448492
zeta_1 <- 0.00400687
zeta_2 <- 3.058122

alpha_1 <- 0.4182294
alpha_2 <- 13327.57
gamma <- 0.3608298
zeta_1 <- 0.00377199
zeta_2 <- 1.601466

alpha_1 <- 0.5093654
alpha_2 <- 1040.626
gamma <- 0.4688967
zeta_1 <- 4.359170e-05
zeta_2 <- 2.979770

D_star <- run_model(data$S[1], data$I[1], data$R[1], data$Dw[1], data$D[1], alpha_1, alpha_2, gamma, zeta_1, zeta_2, t_final)

calc_distance(data$I, D_star$I)

plot(seq(1,t_final,by=1), data$I, pch=16, ylim = c(0,300000), xlab = 'Day', ylab = 'Number of Cases')
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

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

beta <- 0.5086747
gamma <- 0.4683127
D_star_2 <- run_model(data$S[1], data$I[1], data$R[1], beta, gamma, t_final)


source('case_1/sir_sin_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 0.6336706
alpha_2 <- 0.03869089
alpha_3 <- 0.05035212  
gamma <- 0.5795639

D_star_4 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, alpha_3, gamma, t_final)


source('sir_behavior_preamble.R') # load the necessary data
data <- read.csv('../Data/england_SIR_data.csv')
data <- data[-1]

data <- data[data$date >= '2020-08-28' & data$date <= '2021-04-11',]
rownames(data) <- NULL
data <- data[c('date','S','I','cum_R')]
colnames(data)[4] <- 'R'

alpha_1 <- 0.6841056
alpha_2 <- -20.16637
gamma <- 0.6513138

D_star_3 <- run_model(data$S[1], data$I[1], data$R[1], alpha_1, alpha_2, gamma, t_final)



points <- data.frame(Day = seq(1,t_final,by=1),
                     data = data$I,
                     sir_death = D_star$I,
                     sir_behave = D_star_3$I,
                     sir_cos = D_star_4$I,
                     sir = D_star_2$I)

library(ggplot2)
ggplot(points, aes(x=Day)) + 
  geom_point(aes(y = data, color = "Data", shape = "Data"), size = 3) +  
  geom_point(aes(y = sir, color = "SIR Model", shape = "SIR Model"), size = 3.5) +
  geom_point(aes(y = sir_cos, color = "SIR Cosine  Model", shape = "SIR Cosine  Model"), size = 4) +
  geom_point(aes(y = sir_death, color = "SIRD Behavioral  Model", shape = "SIRD Behavioral  Model"), size = 2.5) +
  geom_point(aes(y = sir_behave, color = "SIR Behavioral  Model", shape = "SIR Behavioral  Model"), size = 3) +
  scale_color_manual(values=c("blue",'grey','purple',"green","red")) +
  scale_shape_manual(values=c(16, 19 ,18, 17, 15)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.background = element_blank()) +
  labs(color = " ", shape = " ") 

