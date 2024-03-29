# First we extract the data and do a bit of data cleaning

data <- read.csv('../Data/england_2022-11-03.csv')
data$date <- as.Date(data$date)
data <- data[order(as.Date(data$date, format="%Y/%m/%d")),]
names(data)[names(data) == "newCasesByPublishDate"] <- "daily_I"
data <- data[c('areaName','date', 'daily_I')]
rownames(data) <- NULL     

# Now we need to use the data to create our S, I and R variables
# For the basic SIR model we assume that N=S+I+R is constant and that 
# every infected person recovers in 5 days
# https://www.imperial.ac.uk/news/239213/real-world-study-details-average-duration-infectiousness/

N <- 56000000 # https://datacommons.org/tools/timeline#&place=wikidataId/Q21&statsVar=Count_Person
data$cum_I <- cumsum(data$daily_I)
data$S <- N - data$cum_I
data$I <- 0
data$R <- 0
data$R[6:length(data$I)] <- data$daily_I[1:(length(data$I)-5)]
data$cum_R <- cumsum(data$R)
data$I <- data$cum_I - data$cum_R # same as using data$I <- N - data$S - data$cum_R

write.csv(data,"../Data/england_SIR_data.csv")