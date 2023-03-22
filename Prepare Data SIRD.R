# First we extract the data and do a bit of data cleaning

data <- read.csv('../Data/england_2022-11-03.csv')
data$date <- as.Date(data$date)
data <- data[order(as.Date(data$date, format="%Y/%m/%d")),]
names(data)[names(data) == "newCasesByPublishDate"] <- "daily_I"
data <- data[c('areaName','date', 'daily_I')]
rownames(data) <- NULL   

deaths <- read.csv('../Data/england_deaths.csv')
deaths$date <- as.Date(deaths$date)
deaths <- deaths[order(as.Date(deaths$date, format="%Y/%m/%d")),]
deaths <- deaths[deaths$date <= '2022-11-03',]

# Now we need to use the data to create our S, I and R variables
# For the basic SIR model we assume that N=S+I+R is constant and that 
# every infected person recovers in 5 days
# https://www.imperial.ac.uk/news/239213/real-world-study-details-average-duration-infectiousness/

N <- 56000000 # https://datacommons.org/tools/timeline#&place=wikidataId/Q21&statsVar=Count_Person
data$D <- c(rep(0,35),deaths$cumDeaths28DaysByPublishDate)
data$cum_I <- cumsum(data$daily_I)
data$S <- N - data$cum_I
data$I <- 0
data$R <- 0
data$R[6:length(data$I)] <- data$daily_I[1:(length(data$I)-5)]
data$R <- cumsum(data$R)
data$R <- data$R - data$D
data$I <- N - data$S - data$R - data$D

write.csv(data,"../Data/england_SIRD_data.csv")
