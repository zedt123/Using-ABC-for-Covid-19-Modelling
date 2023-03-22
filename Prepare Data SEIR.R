# First we extract the data and do a bit of data cleaning

data <- read.csv('../Data/england_2022-11-03.csv')
data$date <- as.Date(data$date)
data <- data[order(as.Date(data$date, format="%Y/%m/%d")),]
names(data)[names(data) == "newCasesByPublishDate"] <- "daily_I"
data <- data[c('date', 'daily_I')]
rownames(data) <- NULL     

# Now we need to use the data to create our S, E, I and R variables
# For the SEIR model we assume that N=S+E+I+R is constant and that 
# every infected person recovers in 5 days
# https://www.imperial.ac.uk/news/239213/real-world-study-details-average-duration-infectiousness/

# We also assume that the exposed period lasts 5 days
# https://www.health.harvard.edu/diseases-and-conditions/if-youve-been-exposed-to-the-coronavirus#:~:text=How%20soon%20after%20I%27m,days%20for%20the%20Delta%20variant

N <- 56000000 # https://datacommons.org/tools/timeline#&place=wikidataId/Q21&statsVar=Count_Person
data$daily_E <- c(data$daily_I[6:length(data$date)],0,0,0,0,0)
data$cum_I <- cumsum(data$daily_I)
data$cum_E <- cumsum(data$daily_E) + 2
data$S <- N - data$cum_E
data$I <- 0
data$daily_R <- 0
data$daily_R[6:length(data$I)] <- data$daily_I[1:(length(data$I)-5)]
data$R <- cumsum(data$daily_R)
data$E <- c(0,1,1,1,zoo::rollsum(data$daily_E,5,))
data$I <- N - data$S - data$E - data$R

write.csv(data,"../Data/england_SEIR_data.csv")
