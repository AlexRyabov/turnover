#' Alexey Ryabov 2026
#' examples of using community turnover (overlap) indices with turnover() and turnover_s()

source("turnover.R")
set.seed(1)

##read data
##
data = read.csv( file = "Species.csv")
# Date,X,Y,Species1,Species2,...

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 4:N;


#richness based turnover values
SERr = turnover_s(data[, SpecColumns]) #default parameters
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
U0 = turnover_s(data[, SpecColumns], method = "U0") #alias for SERr
#note that our richness based turnover values are equivalent to the binary distance in R
#so you can get the same result using 
SERr2 = (dist(data[, SpecColumns][, ], upper = TRUE, diag = TRUE, method = "binary"))
SERr2 = as.matrix(SERr2); #convert to a MxM matrix
upper_ind = mat_index(SERr2, "i<j"); #select upper triangular part
#both functions give the same result
plot(SERr$SER, SERr2[upper_ind])

#Horn turnover (overlap) values
Horn = turnover_s(data[, SpecColumns], method = "Horn")
U1 = turnover_s(data[, SpecColumns], method = "U1") #alias for Horn

#abundance based turnover values
SERa = turnover_s(data[, SpecColumns], method = "SERa")

#Hill-number order q = 2 overlap
U2 = turnover_s(data[, SpecColumns], method = "U2")

#plot SERr vs SERs
plot(SERr$SER, SERa$SER)

#richness based turnover characteristics (list of turnover values + other metrics)
turnover_df_r = turnover(data[, SpecColumns], method = "SERr", ext_inv = TRUE) 

#abundance based turnover characteristics (list of turnover values + other metrics) 
turnover_df_a = turnover(data[, SpecColumns], method = "SERa", ext_inv = TRUE)

#plot effective number of extinct species vs turnover value
plot(turnover_df_a$SER, turnover_df_a$S_ext)

#plot effective number of extinct species vs common number of species
plot(turnover_df_a$S_common, turnover_df_a$S_ext)


#include observation dates
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
turnover_df = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates)
head(turnover_df[, c("From", "To", "dateFrom", "dateTo", "TimeIntv", "SER")])
#plot turnover as a function time intervals
plot(turnover_df$TimeIntv/365, turnover_df$SER)

#if you use years then dateFrom and dateTo will contain the year
Year <- as.numeric(format(SampleDates, "%Y")); #convert the input dates from string into class "date" 
turnover_df = turnover(data[, SpecColumns], method = "SERa", dates =  Year)


#plot turnover value as a function of euclidean distance between observations
XY = data[, 2:3]; 
turnover_df = turnover(data[, SpecColumns], method = "SERa", locations =  XY)
plot(turnover_df$Dist, turnover_df$SER)

#plot turnover value as a function of geodesic distance between observations
Longitude = c(1:nrow(data))/nrow(data);
Latitude =  c(1:nrow(data))/nrow(data);
LonLat = data.frame(Longitude, Latitude);

if (requireNamespace("geodist", quietly = TRUE)) {
  turnover_df = turnover(data[, SpecColumns], method = "SERa", locations =  LonLat, measure = "lonlat")
  plot(turnover_df$Dist, turnover_df$SER)
} else {
  message("Skipping measure = 'lonlat' example because package 'geodist' is not installed.")
}

#Group by some factors 
#define stations and areas 
nr = nrow(data)
Area = sample(c("Area 1", "Area 2"), size = nr, replace = TRUE)
Area = sort(Area)
Station = sample(c("A", "B", "C"), size = nr, replace = TRUE)
#
StatArea = data.frame(Area, Station);
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
turnover_df_SA = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea)
#turnover_df_SA is a data frame with turnover within each group 
#show group names
summary(turnover_df_SA)
unique(turnover_df_SA$groupname)
plot(turnover_df_SA$TimeIntv, turnover_df_SA$SER)

#group by area only
turnover_df_A = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea[, "Area"])
summary(turnover_df_A)
unique(turnover_df_A$groupname)
plot(turnover_df_A$TimeIntv, turnover_df_A$SER)

