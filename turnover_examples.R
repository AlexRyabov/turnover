#' Alexey Ryabov 2020
#' examples of using turnover, turnover_s and turnover_g

source("turnover.R")
##read data
##
data = read.csv( file = "Species.csv")
# Date,X,Y,Species1,Species2,...

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 4:N;


#richness based turnover index  
SERr = turnover_s(data[, SpecColumns]) #default parameters
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
#note that our richness based turnover index is equivalent to the binary distance in R
#so you can get the same result using 
SERr2 = (dist(data[, SpecColumns][, ], upper = TRUE, diag = TRUE, method = "binary"))
SERr2 = as.matrix(SERr2); #convert to a MxM matrix
upper_ind = mat_index(SERr2, "i<j"); #select upper triangular part
#both functions give the same result
plot(SERr$SER, SERr2[upper_ind])

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")

#plot SERr vs SERs
plot(SERr$SER, SERa$SER)

#richness based turnover characteristics (list of turnover index + other metrics)
lSERr = turnover(data[, SpecColumns], method = "SERr", ext_inv = TRUE) 

#abundance based turnover characteristics (list of turnover index + other metrics) 
lSERa = turnover(data[, SpecColumns], method = "SERa", ext_inv = TRUE)

#plot effective number of extinct species vs turnover index
plot(lSERa$SER, lSERa$S_ext)

#plot effective number of extinct species vs common number of species
plot(lSERa$S_common, lSERa$S_ext)


#include observation dates
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
lSERa = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates)
#plot turnover as a function time intervals
plot(lSERa$TimeIntv/365, lSERa$SER)

#plot turnover index as a function of euclidean distance between observations
XY = data[, 2:3]; 
lSERa = turnover(data[, SpecColumns], method = "SERa", locations =  XY)
plot(lSERa$Dist, lSERa$SER)

#plot turnover index as a function of euclidean distance between observations
Longitude = c(1:nrow(data))/nrow(data);
Latitude =  c(1:nrow(data))/nrow(data);
LonLat = data.frame(Longitude, Latitude);

lSERa = turnover(data[, SpecColumns], method = "SERa", locations =  LonLat, measure = "lonlat")
plot(lSERa$Dist, lSERa$SER)

#Group by some factors 
#define stations and areas 
nr = nrow(data)
Area = sample(c("Area 1", "Area 2"), size = nr, replace = TRUE)
Area = sort(Area)
Station = sample(c("A", "B", "C"), size = nr, replace = TRUE)
#
StatArea = data.frame(Area, Station);
SampleDates = as.Date(data$Date, format = "%Y-%m-%d"); #convert the input dates from string into class "date"
lSERa_SA = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea)
#lSERa_gr is a data frame with turnover within each group 
#show group names
summary(lSERa_SA)
unique(lSERa_SA$groupname)
plot(lSERa_SA$TimeIntv, lSERa_SA$SER)

#group by area only
lSERa_A = turnover(data[, SpecColumns], method = "SERa", dates =  SampleDates, groupby = StatArea[, "Area"])
summary(lSERa_A)
unique(lSERa_A$groupname)
plot(lSERa_A$TimeIntv, lSERa_A$SER)

