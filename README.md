# Species turnover index
R functions for calculating the species turnover, based on Hillebrand et al. (2018). *Biodiversity change is uncoupled from species richness trends: Consequences for conservation and monitoring.* J Appl Ecol, 55, 169–184. https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12959

Input is an (M x N) table or matrix with M observations of N species and (optionally) observation dates and locations.

`turnover <- function (X, method = "SERr",  combinations = "i<j", 
                      dates = NULL, locations = NULL, measure = "euclidean", 
                      ext_inv = FALSE, groupby = NULL, RowID = NULL) `

 This function returns a data frame lSER of the turnover between rows i and j 
 of matrix X. lSER$From = i, lSER$To = j, lSER$SER is the turnover between these rows. 
 
The function result contains:
`lSER$From` is the row number of 'from' observation in X
`lSER$To` is the row number of 'to' observation
`lSER$SER` is the turnover between these observations
if `ext_inv = TRUE`, then the result includes additionally 
 the effective number of common/extinct/immigrating/total species for 
 observation i and j, `lSER$S_common`, `lSER$S_ext`, `lSER$S_imm`, `lSER$S_total` 
 e.g., `lSER$S_total` is the total number of species in observation i and j
 and  `lSER$S_ext` the number of species observed in i but not in j
 Note that the effective species number is richness based when  `method="SERr"` and 
 abundance (Simpson index) based when method="SERa", 
 
 if `dates` are specified, `lSER$TimeIntv` is the time interval between observation i and j, 
 if `locations` are specified, `lSER$Dist` is the spatial distances between observation i and j, 

 `lSER = turnover(X, method = "SERa", dates = NULL, locations = NULL, measure = "euclidean", groupby)`
 calculates turnover within groups defined by `groupby` dataframe. 
 `groupby` must contain the same number of rows as `X`. Function `turnover(X, ...)` splits `X` 
 into groups and returns a data frame with the same structure as returned by  function `turnover(X, ...)`, 
 plus an additional field "groupname"   

 
 @param X  is  (M x N) data frame or matrix with M observations of N species. 
   X can contain species abundances, frequencies or, only for SERr, presence/absence data 
   if X contains abundance then for calculating SERa the abundances will be normalized, 
   so that the sum in each row equal 1. 
 @param method  must be one of "SERr" and  "SERa". 
 If method = "SERr", then richness based turnover index is calculated 
 SER_ij  = (S_immigrant + S_extinct)/S_total, 
 and the effective numbers of species are richness based
 
 If method = "SERa", then  abundance based turnover index is  calculated
  SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2  sum_k( pik  pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )
 and the effective numbers of species are based on Simpson index
  
 @param combinations	indicates for which combinations of 
 rows i and j of matrix X the turnover must be calculated. The available 
 values are "i<j", "i!=j" and "i,j". 
 For combinations="i<j", the turnover is calculated between each row i and all subsequent 
 rows j of matrix X. This type is useful for calculating temporal turnover, 
 because in this case (if the data are sorted by date) we will calculate 
 turnover from past to future, but not vice versa.  
 For combinations="i!=j", turnover is calculated for all combinations of i and j, 
 except for i=j. 
 For combinations="i,j", the turnover is computed for all combinations 
 of i and j, in which case the turnover for i=j is zero. This option is useful 
 if the results should be converted later into a square matrix T, such that 
 element T[i,j] shows the turnover between rows X[i, :] and X[j, :].
  
 @param dates observation dates (datetime) or any numbers e.g. observation years, 
 
 @param locations spatial coordinates, typically should be Mx1, for 1D gradients, 
 Mx2 matrix if 2 coordinates are know, or Mx3 matrix for a 3D habitats
 or  longitude and latitude (in this order!) in this case measure = "lonlat"
 
 @param measure measure for measuring distance can be any of measures used 
 by dist(X, ..) function, or can be "lonlat" in this case the geodesic distance 
 in meters is measured  using 'geodist' package. 
 This package should be installed before. the 1st column of locations should 
 contain longitude and the second column latitude
 
 @param ext_inv should the information about invaded, extinct, common and total 
 species number be included?
  
 @param groupby  is a column or data frame with M rows  which are used to  
  group the data  (see split(x, ...) function). In this case turnover is 
  calculated within groups defined by "groupby" and resutl contains column 
  groupname.
  
  

 @examples
 `SER = turnover_s(X, method = "SERa")`
 
 `SER = turnover_s(X, method = "SERr")`
 
 `lSER = turnover_s(X, method = "SERa", dates)`
 spatial turnover, when coordinates are longitude and lattitude 
 `lSER = turnover_s(X, method = "SERa", location = LonLat, measure = "lonlat")`
 
 Calculating temporal species turnover separately for multiple 
 stations 
 `lSER_gr = turnover_g(X, method = "SERa", dates = ObservationDates, groupby = StationList)`
 where StationList is a data frame with station names for each observation. 
 Output: `lSER_gr$groupname` contains station names, 
 `lSER_gr$SER` is turnover information calculated within each station, but not between them. 
 




```R
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


```
