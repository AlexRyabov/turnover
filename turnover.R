#' Alexey Ryabov 2020
#' Calculate species turnover index
#' see details in 
#' Hillebrand, H. et al. J Appl Ecol 55, 169-184 (2018).

#' 
#' @description
#' 
#' SER = turnover_s(X, method = "SERa")
#' returns a matrix SER(i, j) with  turnover indexes between the rows i and j of the data in X
#' 
#' 
#' lSER = turnover(X, method = "SERa", dates = NULL, locations = NULL)
#' Return a list lSER of M x M matrices characterizing species turnover between observation i and j of X 
#' measured at given dates and locations. 
#' 
#' The function result contains:
#' lSER$SER is an upper triangle M x M matrix with turnover index. lSER$SER[i,j] turnover between row i and j in X
#' 
#' lSER$S_common, lSER$S_ext, lSER$S_imm, lSER$S_total are an upper triangle MxM matrices  
#' the effective number of common/extinct/immigrating/total species,
#' e.g., lSER$S_total[i,j] is the total number of species in observation i and j
#' and  lSER$S_ext[i,j] the number of species observed in i but not in j
#' Note that the effective species number is richness based when  method="SERr" and 
#' abundance (Simpson index) based when method="SERa", 
#' see details in Hillebrand, H. et al. J Appl Ecol 55, 169-184 (2018).
#' 
#' lSER$TimeIntv[i,j] time intervals between observation i and j
#' 
#' lSER$Dist[i,j]  spatial distances between observation i and j
#' 
#' 
#' lSER_gr = turnover_g(X, method = "SERa", dates = NULL, locations = NULL, groupby)
#' This function calculates turnover within groups defined by groupby dataframe. 
#' The groupby must contain the same number of rows as X. The function splits X 
#' into groups and returns a sheet with two fields: lSER_gr$groupnames contains 
#' unique group names, lSER_gr$turnover contains information on species turnover 
#' within each group. 
#' 
#' 
#' 
#' 
#' 

#' @param X  is  (M x N) table or matrix with M observations of N species. 
#'   X can contain species abundances, frequencies or, only for SERr, presence/absence data 
#'   if X contains abundance then for calculating SERa the abundances will be normalized, 
#'   so that the sum in each row equal 1. 
#' @param method  equals  "SERr" for calculating richness based turnover index 
#' SER_ij  = (S_immigrant + S_extinct)/S_total, effective number of species is richness based
#' 
#' @param method  equals  "SERa" for calculating abundance based index 
#'  SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2  sum_k( pik  pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )
#' 
#'  
#' @param dates observation dates (datetime) or any numbers e.g. observation years, 
#' 
#' @param locations spatial coordinates, typically should be Mx1, for 1D gradients, 
#' Mx2 matrix if 2 coordinates are know, or Mx3 matrix for a 3D habitats
#'  
#' @param groupby (M x ..) a 'factor' in the sense that as.factor(f) defines 
#' the grouping, or a list of such factors in which case their interaction 
#' is used for the grouping.  (see split(x, ...) function)
#'  
#'  

#' @examples
#' SER = turnover_s(X, method = "SERa")
#' 
#' SER = turnover_s(X, method = "SERr")
#' 
#' lSER = turnover_s(X, method = "SERa", dates)
#' lSER = turnover_s(X, method = "SERa", NULL, location)
#' 
#' Calculating temporal species turnover separately for multiple 
#' stations 
#' lSER_gr = turnover_g(X, method = "SERa", ObservationDates, locations = NULL, StationList)
#' where StationList is a data frame with station names for each observation. 
#' Output: lSER_gr$groupnames contains a list of unique station names, 
#' lSER_gr$turnover is a list of turnover information for those stations. 
#' The elements of this list have a structure returned by the function 
#' lSER = turnover(X, ...)
#' For example, the name of the second group 
#' lSERa_gr$groupnames[2] 
#' Plot temporal turnover for the second group 
#' plot(lSERa_gr$turnover[[2]]$TimeIntv, lSERa_gr$turnover[[2]]$SER)
#' 


turnover_s <- function (X, method = "SERr") {
  #short version, returns only turnover index
  #between rows of the X 
  
  #X is  (M x N) table or matrix with M observations of abundances of N species, 
  #M observation -- rows
  #N species -- columns

  #method = 'SERa' output is based on relative species frequencies
  # SER_ij = (sum_k pik^2 + sum_k pjk^2 - 2 * sum_k( pik * pjk) ) / (sum_k pik^2 + sum_k pjk^2 -  sum_k( pik * pjk) )

  #method = 'SERr' output is based on presence/absence data
  #SER_ij  = (S_immigrant + S_extinct)/S_total

  #the output SER(,) is  MxM matrix with turnover rates between different pairs of observations. 
  #SER(i,j) is turnover from observation i to j

  Res = turnover(X, method, NULL)
  return(Res$SER)
}


turnover_g <- function(X, method = "SERr", dates = NULL, locations = NULL, groupby = GroupBy){
  Xsplit = split(X, groupby);
  if (!is.null(dates))     { DatesSplit = split(dates,     groupby)}
  DatesPiece = NULL;
  if (!is.null(locations)) { LocSplit   = split(locations, groupby)}
  LocPiece = NULL;  
  Res <- vector(mode = "list", length = length(Xsplit))
  for(i in 1:length(Xsplit)){
    if (!is.null(dates))      { DatesPiece = DatesSplit[[i]]}
    if (!is.null(locations))  { LocPiece =   LocSplit[[i]]}
    
     Res[[i]] =  turnover(Xsplit[[i]], method, DatesPiece, LocPiece)
  }
  return(list(groupnames = names(Xsplit), 'turnover' = Res))
}


turnover <- function (X, method = "SERr", dates = NULL, locations = NULL) {
  switch(method, 
         SERa={
           # abundance based turnover...
           Res = turnover_a(X, method);
         },
         SERr={
           # richness based turnover...
           Res = turnover_r(X, method);
         },
         {
           stop('turnover method is undefined')
         }
  )  

  if (!is.null(dates)){
    #get number of records 
    M <- nrow(X);
    
    M2 = length(dates)
    if (M != M2){
      stop("The length of the dates array must be equal to the number of rows in X")
    }
    
    #find differences in dates (if not null) or interval between columns
    D = rep(dates, M);
    D1 = matrix(D, M, M, byrow = FALSE);
    D2 = matrix(D, M, M, byrow = TRUE);
    
    #D1 = as.Date(D1, origin ="1970-01-01")  unwrap D1 and D2
    #D2 = as.Date(D2, origin ="1970-01-01")
    #TimeIntv = difftime(D2, D1, units);
    #need to pack to back into a matrix
    #TimeIntv = matrix(TimeIntv, M, M, byrow = FALSE)
    
    TimeIntv = D2 - D1; #interval in days or integer numbers
    #set elements below the main diagonal to zero
    for(i in 2:M){
      for (j  in 1:(i-1)){
        TimeIntv[i,j] = NA;
      }
    }
    Res$TimeIntv=TimeIntv;
  }
  
  if (!is.null(locations)){
    M <- nrow(X);
    M2 = nrow(locations)
    if (M != M2){
      stop("The number of rows in the location matrix must be equal to the number of rows in X.")
    }    
    Dist = as.matrix(dist(locations, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
    for(i in 2:M){
      for (j  in 1:(i-1)){
        Dist[i,j] = NA;
      }
    }
    Res$Dist=Dist;
  }  
  return(Res);
}


turnover_a <- function(X, method) {
  X = as.matrix(X)
  
  dims = dim(X)
  M = dims[1];
  N = dims[2];

    
  #to calculate abundance based turnover index, the data should be normalized
  #each row should contain species frequencies, for row j, \sum_i(p_ij)=1
  s = rowSums(X);
  Sums = matrix(rep(s, N), M, N);
  X = X/Sums;
  
  
  
  SimpsInd = as.matrix(rowSums(X^2));  #sum of squares, Simpson index
  P_sq =rep(SimpsInd, M);
  P1_sq =matrix(P_sq, M, M, byrow = FALSE);
  S1 = 1/P1_sq; #in each row the number of species in row i
  P2_sq =matrix(P_sq, M, M, byrow = TRUE);
  S2 = 1/P2_sq; #in each column the number of species in row i
  
  P12  = X %*% t(X);  #overlap 

  #find invaded, common and etc species 
  S_common  = (P12)/(P1_sq * P2_sq) ;
  S_ext   = S1 - S_common;
  S_imm   = S2 - S_common;
  S_total = S1 + S2 - S_common;

  SER = (S_imm + S_ext)/S_total;
  #SER2 = (P1_sq + P2_sq - 2 * P12)/(P1_sq + P2_sq - P12);
  
  Result = list(SER=SER, S_total = S_total, S_common = S_common, S_imm = S_imm, S_ext = S_ext);
  Result = set2nullLowTriang(Result);
  return(Result)
}



turnover_r <- function(X, method) {
  dims = dim(X)
  M = dims[1];
  N = dims[2];
  
  P1 = as.matrix(X > 0)*1;
  P2 = t(P1);

  S_common = P1 %*% P2;

  S_sample = as.matrix(rowSums(X>0));  #number of species in each sample
  S_sample =rep(S_sample, M);
  S1 =matrix(S_sample, M, M, byrow = FALSE); #put number of species columnwise  
  S2 =matrix(S_sample, M, M, byrow = TRUE);  #put number of species rowwise 
  
  S_total = S1 + S2 - S_common;
  S_imm = ((!P1) + 0) %*% (P2);
  S_ext =  P1 %*% ((!P2) + 0);
  SER = (S_imm + S_ext)/(S_total);
  Result = list(SER=SER, S_total = S_total, S_common = S_common, S_imm = S_imm, S_ext = S_ext);
  Result = set2nullLowTriang(Result);
  return(Result)
}

set2nullLowTriang <- function(Res){
  #Set to NULL all elements below the main diagonal, because this elements correspond  
  #to 'back-in-time' turnover rates
  M <- nrow(Res$SER);
  for(i in 2:M){
    for (j  in 1:(i-1)){
      Res$SER[i,j]      = NA;
      Res$S_total[i,j]  = NA;
      Res$S_common[i,j] = NA;
      Res$S_imm[i,j]    = NA;
      Res$S_ext[i,j]    = NA;
    }
  }
  return (Res);
}