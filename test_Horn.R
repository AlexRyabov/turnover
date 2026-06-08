#Testing Horn turnover (overlap) index in this package vs dissCqN::dissCqN
# test completed, the results match

library(dissCqN)
source("turnover.R")

dat <- read.csv("Species.csv", check.names = FALSE)
X <- as.matrix(dat[, 4:ncol(dat)])
mode(X) <- "numeric"

#Find Horn dissimilarity using dissCqN::dissCqN
system.time({
  D2 <- dissCqN::dissCqN(X, q = 1, pairwise = TRUE)
})

# dissCqN::dissCqN  system time
# User      System verstrichen 
# 17.54        0.53       18.31 
D2 <- if (is.list(D2)) D2[[1]] else D2
D2 <- as.matrix(D2)


#Find Horn dissimilarity using turnover_h
system.time({
  D1 <- turnover_h(X)$SER
})
# turnover_h  system time  
# User      System verstrichen 
# 0.14        0.03        0.18 

##
##  turnover_h is about 100 times faster than  dissCqN::dissCqN
##  for a matrix with 100 samples and 200 species 


#Find Horn dissimilarity using public turnover(..., method = "Horn")
system.time({
  D3_df <- turnover(X, method = "Horn", combinations = "i,j")
})

D3 <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(X))
D3[cbind(D3_df$From, D3_df$To)] <- D3_df$SER

cat("\nall.equal result for turnover_h vs dissCqN:\n")
print(all.equal(D1, D2, tolerance = 1e-10, check.attributes = FALSE))

cat("\nall.equal result for turnover(method = \"Horn\") vs dissCqN:\n")
print(all.equal(D3, D2, tolerance = 1e-10, check.attributes = FALSE))

# Compare lower triangles only
idx <- lower.tri(D1)

v_user <- D1[idx]
v_pkg  <- D2[idx]

keep <- is.finite(v_user) & is.finite(v_pkg)
v_user <- v_user[keep]
v_pkg  <- v_pkg[keep]

# Summary statistics
cor_val <- cor(v_user, v_pkg)
rmse_val <- sqrt(mean((v_user - v_pkg)^2))
max_abs_diff <- max(abs(v_user - v_pkg))

cat("\nComparison summary:\n")
cat("Number of compared pairs:", length(v_user), "\n")
cat("Correlation            :", cor_val, "\n")
cat("RMSE                   :", rmse_val, "\n")
cat("Max abs difference     :", max_abs_diff, "\n")

cat("\nSummary of differences (turnover_h - dissCqN):\n")
print(summary(v_user - v_pkg))

idx3 <- lower.tri(D3)
v_turnover <- D3[idx3]
v_pkg3 <- D2[idx3]
keep3 <- is.finite(v_turnover) & is.finite(v_pkg3)

cat("\nSummary of differences (turnover method Horn - dissCqN):\n")
print(summary(v_turnover[keep3] - v_pkg3[keep3]))

# Scatterplot
plot(
  v_pkg, v_user,
  xlab = "Horn dissimilarity from dissCqN",
  ylab = "Horn dissimilarity from turnover_h",
  main = "Comparison of Horn dissimilarity",
  pch = 16,
  cex = 0.6
)
abline(0, 1, lwd = 2, lty = 2)


#the results are identical 
