# Minimal comparison of U2 non-overlap from turnover() vs dissCqN::dissCqN

library(dissCqN)
source("turnover.R")

dat <- read.csv("Species.csv", check.names = FALSE)
X <- as.matrix(dat[, 4:ncol(dat)])
mode(X) <- "numeric"

# C2 dissimilarity = 1 - C2 = Morisita-Horn dissimilarity
D_C2 <- dissCqN::dissCqN(X, q = 2, pairwise = TRUE)
D_C2 <- if (is.list(D_C2)) D_C2[[1]] else D_C2
D_C2 <- as.matrix(D_C2)

# Convert to U2 non-overlap.
D_U2_dissCqN <- D_C2 / (2 - D_C2)

# turnover(method = "U2") returns U2 overlap, so convert it to non-overlap.
U2_df <- turnover(X, method = "U2", combinations = "i,j")
U2_overlap <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(X))
U2_overlap[cbind(U2_df$From, U2_df$To)] <- U2_df$SER
D_U2_turnover <- 1 - U2_overlap

idx <- lower.tri(D_U2_dissCqN)
x <- D_U2_dissCqN[idx]
y <- D_U2_turnover[idx]
keep <- is.finite(x) & is.finite(y)
x <- x[keep]
y <- y[keep]

msd <- mean((y - x)^2)
cat("Mean square deviation:", msd, "\n")

plot(
  x, y,
  xlab = "U2 non-overlap from dissCqN",
  ylab = "U2 non-overlap from turnover(method = \"U2\")",
  main = paste("U2 non-overlap comparison, MSD =", signif(msd, 4)),
  pch = 16,
  cex = 0.6
)
abline(0, 1, lwd = 2, lty = 2)
