# Community Turnover (Overlap) Indices

R functions for calculating pairwise community turnover (overlap) indices.

The main implementation and detailed API documentation are in [turnover.R](turnover.R). Complete runnable examples using `Species.csv` are in [turnover_examples.R](turnover_examples.R).

References:

- Hillebrand et al. (2018), *Biodiversity change is uncoupled from species richness trends: Consequences for conservation and monitoring*, J Appl Ecol, 55, 169-184.
- Gotelli, Chao & Colwell (2024), *Measuring and Estimating Species Richness, Taxonomic and Phylogenetic Diversity, and Related Biotic (Dis)similarity Indices From Sampling Data*, 314-339.

## Input

Input is an `M x N` table or matrix `X`, where rows are community observations and columns are species abundances, frequencies, or presence/absence values.

Optional vectors or data frames can provide observation dates, spatial coordinates, and grouping variables.

## Main Functions

```r
turnover(
  X,
  method = "SERr",
  combinations = "i<j",
  dates = NULL,
  locations = NULL,
  measure = "euclidean",
  ext_inv = FALSE,
  groupby = NULL
)
```

`turnover()` returns a data frame, usually named `turnover_df` in examples, with pairwise community turnover values and optional metadata such as dates, time intervals, spatial distances, and group names.

```r
turnover_s(X, method = "SERr", combinations = "i<j")
```

`turnover_s()` is a convenience wrapper for species-only calculations. It is equivalent to calling `turnover()` without `dates`, `locations`, `ext_inv`, or `groupby`.

## Methods

| Method | Alias | Meaning |
| --- | --- | --- |
| `method = "SERr"` | `method = "U0"` | Richness-based turnover for presence/absence data; equivalent to Jaccard dissimilarity / `dist(X, method = "binary")`. |
| `method = "Horn"` | `method = "U1"` | Shannon-entropy-based Horn turnover value for Hill number `q = 1`. |
| `method = "SERa"` |  | Abundance-based turnover from Hillebrand et al. (2018), the complement of the Wishart similarity ratio. |
| `method = "U2"` |  | Hill-number order `q = 2` turnover (overlap) value from the Chao `U_q` family. Identical communities have `U2 = 1`. |

See the roxygen comments in [turnover.R](turnover.R) for formulas, returned columns, and method-specific notes.

## Minimal Example

```r
source("turnover.R")

data <- read.csv("Species.csv")
SpecColumns <- 4:ncol(data)
X <- data[, SpecColumns]

SERr <- turnover(X, method = "SERr")
U0 <- turnover(X, method = "U0")
Horn <- turnover(X, method = "Horn")
U1 <- turnover(X, method = "U1")
SERa <- turnover(X, method = "SERa")
U2 <- turnover(X, method = "U2")

SampleDates <- as.Date(data$Date, format = "%Y-%m-%d")
XY <- data[, 2:3]

turnover_df <- turnover(X, method = "SERa", dates = SampleDates, locations = XY)
head(turnover_df)
```

## Full Examples

Complete examples are provided in [turnover_examples.R](turnover_examples.R) and can be run with:

```r
source("turnover_examples.R")
```
