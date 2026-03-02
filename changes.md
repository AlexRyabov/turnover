Changes 02.03.2026 

* **`ext_inv = TRUE`: added `S_imm` to the output data frame**

* **`groupby`: fixed grouping for matrices**
  `split(X, groupby)` behaves incorrectly for a `matrix` (it splits it like a vector). I replaced it with `split(as.data.frame(X), groupby)` to ensure splitting is done **row-wise** for both `data.frame` and `matrix` inputs.

* **`measure = "lonlat"`: fixed the `geodist()` call**
  The original code used `paired = TRUE`, which is meant for *paired distances between x and y* (vectors), while here we need a **full pairwise distance matrix**. According to the documentation, `paired = FALSE, sequential = FALSE` returns a full matrix.
  I also removed `library()` from inside the function and added a `requireNamespace()` check.

* **`SERa`: handled rows with zero total abundance**
  If some rows have total abundance = 0, normalization produces `NaN/Inf`. I added a warning and set those rows to `NA` so that pairs involving them become `NA`, while all other pairs are computed normally.


