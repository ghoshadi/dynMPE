# dynMPE
Dynamic marginal policy effects in thresholding designs, as proposed by Ghosh and Wager (2025).

The development version of this package can be installed using remotes:

```R
remotes::install_github("ghoshadi/dynMPE")
```

Example usage:

```R
library(dynMPE)
# A threshold rule applied to the same units period after period
set.seed(1)
n <- 500; periods <- 8
sim <- do.call(rbind, lapply(seq_len(n), function(i) {
  z <- numeric(periods); y <- numeric(periods); x <- 100 + 9 * rnorm(1)
  for (t in seq_len(periods)) {
    z[t] <- x
    x <- 100 + 0.9 * (x - 100) - (x >= 110) * 0.1 * (x - 100) + 4 * rnorm(1)
    y[t] <- -x
  }
  data.frame(id = i, period = 0:(periods - 1), z = z, y = y)
}))
fit <- dynMPE(sim, "y", "z", "period", "id", threshold = 110, gamma = 0.8)
print(fit)
summary(fit)
```

#### References
Aditya Ghosh and Stefan Wager.
<b>Non-parametric Causal Inference in Dynamic Thresholding Designs.</b> [arXiv preprint arXiv:2512.15244](https://arxiv.org/abs/2512.15244).

#### Funding
This research was supported by the Office of Naval Research under grant number N00014-24-1-2091.
