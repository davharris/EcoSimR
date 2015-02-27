

```r
library("EcoSimR")

set.seed(1)

## create “observed” data matrix ----------------------------------------------

m = matrix(rbinom(50^2, size = 1, prob = .1), ncol = 100)

burn_in = 1000    # number of iterations before first sample is collected
n_samples = 1000   # number of samples to collect
thin = 1000       # number of iterations between samples

n_species = nrow(m)
```

I think this is a correct way to get scaled pairwise C-scores...


```r
scaled_pairwise_c_score = function(x){
  shared = tcrossprod(m)
  sums = rowSums(m)
  upper = upper.tri(shared)
  raw_scores = (sums[row(shared)[upper]] - shared[upper]) *
    (sums[col(shared)[upper]] - shared[upper])
  
  raw_scores / tcrossprod(sums)[upper]
}
```

Record observed c_scores


```r
observed = scaled_pairwise_c_score(m)


## Null --------------------------------------------------------------------
```

Create empty matrix to store null results for all n-choose-2 pairs


```r
null = matrix(NA, nrow = choose(n_species, 2), ncol = n_samples)

message("burning in...")
```

```
## burning in...
```

```r
for(i in 1:burn_in){
  # Update matrix
  m = sim9_single(m)
}

message("sampling...")
```

```
## sampling...
```

```r
for(i in 1:n_samples){
  for(j in 1:thin){
    # Update matrix
    m = sim9_single(m)
  }
  # Fill in the next column of the null distribution
  null[ , i] = scaled_pairwise_c_score(m)  
}


## Calculate P-values ------------------------------------------------------
```

For each observed value, compare it to the null distribution from _ALL_ 
species pairs to calculate p-values


```r
p_lower = sapply(observed, function(x) mean(x <= null))
p_upper = sapply(observed, function(x) mean(x >= null))



## Evaluate P-value coverage -----------------------------------------------


mean(p_upper > .975) # proportion of p-values in upper tail
```

```
## [1] 0.37
```

```r
mean(p_lower < .025) # proportion of p-values in lower tail
```

```
## [1] 0
```

