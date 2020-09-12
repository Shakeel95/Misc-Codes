# Key: simulation results

Three experiments are run, investigating the effects of:

* Changepoint location (`cpt_loc`)
* Changepoint density (`cpt_density`)
* Slope change magnitude (`slope_change_mag`)

Within each 4 dependence structures are considered:

1. Identity dispersion matrix with diagonals drawn uniformly from [1,3]
2. Toeplitz (autoregressive) spatial dependence with variances drawn uniformly from [1,3]
3. Same as 2 but with T.3 as opposed to Gaussian errors
4. Same as 2 but with autoregressive serial dependence of order 1 (decay = 0.5)
