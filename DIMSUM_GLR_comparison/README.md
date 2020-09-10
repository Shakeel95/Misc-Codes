# DIM-SUM GLR comparison

For the canonical multivariate changepoint problem with a single change it been show ([Aston and Kirch (2018)][1]) that projecting onto the vector of differences between pre and post change means will result in no loss of information about the changepoint location. A naïve estimator for the vector is the difference in sample means around the changepoint location being estimated; Piotr calls this the Difference-In-Means Scaled cuSUM (DIM-SUM).

These simulations show that the DIM-SUM statistic is good at detecting dense changes in a panel, but breaks down when changes occur near the start or end of the panel. Two solutions are proposed:

* Re-scale to arrive at (multivariate) GLR statistic
* Use naïve estimator for projection direction normalized to have l2 norm 1. 

[1]: https://projecteuclid.org/euclid.ejs/1528941678
