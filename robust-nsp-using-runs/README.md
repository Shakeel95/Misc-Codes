# Robust NSP using runs

The narrowest significance pursuit algorithm [Fryzlewicz (2021a)][1] works by applying local tests for change points to subsets of a data vector. The specific local test used is a residual based test based on the multi-resolution sup norm. In the robust setting, a similar test based on the multi-resolution norm of residual signs is proposed by [Fryzlewicz (2021b)][2] however the test statistic is only efficiently computable in the "piecewise constant" setting and as a result unlike the original NSP algorithm the robust variant cannot be used to test for change points in higher order polynomials.

Here I investigate replacing the local test by a procedure which involves computing a confidence set for the regression function by inverting the runs test - see [Davies and Kovac (2001)][3] - then checking whether the confidence set contains a plausibly change point free regression function. Since the method for obtaining the confidence set is non-parametric the test easy to apply to change points in higher order polynomials.

[1]: https://stats.lse.ac.uk/fryzlewicz/nsp/nsp.pdf
[2]: https://stats.lse.ac.uk/fryzlewicz/nsp/rnsp.pdf
[3]: https://projecteuclid.org/journals/annals-of-statistics/volume-29/issue-1/Local-Extremes-Runs-Strings-and-Multiresolution/10.1214/aos/996986501.full
