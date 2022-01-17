# Power curve simulations for local tests employed by NSP variants

The narrowest significance pursuit algorithm Fryzlewicz ([2021a](https://stats.lse.ac.uk/fryzlewicz/nsp/nsp.pdf), [2021b](https://stats.lse.ac.uk/fryzlewicz/nsp/rnsp.pdf)) works by applying local test for change points to subsets of a data vector. The specific local test used is a residual based test based on the multi-resolution sup norm; however, other tests can also be used. One option is construct non-parametric confidence sets for a supposed regression function, then use these sets to test for a change point locally. See for example [here](https://github.com/Shakeel95/CMStatistics-2021) and [here](https://github.com/Shakeel95/Stats-Experiments/tree/master/robust-nsp-using-runs).

Here I investigate the power of four local tests which may be used by the NSP algorithm, and which are particularly useful when the analyse would like to test for change points "robustly". The specific tests are: 

* The test in Fryzlewicz (2021b) based on the multi-resolution norm of residual signs. 
* A test based on confidence sets obtained by inverting the runs test for noise. 
* A test based on confidence sets obtained by inverting the multi-resolution norm of residual signs. 
* Same as above, with a scale correction as introduced by [Dumbgen and Spokoiny (2001)](https://projecteuclid.org/journals/annals-of-statistics/volume-29/issue-1/Multiscale-Testing-of-Qualitative-Hypotheses/10.1214/aos/996986504.full)