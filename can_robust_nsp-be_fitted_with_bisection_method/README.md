# Can robust NSP be fits be performed with the bisection method?

The [robust Narrowest Significance Pursuit algorithm](https://stats.lse.ac.uk/fryzlewicz/nsp/rnsp.pdf) involves fitting a signal to data by minimizing the multi-resolution norm of empirical residual signs. For a vector having length n there are 2n + 1 (known) candidate constant functions which lead to distinct patterns of residual signs. The approach to function fitting in the paper is to try all candidate fits, using restricted multiresolution norm consisting of O(n) as opposed to O(n^2) intervals.

I claim it is only necessary to consider O(log n) candidate fits, since the function mapping from a candidate fit to the norm of residuals in weakly decreasing up to the arg-min and weakly increasing thereafter. The R script provided gives some empirical evidence for this claim. I also have a theoretical proof of the claim. 
