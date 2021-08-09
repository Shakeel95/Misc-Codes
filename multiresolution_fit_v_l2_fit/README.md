# Fitting linear models under the l2 loss and multi-resolution sup-norm loss

[Fryzlewicz (2020)](https://stats.lse.ac.uk/fryzlewicz/nsp/nsp.pdf) introduces a test for changepoints in linear models by comparing the "size" of residuals from a fit under the multi-resolution sup norm. How much do the fitted parameters differ from those which would have been obtained under the l2 loss?

This experiment confirms that even in the absence of noise the fitted parameters under the two loss functions will always be always O(1) apart. The only exception is the piecewise constant setting with a single jump, where the estimated parameters are identical. 
