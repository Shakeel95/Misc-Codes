# Feature Alignment in High Dimensions

At each time step members of a high dimensional panel (n >> T) either contain a feature - mean shift, slope shift, etc. - or do not.
The code simulates the presence (1) or absence (0) of features directly; feature alignment is then investigated.

Two models are considered: a **random** model in which each time series has N features drawn independently and uniformly at random from {1,...,T} and a **factor** model in which each time series again has N features but the occur on a restricted set of common locations. Plots are produced to investigate, in particular, the maximum number of aligned features.
