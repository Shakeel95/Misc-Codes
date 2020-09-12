# Most Prominent Trend Change

In a high dimensional panel consisting of piecewise linear signals with possibly aligned kink discontinuities (changepoints), it is useful to have a method for ranking each changepoint by its *importance*.

## Definition of importance

Take a changepoint detection statistic T acting on an interval [s,e], changepoints could be ranked by importance as follows.

Imagine each changepoint has been isolated to its own suitably large interval *. For each interval use T to test the null of no changepoint against the alternative of a single changepoint. The most prominent / important changepoint is the one whose test provides the most evidence against the null, and so on...

## Choice of test

How to choose T? For a univariate signal [Baranowski et al (2019)][1] use the generalized likelihood ratio statistic (GLR) as a building block. To extend to the panel data setting I investigate aggregating GLRs from each channel of the panel using:

1. L-2 aggregation following [Horvath and Huskova (2012)][2]
2. L-infinity aggregation following [Jirak (2015)][3]
3. 'Scan' aggregation following [Enikeeva and Harachoui (2013)][4]
4. L-2 aggregation with self-normalization following [Shao and Zhang (2010)][5]

## Simulations

Simulations are run to check whether each aggregation method performs well (✔) or poorly (❌) in the presence of the features which are know to make multivariate changepoint detection tricky:  

| aggregation method |  changepoint location | change density | size of change | spatial dependence | serial dependence |
|--------------------|-----------------------|----------------|----------------|--------------------|-------------------|
| L-2                |                       |                |                |                    |                   |
| L-Infinity         |                       |                |                |                    |                   |
| 'Scan'             |                       |                |                |                    |                   |
| Self-Normalization |                       |                |                |                    |                   |

\* To construct the most balanced intervals one possible approach is to isolate each changepoint in an interval so that the sum of squared lengths of the intervals in minimal.

[1]: https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12322
[2]: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9892.2012.00796.x
[3]: https://projecteuclid.org/euclid.aos/1444222081
[4]: https://arxiv.org/abs/1312.1900
[5]: http://www.stat.tamu.edu/~zhangxiany/JASA-2010.pdf
