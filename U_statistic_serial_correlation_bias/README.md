# Asymptotic Bias of Two Sample U-Statistic Under Serial Correlation

[He et al. (2020)][1] introduce a series of asymptotically independent U-statistics for high dimensional testing problems. The two sample difference in means test is particularly useful for changepoint problems. However, as shown by [Wang and Shao][2] without trimming these U-statistics are not appropriate for time series data as serial correlation between observations causes the statistics to be biased.

No indication of the degree of bias is given. In my theoretical investigations I found that for an a-th order U-statistic will have bias of order (n/T^{2a-1}),

* T = total sample size
* n = number of channels

meaning that the bias is asymptotically negligible, and the statistic is still usable, provided n grows sufficiently slowly compared to T. Simulations are run to confirm this.

[1]: https://arxiv.org/pdf/1809.00411.pdf
[2]: https://pdfs.semanticscholar.org/346a/fef94a0c0d27aec2b6f4dc3ba7a781ad1620.pdf
