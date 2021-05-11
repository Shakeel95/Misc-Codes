# Eigenvalue Test for Misaligned Changepoints in a Panel

Most algorithms for changepoint detection in multivariate settings assume aligned changes, but this assumption is unlikely to hold in practice: in real data settings channel-wise changes tend to occur at points in time which are close but not identical. Here a simulation study is presented to investigate whether existing state of the art methods are able to perform well in this realistic setting.

Consider the canonical setting where data follows a signal + noise model and the signal component is a piecewise constant function. The aim will be to use a changepoint algorithm to fully eliminate the signal by estimating the changepoint locations and de-meaning between these locations. To check whether the signal has been fully eliminated the first eigenvalue of the sample covariance matrix for the de-meaned data can be inspected. This is because the presence of a mean shift in multivariate data is leads to a rank one update of the covariance matrix; if the signal component has not been fully eliminated we expect to see a spike in the first eignevalue. 

The following methods, which have CRAN implementations which allow for fast multiple changepoint detection in multivariate data are considered:

* [Sparse projection](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12243?campaign=wolearlyview)
* [Double CUSUM binary segmentation](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-2/Change-point-detection-in-panel-data-via-double-CUSUM-statistic/10.1214/16-EJS1155.full)
* [Sparsified binary segmentation](https://www.jstor.org/stable/24774746?seq=1#metadata_info_tab_contents)
* [Geometrically inspired mapping](https://link.springer.com/article/10.1007/s11222-020-09940-y)
