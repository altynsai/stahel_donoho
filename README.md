# Huberized Stahel-Donoho estimator
Most statistical models assume that there is some underlying distribution of the data and often give misleading results when some observations show anomalous behavior. Before starting using these models, statisticians need to deal with the outliers causing bias in the estimates. Stahel-Donoho estimator is one of the developed tools with great robustness properties aimed to detect outlying observations. The main idea is weighting observations depending on their maximum outlyingness on some univariate projections. However, the estimator is not effective when most observations are contaminated in at least one of the variables. 

S. Van Aelst, E. Vandervieren, and G. Willemsa in their article ‘A Stahel–Donoho estimator based on huberized outlyingness’ (2012) researched the possibilities and limitations of the SD estimator and suggested huberization of the data to overcome drawbacks of the original estimator. They suggest hubersization of the data before caculating the SD estimator which means bounding components value to reduce the influence of outliers.

In this research we are running similar simultation study. To illustrate the efficiency of both estimators we are going to generate some toy datasets the same way as described in the article. For convenience observations are sampled from multivariate normal distributions with mean 0 and variance 1. In this section we have 5-dimensional data with 2 contaminated variables.

Performance of the estimators might depend on the characteristics of the data. First, we want to compare robustness and precision of the estimators for independent (ρ = 0) and highly correlated data (ρ = 0. 9). Second, we change the percentage of outliers in the component (𝜖 ∈ {0. 1, 0. 2}). The outlying values are generated from univariate normal distribution with mean 𝑘/ 𝑑, where 𝑘 is outlying distance and 𝑑 is the number of contaminated variables (in our case, 2). Standard deviation is equal to 0.1. We randomly replace the given percent of original data with outlying values. With our chosen parameters, the maximum number of outliers we can have in the data is 40% which is smaller than the breakpoint of the SD estimator (50%).

After the dataset is formed, we calculate the classic Stahel-Donoho estimator. In accordance with the article suggestion we sample 200 * 𝑝 = 200 * 5 = 1000 directions in order to calculate the outlyingness of the observations. Estimates of mean and variance-covariance matrix are computed using the Huber-type weight function with cut-off value 𝑐 = 𝑚𝑖𝑛( χ . The same procedure takes place 𝑝 2(0. 5), 4) for the huberized dataset. For comparison of the estimators we compute the ratio of mean squared errors for location and scale estimators. In order to get consistent results, we sample 500 observations and take the ratio of average MSE of both estimators. MSE for all components and only contaminated ones are calculated. Off-diagonal elements of the variance-covariance matrix are divided into additional groups whether correlation between one or two contaminated variables is estimated. 

In Table 1, we can see in which cases the Huberized SD-estimator is more precise than the classic one:
<img width="500" alt="table 1" src="https://user-images.githubusercontent.com/101756813/203383234-83523ad1-b1a5-4cbd-a1d7-c9f40bde5443.png">

No doubt the modified estimator cannot deal well with highly correlated data. With an increase of outliers percentage from 10% to 20% we observe improvement of the HSD-estimator. This improvement is especially visible for big outlying distances and non-correlated data. The MSE ratios are slightly smaller when we only consider contaminated variables. Overall, there is a small difference between the estimators for data with a low percentage of outliers, but the difference grows with this fraction.

In Figure 1, we plot the mean absolute errors for location estimators. We compare them for contaminated and non-contaminated components. The top boxplots show results for 10% of outliers, while bottom ones - for 20%. In general, we do not observe large errors for the estimators and they have great performance for low dimensional data with limited number of contaminated components.

<img width="500" alt="image" src="https://user-images.githubusercontent.com/101756813/203383351-0eb594fb-6438-442a-8282-4eae3fd6281a.png">

Our results are completely in line with conclusions made in the article. No significant improvement takes place as the fraction of outlying data is smaller than the breakdown point. We expect that difference between the estimators increases in the next section with a bigger number of contaminated components.

Following the same procedure as in Table 1, we analyze the results of our Huberized estimator with a contamination rate of 20%, varying levels of correlation in the data, and variable distance of our contaminated components from our uncontaminated data. We have seven dimensions, five of which were contaminated in the previously described manner. In this example we intentionally create data with outlier propagation to exceed the breakdown point of 50%.

We note that for large distance outliers (k=64) Huberized estimator is more effective than standard one. Such improvement is observed for both correlated and independent data. However, it is worth noting that for independent data the ratio is smaller, and this means that for such data MSE of HSD-estimator drops faster than of standard SD-estimator.

For small distance outliers (k=6) the estimators do not differ significantly. Overall, the modified estimator is strictly superior to the original SD estimator by this table in nearly all cases.

In the following table we summarize the results:

<img width="500" alt="image" src="https://user-images.githubusercontent.com/101756813/203384058-d0d6547b-b876-4fbc-93c6-c598bdf58146.png">
