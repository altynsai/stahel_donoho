# Huberized Stahel-Donoho estimator
Most statistical models assume that there is some underlying distribution of the data and often give misleading results when some observations show anomalous
behavior. Before starting using these models, statisticians need to deal with the outliers causing bias in the estimates. Stahel-Donoho estimator is one of the
developed tools with great robustness properties aimed to detect outlying observations. The main idea is weighting observations depending on their maximum outlyingness on some univariate projections. However, the estimator is not effective when most observations are contaminated in at least one of the variables. 
S. Van Aelst, E. Vandervieren, and G. Willemsa in their article ‘A Stahel–Donoho estimator based on huberized outlyingness’ (2012) researched the possibilities and limitations of the SD estimator and suggested huberization of the data to overcome drawbacks of the original estimator. They suggest hubersization of the data before caculating the SD estimator which means bounding components value to reduce the influence of outliers.


