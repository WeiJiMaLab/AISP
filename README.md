# Point estimate observers
This is the code and data to accompany our manuscript on point estimate observers.

The code here is MATLAB code, which depends only on the code linked into this repository. Simply downloading should be sufficient.

To run the code you will need the data and parameter fits, which are available at [osf](https://osf.io/x8q6j/) . To use, extract the data and parameter files into their respective code folders.
In each code folder you find `*_simulate.m` files that implement the different observer models, `ibslike` functions that compute likelihoods using inverse binomial sampling, `fit_cluster_ibs` that run BADS to fit parameters and diverse plotting functions that illustrate results.
