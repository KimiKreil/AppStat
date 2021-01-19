The Functions in external are:
------------------------------------------------------------------

lprint(*args, **kwargs)
------------------------------------------------------
For printing arguments with latex. 

Class: BinnedLH
------------------------------------------------------
__init__(self, f, data, bins=40, weights=None, weighterrors=None, bound=None, badvalue=1000000, extended=False, use_w2=False, nint_subdiv=1)


Class: UnbinnedLH()
------------------------------------------------------
__init__(self, f, data, weights=None, bound=None, badvalue=-100000, extended=False, extended_bound=None, extended_nint=100)

Class: Chi2Regression()
------------------------------------------------------
__init__(self, f, x, y, sy=None, weights=None, bound=None)


profile_x(x, y, bins=(50, 50), xyrange=[(0, 50), (-1,1)])
------------------------------------------------------
Gets the values and poisson errors from a 2d histogram. Used for calibration. 
Returns: x_center[mask], mean, std

Corr(mu1, sig1, mu2, sig2, rho12)
------------------------------------------------------
Generates a set of correlated randomly distributed gaussian numbers (from a seed of 42). 
Returns: x, y

calc_separation(x, y)
------------------------------------------------------
Calculating the seperation between two lists of numbers. This is a simple and quantitative way of doing so, and may be less precise than finding alpha and beta. 
Returns: d

FDA(spec_A, spec_B, labels=None)
------------------------------------------------------
Calculating the Fisher discriminant for seperating two sets of numbers, along with the Fisher Weights.
Returns: fisher_data_A, fisher_data_B, wf

calc_ROC(hist1, hist2):
------------------------------------------------------
Calculating the False Positive Rate and the True positive Rate, which can be used for a ROC Curve. 
Returns: FPR, TPR

plot_ROC(FPR, TPR, labels=None, colors=None, figsize=(8,6), ax=None):
------------------------------------------------------
Uses the arguments from the function before and plots it. This has been changed from ax. 
Returns: ax
