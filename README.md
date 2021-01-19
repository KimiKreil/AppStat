# AppStat
General functions to use in the course Applied Statistics.
Developed during the course by: Kimi Kreilgaard, Jonathan Baltazar Steensgaard.

The ExternalFunctions module is the module provided by the course, with small adjustments and inclutions of functions found in the exercises. See the README_ext for more information.



Contains the functions below:


lprint(*args,**kwargs) 
--------------------------------------------------------------------------------------------------- 
 
Function provided by course to display symy results in latex 



weighted_avg(val, err, plot=False, title=None)
---------------------------------------------------------------------------------------------------
For problems asking us for the weighted average of a few values with respective errors, we need to 
- 1) calculate the average mean 
- 2) calculate its error 
- 3) make sure they are from a uniform distribution by calculating a chi2 value with a meaningful p-value. This is done by fitting with a constant function mu which has one parameter. The degrees of freedom will then be sample_size - 1.
- 4) If the p value is reasonable, our values agree.


val_err_contr(expr, **kwargs)
--------------------------------------------------------------------------------------------------

A function that takes in a paramater dependent on a set of variables and calculates the combined error analytically. For this to work we assume all of our errors are gaussian which is likely because of the central limit theorem. It also returns the contributions for each variable (the terms under the square root)

Uses the helper function: value_error_contribution_func_gen


binom_prob(r, n , p, plot=True)
---------------------------------------------------------------------------------------------------
This is a funcion that considers a problem following the binomial distribution, i.e. something that given a number of identical trials $n$ with two possible outcomes. One outcome is considered a succes with. probability $p$. The functions finds the probability of $r$ successes from $n$ trials given the probability $p$ of succes.


poisson_prob(r, lamb, plot=True)
---------------------------------------------------------------------------------------------------
This is a funcion that considers a problem following the poisson distribution, given the number of desired succeses r and the poisson parameter lambda, the function finds the probability of r successes occuring.


binom_trials(r, p, prob_r, test_range, plot=True)
---------------------------------------------------------------------------------------------------
If we know how sure we want to be of obtaining r number of succeses. How many trials are needed for the binomial distribution? Notice that both the probability mass function for the binomial AND the parameter n are discrete. This means we can't use a root finding algorithm for this problem, instead we use a binary search to narrow the windows of n's we will try until we find the solution.


poisson_trials(r, p, prob_r, guess = 1, plot=True)
---------------------------------------------------------------------------------------------------
If we know how sure we want to be of obtaining r number of succeses. How many trials are needed for the poisson distribution? Notice that both the probability mass function is discrete, however the parameter $\lambda$ from which we will find $n$ number of trials by $n = \frac{\lambda}{p}$. This means we can use a root finding algorithm to find the number of trials that provides a statistical guarantee of obtaining r succeses.


chauvenet_iterative(data, crit = 1/2)
---------------------------------------------------------------------------------------------------
If we are presented with real data, it is likely that it contains mismeasurement or outliers. Putting the data through this filter based on Chauvenets Criterium, we can identify the outliers and remove them from the data set. However, mismeasurements I believe have to be removed qualitatively.

More on the subject can be found here (also notes on which chauvenets factor to use): https://www.nbi.dk/~petersen/Teaching/Stat2017/Notes/StatNote_RejectingData.pdf

Essentially we take the point with the highest Z-score / furthest away from the mean, and see how many measurements we expect at that position, if we have under 1/2 (or another defined criteria) of what we expect, most likely we are working with the wrong probability distribution and we discard the data. We now repeat this proces until either no more points can be discarded with the criteria or if we reach a maximum number of points we want to discard. Say for example we will only discard 5% of the data, we stop if this is reached.


chauvenet_mask(data, crit=1/2)
---------------------------------------------------------------------------------------------------
If we are presented with real data, it is likely that it contains mismeasurement or outliers. Putting the data through this filter based on Chauvenets Criterium, we can identify the outliers and remove them from the data set. However, mismeasurements I believe have to be removed qualitatively.

More on the subject can be found here (also notes on which chauvenets factor to use): https://www.nbi.dk/~petersen/Teaching/Stat2017/Notes/StatNote_RejectingData.pdf

This function will check all data points at once against the criteria, and create a mask with False/0 for points that need to be discarded and True/1 for the points we want to keep.

The amount of measurement we expect a distance from the mean is:
$$ Prob = N \cdot \text{erfc} \left(\frac{|x_i - \mu|}{\sigma} \right) $$

Here we assume Gaussian, and thus use the error function.


sigma_cut(data, x)
---------------------------------------------------------------------------------------------------
This function will calculate the mean and std of a sample (assuming gaussian) and will then set a lower and upper limit at +/- x sigma away from the mean. It takes in the original data sample as an array and a number of sigma and returns an array with accepted data and one with the rejected data.


Sturges_bins(data)
---------------------------------------------------------------------------------------------------
A rule of thumb to estimate an appropriate number of bins. Sturges formula works well for normally distributed data. The formula can be found in https://en.wikipedia.org/wiki/Histogram#Sturges'_formula


Doanes_bins(data)
---------------------------------------------------------------------------------------------------
A two rule of thumb to estimate an appropriate number of bins. Doanes formula works well for skewed data. The formula can be found in https://en.wikipedia.org/wiki/Histogram#Sturges'_formula


KS_plotter(data, fit_function, args=(), zoom=True)
---------------------------------------------------------------------------------------------------
Calculates the Kolmogorov-Smirnov test and the associated p-value for a data set against a scipy stats function that has been fitted. An alternative to finding a goodness of fit that does not consider the number of bins.


KS_2samp_plotter(data1, data2, zoom=True, label1 = 'Cumulative Data 1', label2 = 'Cumulative Data 2')
---------------------------------------------------------------------------------------------------
Calculates the Kolmogorov-Smirnov test and the associated p-value for whether two data samples are drawn from the same distribution. An alternative to finding a goodness of fit that does not consider the number of bins.


chi2_2samp(data1, data2, N_bins=50, label1='Data 1', label2='Data 2', xlabel=None, txtcoord=(0.02, 0.95)
---------------------------------------------------------------------------------------------------
Performs the chi2 test on two data samples, comparing the two without fitting, the number of degrees of freedom is 
therefor the number of bins. Chi square is found by: np.sum( (counts1 - counts2)^2  / ( error1^2 + error2^2 ) )


find_C(fx_expr, xmin, xmax, all_sol = False)
---------------------------------------------------------------------------------------------------
Finds the normalising constant C for a function that we want to create monte carlo simulated data from. C can be either in the function, or one of the limits where the function is defined, for example x in [0,C].


find_invF(fx_expr_no_C, C_val=None, xmin=-oo, all_sol = False )
---------------------------------------------------------------------------------------------------
Finds the inverse function to a function f(x) we want to create our simulated monte carlo data from. This is for the transformation method.

accept_reject(func, N_points, xbound, ybound)
---------------------------------------------------------------------------------------------------
Performs the accept reject method of generating random numbers according to a distribution. The errors of the efficiency are explained in the link: http://www.pp.rhul.ac.uk/~cowan/atlas/efferr.ps

get_corr(x_data, y_data)
---------------------------------------------------------------------------------------------------
Calculates the linear correlation of two variables, given a data set with values for each variable.

error_rates(species_A, species_B, cut, N_bins = 50, plot=True, alp_coord=(0.1, 0.52), bet_coord=(0.1, 0.52), labelA = 'Species A', labelB='Species B')
---------------------------------------------------------------------------------------------------
Calculates the type I and II errors given a data set A and data set B and a cut to seperate the two.


gauss_hist(data, label=None, N_bins=None, ax=None)
---------------------------------------------------------------------------------------------------
Fits (chi2) a Gaussian to a histogram of the data and plots fit.


studentT_hist(data, label=None, N_bins=None, ax=None)
---------------------------------------------------------------------------------------------------
Fits (chi2) a Student T to a histogram of the data and plots fit.


double_gauss_hist(data, label=None, N_bins=None, ax=None, sub_pdfs = True, p0=())
---------------------------------------------------------------------------------------------------
Fits (chi2) a double gaussian to a histogram of the data and plots fit, along with the two sub gaussians.






