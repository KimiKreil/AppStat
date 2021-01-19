#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:22:28 2020

@author: Kimi Kreilgaard, Jonathan Baltazar Steensgaard and Amalie Beate Albrechtsen

"""

import numpy as np                                    
import matplotlib.pyplot as plt

import seaborn as sns                                  
from iminuit import Minuit                             
import sys                                             
from scipy import stats, optimize
from scipy.optimize import minimize

from sympy import * 
from sympy import sympify

from ExternalFunctions import *

from IPython.display import display
from IPython.core.display import Latex

from matplotlib.patches import ConnectionPatch

def lprint(*args,**kwargs):
    """Pretty print arguments as LaTeX using IPython display system 
    
    Parameters
    ----------
    args : tuple 
        What to print (in LaTeX math mode)
    kwargs : dict 
        optional keywords to pass to `display` 
    """
    display(Latex('$$'+' '.join(args)+'$$'),**kwargs)

def weighted_avg(val, err, plot=False, title=None):
    """
    INPUT:
    val = values, array_like
    err = erros, array_like
    plot = option to plot or not
    
    """
    
    # Calculate the avg according to Barlow (4.6)
    avg = np.sum( (val / err**2) / np.sum( 1 / err**2 ) )
    
    # Calculate the error
    avg_sig = np.sqrt( 1 / np.sum(1 / err**2) ) 
    
    # Find degrees of freedom (-1 )
    N_dof = len(val) - 1
    
    # Calculate chi_square
    chi2 = np.sum( (val - avg)**2 / err**2 )
    
    # Calculate p-value (the integral of the chi2-distribution from chi2 to infinity)
    p = stats.chi2.sf(chi2, N_dof)
    
    # Option to plot the fitted line
    if plot:
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12,6))
        
        # X values are measurement number
        x = np.arange(len(val))+1
        
        # Plot values and errorbars
        ax.scatter(x, val)
        ax.errorbar(x, val, err, fmt='ro', ecolor='k', elinewidth=1, capsize=2, capthick=1)
        
        #Plot the weighted average line
        ax.hlines(avg, 0, len(val)+0.5, colors = 'r', linestyle = 'dashed')
        
        # Nice text
        d = {'mu':   avg,
             'sigma_mu': avg_sig,
             'Chi2':     chi2,
             'ndf':      N_dof,
             'Prob':     p,
            }

        text = nice_string_output(d, extra_spacing=2, decimals=3)
        add_text_to_ax(0.02, 0.95, text, ax, fontsize=14)
        ax.set_title(title, fontsize=18)
        fig.tight_layout()

    return avg, avg_sig, chi2, p




def value_error_contribution_func_gen(expr, variables):
    """
    expr = takes in a math expression in a string of type 'a+b'
    var = takes in a tuple of variables strings, fx ('a', 'b') 
    """
    # Convert expression into a sympy expression
    expr = parse_expr(expr)
    
    # Define sympy symbols for the parameters (the tuple variables) and the standard deviations
    var_symbols = symbols(variables)
    err_symbols = symbols( tuple("sigma_" + k for k in variables) )
    
    # Find expressions for each contributions
    contributions = [expr.diff(var) ** 2 * err**2 for var, err in zip(var_symbols, err_symbols)]
    
    # Convert contributions to numerical functions
    f_contributions = [ lambdify(var_symbols + err_symbols, expression) for expression in contributions ]

    # Find the error propagation expression to be evaluated, and display
    expr_sig = sqrt( sum(contributions) )
    display(expr_sig)
    
    # Convert the expression for the value and the error into numerical functions
    f_val = lambdify(var_symbols, expr)
    f_err = lambdify(var_symbols + err_symbols, expr_sig)
    
    def func(**kwargs):
        """
        Define a function that will take in keywordarguments **kwargs which is a dictionary of type: 
        {'a':(1,0.1), 'b':(2,0.3)}. Kwargs.values calls the two tuples as one list [(1,0.1),(2,0.3)].
        From there an array of variables and an array of errors can be extracted and the numerical
        functions found above can be used.
        
        """
        # Create tuple of values of variables
        v = tuple(v[0] for v in kwargs.values())
        
        # Create tuple of errors of variables
        s = tuple(v[1] for v in kwargs.values())
        
        # Calculate value and error
        value, error = f_val(*v), f_err(*v, *s)

        # Calculate contribution from each variable
        contr_list = [ function(*v,*s) for function in f_contributions ]
        
        #Return value and analytical error
        return value, error, contr_list
    
    # Return the main function that we set out to generate
    return func


# Define function that gets variables from **kwargs and uses the function above to return value and error
def val_err_contr(expr, **kwargs):
    """
    INPUT:
    expr = takes in a math expression in a string of type 'a+b'
    **kwargs = variable names = (value, error) of type a=(3, 0.3)
    
    Note that if the relation depends on constant, type those in as variables with sigma = zero.
    
    
    OUTPUT:
    value = integer
    error = integer
    contributions = array_like with contributions from each variable in the same order as in the input
    """
    
    return value_error_contribution_func_gen(expr, tuple(kwargs))(**kwargs)




def binom_prob(r, n , p, plot=True):
    """
    INPUT: r = number of succeses, n = number of trials, p = probability of succes.
    
    OUTPUT: prob_r = Probability of r
    """
    
    # Calculate probability of r succeses, and print
    prob_r = stats.binom.pmf(r, n, p)
    print(f'The probability of {r:.0f} succeses is {prob_r:.4f}')
    
    # Option to plot
    if plot:
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Make an array of r's to evaluate the probability of (centered around the mean)
        mean = n*p
        variance = n*p*(1-p)
        
        # Plot
        xaxis = np.linspace( mean - variance, mean + variance, 1000)
        yaxis = stats.binom.pmf(np.round(xaxis), n, p)
        ax.plot(xaxis, yaxis, '-', label = 'Binomial pmf', zorder=1, color='k')
        
        # Nice text
        d = {'Number of trials': n,
             'Prob of succes': p,
             'Mean': mean,
             'Std': np.sqrt(variance),}
        
        text = nice_string_output(d, extra_spacing=2, decimals=3)
        add_text_to_ax(0.02, 0.87, text, ax, fontsize=12)
        ax.set(title='Binomial Probability', xlabel='Number of succeses r', ylabel='Probability')
        
        ax.scatter(r, stats.binom.pmf(r,n,p), color = 'r', label=f'Probability of {r} successes', zorder=2)
        ax.legend(loc='upper left')

    return prob_r


def poisson_prob(r, lamb, plot=True):
    """
    INPUT: r = number of successes, lamb = poisson parameter (n*p).
    
    OUTPUT: prob_r = Probability of r
    """
    
    # Calculate probability of r succeses, and print
    prob_r = stats.poisson.pmf(r, lamb)
    print(f'The probability of {r:.0f} succeses is {prob_r:.4f}')
    
    # Option to plot
    if plot:
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Make an array of r's to evaluate the probability of (centered around the mean)
        mean = lamb
        variance = lamb
        
        # Plot
        xaxis = np.linspace( mean - variance, mean + variance, 1000)
        yaxis = stats.poisson.pmf(np.round(xaxis), lamb)
        ax.plot(xaxis, yaxis, '-', label = 'Poisson pmf', zorder=1, color='k')
        
        # Nice text
        d = {'Lambda/Mean': mean,
             'Std': np.sqrt(variance),}
        
        text = nice_string_output(d, extra_spacing=2, decimals=3)
        add_text_to_ax(0.02, 0.87, text, ax, fontsize=12)
        ax.set(title='Poisson Probability', xlabel='Number of succeses r', ylabel='Probability')
        
        ax.scatter(r, stats.poisson.pmf(r, lamb), color = 'r', label=f'Probability of {r} successes', zorder=2)
        ax.legend(loc='upper left')

    return prob_r


    
        
def binom_trials(r, p, prob_r, test_range, plot=True):
    """
    INPUT: 
    r = number of succeses we want guaranteed, 
    p = probability of succes in each trial
    prob_r = statistically guaranteed probability of getting at least r number of succeses in n trials
    test_range = (low, high), limits of the arrays of n's we will test
           
    OUTPUT: 
    n = number of trials to ensure the prob_r
    """
    
    # Create an array of n trials to test
    test_array = np.arange(*test_range, 1)
    
    # We will reduce the array to something bigger than len(3) through binary search
    while len(test_array) > 3:
        
        # Choose the median of the array 
        middle_index = np.int( np.floor( len(test_array)/2 ) ) + 1
        n_median = test_array[middle_index]
                              
        # Test the survival function of this n. Notice we use r-1 since we want all successes > r-1
        sf = stats.binom.sf(r-1, n_median, p)
        
        # If the value is under 90% we look in the second half of the array since we need more trials
        # Else we look in the lower half of the array
        if sf < prob_r: 
            test_array = test_array[middle_index:]
            
        if sf > prob_r: 
            test_array = test_array[:middle_index]
    
    # We will now start from the first element (the lowest n) in the test array and calculate the sf
    # the moment it surpasses prob_r we return that n
    n = test_array[0]
    sf = stats.binom.sf(r-1, n, p)

    while sf < prob_r:
        n += 1
        sf = stats.binom.sf(r-1, n ,p)
    
    # If we obtain the guaranteed succesrate, return the given n
    if sf >= prob_r:
        
        # Print result
        print(f'To be {prob_r*100:.0f}% sure of {r:.0f} success(es), we need {n} trials')
        
        # The option to plot
        if plot:
            
            fig, ax = plt.subplots(ncols=2, figsize=(15,6), sharex=True, 
                                   gridspec_kw = {'width_ratios': [1.2, 1.5], 'wspace': 0.2})
            fig.suptitle(f'Binomial Distribution with n={n} and p={0.054}', fontsize=18)
            
            ### PLOT BINOMIAL PMF ###
            
            # Make an array of r's to evaluate the probability of (centered around the mean)
            mean = n*p
            variance = n*p*(1-p)
        
            xaxis = np.linspace( mean - variance, mean + variance, 1000)
            yaxis = stats.binom.pmf(np.round(xaxis), n, p)
            
            ax[0].plot(xaxis, yaxis, '-', label = 'Binomial pmf', zorder=1, color='k')
            ax[0].scatter(r-1, stats.binom.pmf(r-1,n,p), color='r', zorder=2)
            ax[0].set(title='Probability Mass Function', 
                      xlabel='Number of succeses r', ylabel='Probability')
        
            # Define the vertical line that seperates the bins we dont want to end with the ones we want
            ax[0].fill_between(x = [0, r-0.52], y1 = 0, y2 = 1, color='red', alpha=0.1)
            ax[0].fill_between(x = [r-0.5, xaxis[-1]], y1 = 0, y2 = 1, color='green', alpha=0.1)
            
            ### PLOT SURVIVAL FUNCTION ###
            ax[1].plot(xaxis, stats.binom.sf(np.round(xaxis), n, p), zorder=1, color='k')
            ax[1].scatter(r-1, stats.binom.sf(r-1,n,p),color='r',zorder=2)
            ax[1].set(title='Survival Function', xlabel='Number of succeses r', 
                      ylabel='Probability of r or more successes')
            
            # Define horisontal line that seperates the bins of confidence we are interested in
            ax[1].fill_between(x = [0, xaxis[-1]], y1 = 0.9, y2 = 1, color='red', alpha=0.1)
            ax[1].fill_between(x = [0, xaxis[-1]], y1 = 0, y2 = 0.898, color='green', alpha=0.1)
            
        return n
    
    
    
def poisson_trials(r, p, prob_r, guess = 1, plot=True):
    """
    INPUT: 
    r = number of succeses we want guaranteed, 
    p = probability of succes in each trial
    prob_r = statistically guaranteed probability of getting at least r number of succeses in n trials
    test_range = (low, high), limits of the arrays of n's we will test
           
    OUTPUT: 
    n = number of trials to ensure the prob_r
    """
    
    # Solve for lambda that gives at least prob_r success 
    #Notice we use r-1 since we want all successes > r-1
    root = optimize.root(lambda lamb: stats.poisson.sf(r-1, lamb) - prob_r, guess)
    
    # If succesful, continue
    if root.success:
        
        # Extract solution
        lamb = root.x[0]
        
        # Calculate number of trials, and print. We ceal since we need the next integer.
        n = np.int ( np.ceil(lamb / p) )
        print(f'To be {prob_r*100:.0f}% sure of {r:.0f} success(es), we need {n} trials')
        
        # Option to plot
        if plot:
            
            fig, ax = plt.subplots(ncols=2, figsize=(15,6), sharex=True, 
                                   gridspec_kw = {'width_ratios': [1.2, 1.5], 'wspace': 0.2})
            fig.suptitle(f'Poisson Distribution with lambda={lamb:.1f}, p={0.054} and thus n={n}', fontsize=18)
            
            ### PLOT POISSON PMF ###
            
            # Make an array of r's to evaluate the probability of (centered around the mean)
            mean, variance = lamb, lamb
        
            xaxis = np.linspace( mean - variance, mean + variance, 1000)
            yaxis = stats.poisson.pmf(np.round(xaxis), lamb)
            
            ax[0].plot(xaxis, yaxis, '-', label = 'Poisson pmf', zorder=1, color='k')
            ax[0].scatter(r-1, stats.poisson.pmf(r-1,lamb),color='r', zorder=2)
            ax[0].set(title='Probability Mass Function', 
                      xlabel='Number of succeses r', ylabel='Probability')
            
            # Define the vertical line that seperates the bins we dont want to end with the ones we want
            ax[0].fill_between(x = [0, r-0.52], y1 = 0, y2 = 1, color='red', alpha=0.1)
            ax[0].fill_between(x = [r-0.5, xaxis[-1]], y1 = 0, y2 = 1, color='green', alpha=0.1)
            
            ### PLOT SURVIVAL FUNCTION ###
            ax[1].plot(xaxis, stats.poisson.sf(np.round(xaxis), lamb), zorder=1, color='k' )
            ax[1].scatter(r-1, stats.poisson.sf(r-1, lamb),color='r',zorder=2)
            ax[1].set(title='Survival Function', xlabel='Number of succeses r', 
                      ylabel='Probability of r or more successes')
            
            # Define horisontal line that seperates the bins of confidence we are interested in
            ax[1].fill_between(x = [0, xaxis[-1]], y1 = 0.9, y2 = 1, color='red', alpha=0.1)
            ax[1].fill_between(x = [0, xaxis[-1]], y1 = 0, y2 = 0.898, color='green', alpha=0.1)
            
        return n
    
    # If unsuccesful we probably need another guess
    else: print('A solution could not be obtained. Try another guess.')
    



def chauvenet_iterative(data, crit = 1/2):
    
    """
    INPUT:
    data = 1d array with all data
    crit = float, with crit, default = 1/2, "equal quantity of points close to and far from mean"
    
    OUTPUT:
    data = 1d array with new cleaned data
    rejected_data = 1d array with rejected data
    """
    
    # Define an initial value of Prob, to create a while loop over
    prob = 0.01
    
    # Empty list to save the rejected data
    rejected_data = []
    
    # Start while loop which will rejects points under the criteria
    while prob < crit:
        
        # Calculate the mean of the sample, and the standard deviation with Bessels Correction
        mean, sigma = np.mean(data), np.std(data, ddof=1)

        # Calculate the absolute Z_scores
        Z_abs = abs(data - mean) / sigma
        
        # Find highest Z score - the point furthest away that we will investigate first
        index, Z = np.argmax(Z_abs), np.max(Z_abs)
        
        prob = len(data) * erfc(Z)
        
        # If the point is rejected, remove it from the data set and save it to rejected data
        if prob < crit: 
            rejected_data = np.append( rejected_data, data[index] )
            data = np.delete( data, index )
            
    return data, rejected_data




def chauvenet_mask(data, crit=1/2):
    
    """
    INPUT:
    data = 1d array with data
    crit = chauvenets criterion
    """
    
    # Find mean and standard deviation
    mean, sigma = np.mean(data), np.std(data, ddof=1)
    
    # Define criterion
    criterion = crit
    
    # Z_score: number of standard deviations the point is away from the mean
    Z_abs = abs(data - mean) / sigma 
    
    # Area under gaussian
    prob = len(data) * erfc(Z_abs)
    
    # Return a mask that is true for values we keep
    return prob > criterion




def sigma_cut(data, x):
    
    """
    INPUT:
    data = 1d array with all data
    x = float, number of sigma we want the cut off to be placed at
    
    OUTPUT:
    new_data = 1d array with accepted data
    rejected_data = 1d array with rejected data
    """

    # Find the mean and  standard deviation
    mean, sigma = np.mean(data), np.std(data, ddof=1)

    # Find interval we wish to keep
    low = mean - x*sigma
    high = mean + x*sigma

    # Cut the array to find new data
    new_data = data[ data >= low]
    new_data = new_data[ new_data <= high]
    
    # Collect rejected data
    rejected_data_low = data[ data < low ]
    rejected_data_high = data[ data > high ]
    
    rejected_data = np.concatenate( (rejected_data_low, rejected_data_high) )
    
    return new_data, rejected_data




def Sturges_bins(data):
    
    """
    INPUT:
    data = 1d array of all data
    
    OUTPUT:
    k = number of bins to use
    """
    
    # Number of data
    N = len(data)
    
    # Number of bins
    k = np.int( np.ceil( np.log2(N) ) + 1 )
    
    return k



def Doanes_bins(data):
    
    """
    INPUT:
    data = 1d array of all data
    
    OUTPUT:
    k = number of bins to use
    """
    
    # Number of data
    N = len(data)
    
    # Skewness
    g = 1 / (N * np.std(data)) * np.sum( (data - np.mean(data) )**3 )
    
    sig_g = np.sqrt( (6*(N-2)) / ((N+1)*(N+3)) )
    
    # Number of bins
    k = np.int( np.ceil( 1 + np.log2(N) + np.log2(1 + abs(g) / sig_g) ))
    
    return k



def KS_plotter(data, fit_function, args=(), zoom=True):
    """
    So far this function only works to compare empirical data to a scipy stats function,
    and has not yet been expanded to include user fit functions.

    INPUT:
    data = empirical data to test against a fit
    fit_function = Scipy stats function that is fitted with in a string, fx 'norm'
    args = arguments to give fit function in default order
    zoom = option to zoom in on supremum

    OUTPUT:
    ks_stat = ks statistics value
    pval = p value of KS test
    """
    
    # Create figure
    fig, ax = plt.subplots(nrows=2, figsize=(12,8), sharex=True, 
                       gridspec_kw = {'height_ratios': [1.5, 0.4], 'hspace': 0})
    
    # Plot data
    xaxis = np.sort(data)
    yaxis = np.arange(0, len(data), 1)
    ax[0].plot(xaxis, yaxis, color='red', label='Cumulative Data')
    
    # Plot fitted function's cdf
    cdf = getattr(stats.distributions, fit_function).cdf
    cdf_axis = len(data) * cdf(xaxis, *args)
    ax[0].plot(xaxis, cdf_axis, color='k', label = fit_function + '' + 'CDF')
    ax[0].legend(loc='upper left', prop={"size":14} )
    
    # Plot residuals
    resi = yaxis - cdf_axis
    ax[1].plot(xaxis, resi, color='k', linewidth=0.6, label = 'Residuals (Data $-$ CDF)')
    ax[1].hlines(0, xaxis[0], xaxis[-1], linestyle='dashed', color='r')
    
    ax[1].set_ylim(min(resi)+0.1*min(resi), max(resi)+0.1*max(resi))
    ax[1].legend(loc='lower right', prop={"size":12})

    # Plot zoom
    if zoom:
        
        # Create extra axis
        ax1 = fig.add_axes([0.65, 0.35, 0.2, 0.25]) # add_axes([x0, y0, width, height])
        ax1.plot(xaxis, yaxis, color='red')
        ax1.plot(xaxis, cdf_axis, color='k')
        
        # Find and supremum
        index = np.argmax(abs(resi))
        
        # Adjust limits
        if resi[index] > 0: # Empirical data highest
            ymin, ymax = cdf_axis[index] - resi[index], yaxis[index] + resi[index]
            
        if resi[index] < 0: # Fit is highest
            ymin, ymax = yaxis[index] + resi[index], cdf_axis[index] - resi[index]
        
        xmin, xmax = xaxis[index] - 0.001 * xaxis[index], xaxis[index] + 0.001 * xaxis[index]
        
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin, ymax)
        
        # Mark supremum
        supremum = ConnectionPatch(xyA=(xaxis[index], cdf_axis[index]), xyB=(xaxis[index], yaxis[index]), coordsA=ax1.transData, 
                           arrowstyle='<->', color='b')
        fig.add_artist(supremum)

        ax1.set_title('$D_n = sup_x |F_n(x)-F(x)| $', color='b', fontsize=14)
        
        # Add zoom lines
        con1 = ConnectionPatch(xyA=(xmin, ymin), coordsA=ax[0].transData, xyB=(xmin, ymin), coordsB=ax1.transData, alpha=0.3)
        con2 = ConnectionPatch(xyA=(xmax, ymax), coordsA=ax[0].transData, xyB=(xmax,ymax), coordsB=ax1.transData, alpha=0.3)

        sq1 = ConnectionPatch(xyA=(xmin, ymin), xyB=(xmax, ymin), coordsA=ax[0].transData, alpha=0.3)
        sq2 = ConnectionPatch(xyA=(xmin, ymax), xyB=(xmax, ymax), coordsA=ax[0].transData, alpha=0.3)
        sq3 = ConnectionPatch(xyA=(xmin, ymin), xyB=(xmin, ymax), coordsA=ax[0].transData, alpha=0.3)
        sq4 = ConnectionPatch(xyA=(xmax, ymin), xyB=(xmax, ymax), coordsA=ax[0].transData, alpha=0.3)

        fig.add_artist(con1)
        fig.add_artist(con2)
        fig.add_artist(sq1)
        fig.add_artist(sq2)
        fig.add_artist(sq3)
        fig.add_artist(sq4)
        
    # Perform ks test
    ks_stat, pval = stats.kstest(data, fit_function, args=args )

    # Adding fit results to plot:
    d = {'KS stat':     ks_stat, 'Prob':     pval,}
    
    text = nice_string_output(d, extra_spacing=2, decimals=3)
    add_text_to_ax(0.02, 0.82, text, ax[0], fontsize=14)

    plt.show()
    
    return ks_stat, pval


def KS_2samp_plotter(data1, data2, zoom=True, label1 = 'Cumulative Data 1', label2 = 'Cumulative Data 2'):
    
    """
    INPUT:
    data1 = 1d arraylike, with data sample 1
    data2 = 1d arraylike, with data sample 2
    
    OUTPUT:
    ks_stat = ks statistics value from a two sample test comparing the two data sets
    pval = p value of KS test
    """
    
    # Get data stuff to be compatible (interpolate to same size) and make cdf array
    data1 = np.sort(data1)
    data2 = np.sort(data2)
    data_all = np.sort( np.concatenate([data1, data2]) )
    cdf1 = np.searchsorted(data1, data_all, side='right')
    cdf2 = np.searchsorted(data2, data_all, side='right')

    # Create figure
    fig, ax = plt.subplots(nrows=2, figsize=(12,8), sharex=True, 
                       gridspec_kw = {'height_ratios': [1.5, 0.4], 'hspace': 0})
    
    # Plot the two cdfs
    ax[0].plot(data_all, cdf1, color='r', label=label1)
    ax[0].plot(data_all, cdf2, color='b', label=label2)
    ax[0].legend(loc='upper left', prop={"size":14} )
    
    # Calculate and plot residuals
    resi = cdf1 - cdf2
    ax[1].plot(np.sort(data_all), resi, color='k', linewidth=0.6, label='Residuals (Data1 - Data2)')
    ax[1].hlines(0, data_all[0], data_all[-1], linestyle='dashed', color='k')
    ax[1].legend(loc='lower left', prop={"size":12})
    
    # Plot zoom
    if zoom:
        
        # Create extra axis
        ax1 = fig.add_axes([0.65, 0.35, 0.2, 0.25]) # add_axes([x0, y0, width, height])
        ax1.plot(data_all, cdf1, color='r')
        ax1.plot(data_all, cdf2, color='b')
        
        # Find and supremum
        index = np.argmax(abs(resi))
    
        # Adjust limits
        if resi[index] > 0: # Data1 is highest at supremum
            ymin, ymax = cdf2[index] - resi[index], cdf1[index] + resi[index]
            
        if resi[index] < 0: # Data2 is highest at supremum
            ymin, ymax = cdf1[index] + resi[index], cdf2[index] - resi[index]
        
        xmin, xmax = data_all[index] - 0.05 * data_all[index], data_all[index] + 0.05 * data_all[index]
        
        ax1.set_xlim(xmin, xmax)
        ax1.set_ylim(ymin, ymax)
        
        # Mark supremum
        supremum = ConnectionPatch(xyA=(data_all[index], cdf1[index]), xyB=(data_all[index], cdf2[index]), coordsA=ax1.transData, 
                           arrowstyle='<->', color='k')
        fig.add_artist(supremum)

        ax1.set_title('$D_n = sup_x |F_n(x)-F(x)| $', color='k', fontsize=14)
        
        # Add zoom lines
        con1 = ConnectionPatch(xyA=(xmin, ymin), coordsA=ax[0].transData, xyB=(xmin, ymin), coordsB=ax1.transData, alpha=0.3)
        con2 = ConnectionPatch(xyA=(xmax, ymax), coordsA=ax[0].transData, xyB=(xmax,ymax), coordsB=ax1.transData, alpha=0.3)

        sq1 = ConnectionPatch(xyA=(xmin, ymin), xyB=(xmax, ymin), coordsA=ax[0].transData, alpha=0.3)
        sq2 = ConnectionPatch(xyA=(xmin, ymax), xyB=(xmax, ymax), coordsA=ax[0].transData, alpha=0.3)
        sq3 = ConnectionPatch(xyA=(xmin, ymin), xyB=(xmin, ymax), coordsA=ax[0].transData, alpha=0.3)
        sq4 = ConnectionPatch(xyA=(xmax, ymin), xyB=(xmax, ymax), coordsA=ax[0].transData, alpha=0.3)

        fig.add_artist(con1)
        fig.add_artist(con2)
        fig.add_artist(sq1)
        fig.add_artist(sq2)
        fig.add_artist(sq3)
        fig.add_artist(sq4)
        
    # Perform ks test
    ks_stat, pval = stats.kstest(data1, data2, alternative='two-sided' )

    # Adding fit results to plot:
    d = {'KS stat':     ks_stat, 'Prob':     pval,}
    
    text = nice_string_output(d, extra_spacing=2, decimals=3)
    add_text_to_ax(0.02, 0.82, text, ax[0], fontsize=14)

    plt.show()
    
    return ks_stat, pval


def chi2_2samp(data1, data2, N_bins=50, label1='Data 1', label2='Data 2', xlabel=None, txtcoord=(0.02, 0.95)):
    
    """
    INPUT:
    data1 = 1d arraylike, with data sample 1
    data2 = 1d arraylike, with data sample 2
    
    OUTPUT:
    chi2 = chi2 value
    prob = chi2 probability of the two samples matching
    
    """
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12,6))
    
    # Plot tranformed data in the same range
    rang = ( min(*data1, *data2), max(*data1, *data2) )
    
    counts1, bin_edges1, _ = ax.hist(data1, bins=N_bins, range=rang, color='r', histtype='step', linewidth=2, label=label1)
    counts2, bin_edges2, _ = ax.hist(data2, bins=N_bins, range=rang, color='b', histtype='step', linewidth=2, label=label2)

    x1 = (bin_edges1[1:] + bin_edges1[:-1])/2
    sy1 = np.sqrt(counts1)

    x2 = (bin_edges2[1:] + bin_edges2[:-1])/2
    sy2 = np.sqrt(counts2)

    ax.errorbar(x1, counts1, yerr=sy1, fmt='.r',  ecolor='r', elinewidth=1, capsize=1, capthick=1, alpha=0.5)
    ax.errorbar(x2, counts2, yerr=sy2, fmt='.b',  ecolor='b', elinewidth=1, capsize=1, capthick=1, alpha=0.5)

    # Chi2
    chi2 = np.sum( (counts1 - counts2)**2  / ( sy1**2 + sy2**2 ) )
    ndof = len(counts1) #no degrees of freedom
    prob = stats.chi2.sf(chi2, ndof)

    # Adding fit results to plot:
    d = {'Chi2':     chi2,
         'ndof':     ndof,
         'Prob':     prob,
        }

    text = nice_string_output(d, extra_spacing=2, decimals=3)
    add_text_to_ax(*txtcoord, text, ax, fontsize=18)

    ax.legend()
    
    return chi2, prob


def find_C(fx_expr, xmin, xmax, all_sol = False):
    
    """
    INPUT:
    fx_expr = the expression for f(x) in a string, fx '1+x' or 'C*(1+x)'
    xmin = the lower limit for which x is defined (closed interval)
    xmax = the upper limit for whcih x is defined (closed interval)
    
        if one of the limmits is C, this should be passed as xmin = symbols('C')
        
    all_sol = Print all solutions, also the ones discarded
    
    OUTPUT:
    C = the value of C in a sympy representation to obtain a normalised function
    """
    
    # Load the sympy variables we usually use for Monte Carlo
    x, f, C = symbols("x, f, C")
    
    # Sympify expression
    fx_expr = sympify(fx_expr)
    
    # Find the integrated function F
    integral = integrate(fx_expr, (x, xmin, xmax))
    
    # Find all solutions of C so the integral equals 1
    solutions = solve(integral - 1, C)
    
    # Find the real, positive solution
    for sol in solutions:
        if sol.is_real and sol.is_positive: C_val = sol
    
    # Raise a warning is solutions are discarded
    
    if not all_sol and len(solutions) > 1:             
        print('Non-real and negative solutions were discarded. Set all_sol = True, to see all solutions')
    
    if all_sol and len(solutions) > 1:
        for expr in solutions: lprint(latex(Eq(symbols('C'), expr)))
        
    return C_val




def find_invF(fx_expr_no_C, C_val=None, xmin=-oo, all_sol = False ):
    
    """
    DISCLAIMER: This is not a pleasant function. We should really do this more 
    generally by creating a function generator and another function using that function.
    This will allow us to define the sympy symbols before referencing them in the arguments.
    The arguments would thus be keywordarguments which can be read as a dictionary. See the error 
    propagator for example. But im tired, so this is what you get.
    
    OBS If there are more real solutions to the inverse function (is it possible?) 
    They will be overwritten when invF_expr = sol
    
    INPUT:
    fx_expr_noC = the expression for f(x) in a string, leaving out a potential C factor: 
                  for example '5+x**2' even thought the right function is C*(5+x**2)
    
    C_val = the value of C, if theres is one in the function. If it is in the limits, ignore it.
    xmin = the lower limit where the function is defined
    all_sol = Print all solutions, also the ones discarded
    
    OUTPUT:
    inverse_function = python function representing the inverted function displayed in form f(u),
                       where the input u should be uniform numbers between 0 and 1.
    
    """
    
    # Gather expression if C is defined (I'm sorry)
    if C_val != None: fx_expr = str(C_val) + '*' + '(' + fx_expr_no_C + ')'
    if C_val == None: fx_expr = fx_expr_no_C
    
    # Load the sympy variables we usually use for Monte Carlo inverse function
    x, f, F, C, u = symbols("x, f, F, C, u")
    
    # Sympify expression for the function
    fx_expr = sympify(fx_expr)
    
    # Integrate the function from -inft to x, insert xmin if the function is not defined before
    F_expr = integrate( fx_expr, (x, xmin, x) )
    
    # Print the anti derivative
    lprint(latex(Eq(symbols('F(x)'), F_expr)))
    
    # Inverse function: set F(x) equal to u (uniform variable between 0,1) and isolate for x
    # Find all solutions to this
    solutions = solve(Eq(F_expr, u), x)

    # Find solution and display
    for sol in solutions:
        if 'I' not in str(sol): 
            invF_expr = sol   
            lprint(latex(Eq(symbols('F^{-1}(u)'), invF_expr)))
            
    # Raise a warning is solutions are discarded
    if not all_sol and len(solutions) > 1:             
        print('Non-real and negative solutions were discarded. Set all_sol = True, to see all solutions')
    
    if all_sol and len(solutions) > 1:
        for expr in solutions: lprint(latex(Eq(symbols('C'), expr)))
          
    # Lambdify the inverse function
    inverse_function = lambdify(u, invF_expr)
    
    return inverse_function


def accept_reject(func, N_points, xbound, ybound):
    
    """
    Notice that we look at the accept rejection as a binomial distribution 
    (either you accept (success) or reject (failure), where the efficiency 
    is the probability p.) The error on the efficiency comes from this approach,
    and is further described in http://www.pp.rhul.ac.uk/~cowan/atlas/efferr.ps
    
    INPUT:
    func = function the generated numbers should follow
    N_points = number of points to be generated
    xbound = (xmin, xmax) where xmin, xmax defines the x limits for the function
    ybound = (ymin, ymax) where ymin and ymax are the lower and upper bound for y
    
    OUTPUT:
    x_accept = 1d arraylike with all the accepted x-values
    [eff, eff_err] = tuple with efficiency and efficiency error
    [integral, integral_err] = tuple with integral value and itegral error
    
    """
    
    # Set counter of tries to zero
    N_try = 0 
    
    # Empty array to store the accepted values
    x_accept = np.zeros(N_points)
    
    # Generate N_points random numbers according to the distribution func
    for i in range(N_points):
        
        # Loop until we find an accepted value
        while True: 
            
            # Counts attempts
            N_try += 1 
            
            # Create x,y pair
            x_test, y_test = np.random.uniform(*xbound), np.random.uniform(*ybound)
            
            # Test that the point is within the distribution, if yes break loop
            if y_test < func(x_test): break
                
        x_accept[i] = x_test
        
    # Calculate efficiency
    eff = N_points / N_try  

    # Error on efficiency (binomial, see note above)
    eff_err = np.sqrt(eff * (1-eff) / N_try) 
            
    # Estimate integral, as the efficiency times the area of the bounding box
    integral = eff * np.diff(xbound)[0] * np.diff(ybound)[0]
    
    # Integral error
    integral_err = eff_err * np.diff(xbound)[0] * np.diff(ybound)[0]
    
    return x_accept, [eff, eff_err], [integral, integral_err]


def get_corr(x_data, y_data):
    """
    Calculates the linear correlation of two variables, given a data set with values for each variable
    
    INPUT:
    x_data = 1d arraylike with all x values of a data set where x is a variable
    y_data = 1d arraylike with all corresponding y values of a data set
    
    OUTPUT:
    rho_xy = the linear correlation
    """
    
    # Covariance matrix
    cov_mat = np.cov(x_data, y_data)
    
    # extract covariance and variance for each variable
    cov, var_x, var_y = cov_mat[0,1], cov_mat[0,0], cov_mat[1,1]
    
    # Find correlation cov / sigma_x * sigma_y. Remember sigma = sqrt(Var)
    rho_xy = cov / (np.sqrt(var_x)*np.sqrt(var_y))
    
    return rho_xy



def error_rates(species_A, species_B, cut, N_bins = 50, plot=True, alp_coord=(0.1, 0.52), bet_coord=(0.1, 0.52), labelA = 'Species A', labelB='Species B'):
    """
    INPUT:
    species_A : 1d arraylike containing data for one variable of species A. 
                Species A should be below the cut.
    species_B : 1d arraylike containing data for same variable of species B. 
                Species B should be above or equal to the cut.
    cut = float, defining the cut
    N_bins = number of bins to plot for each species
    plot = option to plot or not
    alp_coord = coordinates to place the axtext alpha on the form (x,y) where 0 < x,y < 1
    bet_coord = coordinates to place the axtext beta on the form (x,y) where 0 < x,y < 1
    
    OUTPUT:
    alp = type I error rate
    bet = type II error rate
    """

    # Separate data acoording to cut
    A_sup_cut = species_A[species_A >= cut]
    B_sub_cut = species_B[species_B < cut]
    
    # Calculate error rates
    alp = len(B_sub_cut) / len(species_B)
    bet = len(A_sup_cut) / len(species_A)

    # Create figure
    if plot:
        fig, ax = plt.subplots(ncols=2, figsize=(12,6))
        
        # Define range
        rang = ( np.min(np.concatenate((species_A, species_B))), np.max( np.concatenate((species_A, species_B))) )
        
        # Type I errors (alpha)---------------------------------------------------------------------------
        countA, _, _ = ax[0].hist(species_A, bins=N_bins, range=rang, histtype='step', color='red', linewidth=2, label=labelA)
        countB, _, _ = ax[0].hist(species_B, bins=N_bins, range=rang, histtype='step', color='blue', linewidth=2, label=labelB)

        ax[0].vlines(cut, 0, np.max((countA, countB))+10, color='k', linestyle='dashed', label='Cut=9$\mu$m')

        # Mark the area under
        ax[0].hist(B_sub_cut, bins=N_bins, range=rang, histtype='stepfilled', color='blue', alpha=0.5,linewidth=2, label='Beta')

        d1 = {'Alpha':     alp}  
        text1 = nice_string_output(d1, extra_spacing=2, decimals=3)
        add_text_to_ax(*alp_coord, text1, ax[0], fontsize=14, color='blue')

        ax[0].set_title('Type I Error', fontsize=14)
        ax[0].legend()
        
        # Type II errors (beta) ---------------------------------------------------------------------------
        ax[1].hist(species_A, bins=N_bins, range=rang, histtype='step', color='red', linewidth=2, label=labelA)
        ax[1].hist(species_B, bins=N_bins, range=rang, histtype='step', color='blue', linewidth=2, label=labelB)

        ax[1].vlines(cut, 0, np.max((countA, countB))+10, color='k', linestyle='dashed', label=f'Cut={cut}')

        # Mark the area under
        ax[1].hist(A_sup_cut, bins=N_bins, range=rang, histtype='stepfilled', color='red', alpha=0.5, linewidth=2, label='Alpha')

        d2 = {'Beta':     bet}  
        text2 = nice_string_output(d2, extra_spacing=2, decimals=3)
        add_text_to_ax(*bet_coord, text2, ax[1], fontsize=14, color='red')
    
        ax[1].set_title('Type II Error', fontsize=14)
        ax[1].legend()

        plt.show()
        
    return alp, bet


def gauss_hist(data, label=None, N_bins=None, ax=None):
    
    """
    Fits a Gaussian distribution to a histogram of the data

    INPUT:
    data = 1d array of all data points
    label = label to put on xaxis
    N_bins = number of bins, if not stated it will use Sturges formula
    ax = optional ax to plot on, if not stated will create a new figure
    
    OUTPUT:
    mean = (mean of data, error on mean from fit)
    sigma = (std of data, error on std from fit)
    Prob_gauss = probability of chi2 fit
    """
    
    # Create a figure if needed
    if ax is None: fig, ax = plt.subplots(figsize=(12,6))

    # Extract values from histogram and outline data
    if N_bins is None: N_bins = Sturges_bins(data)
    counts, bin_edges, _ = ax.hist(data, bins=N_bins, alpha=0.3, label='Histogram of data')
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
    binwidth = bin_edges[1] - bin_edges[0]

    # Poisson errors on the count in each bin
    s_counts = np.sqrt(counts)
    
    # We remove any bins, which don't have any counts in them:
    x = bin_centers[counts>0]
    y = counts[counts>0]
    sy = s_counts[counts>0]

    # Plot data with error
    ax.errorbar(x, y, yerr=sy, fmt='.k',  ecolor='k', elinewidth=1, capsize=1, capthick=1, label='Counts with Poisson errors')
    
    # Perform Gaussian fit ------------------------------------------------------------------------
    
    # Define Gauss function
    def func_gaussian(x, N, mu, sigma) :
        return N * binwidth * stats.norm.pdf(x, mu, sigma)

    # Perform fit by minimizing chi2
    chi2_gaussian = Chi2Regression(func_gaussian, x, y, sy)
    minuit_gaussian = Minuit(chi2_gaussian, pedantic=False, N=len(data), mu=np.mean(data), sigma=np.std(data, ddof=1))  
    minuit_gaussian.migrad(); 
    
    # Extract chi2 values
    chi2_gauss = minuit_gaussian.fval
    Ndof_gauss = len(x) - minuit_gaussian.narg
    Prob_gauss = stats.chi2.sf(chi2_gauss, Ndof_gauss)
    
    #Extract values from fit
    mean = [minuit_gaussian.values['mu'], minuit_gaussian.errors['mu']]
    sigma = [minuit_gaussian.values['sigma'], minuit_gaussian.errors['sigma']]

    # Plot Gaussian fit
    xaxis = np.linspace(min(data), max(data), 100)
    y_gauss = func_gaussian(xaxis, *minuit_gaussian.args)
    ax.plot(xaxis, y_gauss, linewidth=2, color='r')
    
    d = {'N': [minuit_gaussian.values['N'], minuit_gaussian.errors['N']],
          'mu': mean,
          'sigma': sigma,
          'Ndof': Ndof_gauss,
          'Chi2': chi2_gauss,
          'Prob': Prob_gauss}

    text = nice_string_output(d, extra_spacing=3, decimals=3)
    add_text_to_ax(0.02, 0.9, text, ax, fontsize=12, color='red')
    ax.text(0.02, 0.97, 'Gaussian Fit', weight='heavy', fontsize=15, transform=ax.transAxes, 
            verticalalignment='top', color='r')

    ax.legend()
    ax.set_xlabel(label, fontsize=14)
    
    return mean, sigma, Prob_gauss


def studentT_hist(data, label=None, N_bins=None, ax=None):
    
    """
    Fits a student t distribution to a histogram of the data.

    INPUT:
    data = 1d array of all data points
    label = label to put on xaxis
    N_bins = number of bins, if not stated it will use Sturges formula
    ax = optional ax to plot on, if not stated will create a new figure

    OUTPUT:
    loc = (mean of the data, error on mean from fit)
    scale = (corresponding to the standard deviation of the data, error on scale from fit)
    Prob_t = probability of chi2 fit
    """
    
    # Create a figure
    if ax is None: fig, ax = plt.subplots(figsize=(12,6))

    # Extract values from histogram and outline data
    if not N_bins: N_bins = Sturges_bins(data)
    counts, bin_edges, _ = ax.hist(data, bins=N_bins, alpha=0.3, label='Histogram of data')
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
    binwidth = bin_edges[1] - bin_edges[0]
    
    # Poisson errors on the count in each bin
    s_counts = np.sqrt(counts)
    
    # We remove any bins, which don't have any counts in them:
    x = bin_centers[counts>0]
    y = counts[counts>0]
    sy = s_counts[counts>0]

    # Plot data with error
    ax.errorbar(x, y, yerr=sy, fmt='.k',  ecolor='k', elinewidth=1, capsize=1, capthick=1, label='Counts with Poisson errors')
    
    # Perform Students T fit ------------------------------------------------------------------------
    
    # Defining the student t distribution (normalised)
    def func_t(x, N, df, loc, scale) :
        return N * binwidth * stats.t.pdf(x, df, loc, scale)

    # Perform fit by minimizing chi2
    chi2_t = Chi2Regression(func_t, x, y, sy)
    minuit_t = Minuit(chi2_t, pedantic=False, N=len(data), df=len(data)-4, loc=np.mean(data), scale=np.std(data))  
    minuit_t.migrad();     
    
    # Extract chi2 values
    chi2_t = minuit_t.fval
    Ndof_t = len(x) - minuit_t.narg
    Prob_t = stats.chi2.sf(chi2_t, Ndof_t)
    
    #Extract values from fit
    loc = [minuit_t.values['loc'], minuit_t.errors['loc']]
    scale = [minuit_t.values['scale'], minuit_t.errors['scale']]
    
    # Plot student t fit
    xaxis = np.linspace(min(data), max(data), 100)
    y_t = func_t(xaxis, *minuit_t.args)
    ax.plot(xaxis, y_t, linewidth=2, color='r')

    d = {'N': [minuit_t.values['N'], minuit_t.errors['N']],
         'df': [minuit_t.values['df'], minuit_t.errors['df']],
         'loc': loc,
         'scale': scale,
         'Ndof': Ndof_t,
         'Chi2': chi2_t,
         'Prob': Prob_t}

    text = nice_string_output(d, extra_spacing=3, decimals=3)
    add_text_to_ax(0.02, 0.9, text, ax, fontsize=12, color='r')
    ax.text(0.02, 0.97, 'Student T Fit', weight='heavy', fontsize=15,
            transform=ax.transAxes, verticalalignment='top', color='r')

    ax.legend()
    ax.set_xlabel(label, fontsize=14)
    
    return loc, scale, Prob_t





def double_gauss_hist(data, label=None, N_bins=None, ax=None, sub_pdfs = True, p0=()):

    """
    Fits a double Gaussian distribution to a histogram of the data.

    INPUT:
    data = 1d array of all data points
    label = label to put on xaxis
    N_bins = number of bins, if not stated it will use Sturges formula
    ax = optional ax to plot on, if not stated will create a new figure
    sub_pdfs = option to mark the subset of gaussians
    p0 = initial guesses on the form (frac, mu1, sigma1, mu2, sigma2)
    
    OUTPUT:
    mean1 = (mean of first Gaussian fit, error on mean from fit)
    sigma1 = (std of first Gaussian fit, error on std from fit)
    mean2 = (mean of second Gaussian fit, error on mean from fit)
    sigma2 = (std of second Gaussian fit, error on std from fit)
    Prob_double = probability of chi2 fit
    """

    # Create a figure
    if ax is None: fig, ax = plt.subplots(figsize=(12,6))

    # Extract values from histogram and outline data
    if not N_bins: N_bins = Sturges_bins(data)
    counts, bin_edges, _ = ax.hist(data, bins=N_bins, alpha=0.3, label='Histogram of data')
    bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
    binwidth = bin_edges[1] - bin_edges[0]
    
    # Poisson errors on the count in each bin
    s_counts = np.sqrt(counts)
    
    # We remove any bins, which don't have any counts in them:
    x = bin_centers[counts>0]
    y = counts[counts>0]
    sy = s_counts[counts>0]

    # Plot data with error
    ax.errorbar(x, y, yerr=sy, fmt='.k',  ecolor='k', elinewidth=1, capsize=1, capthick=1, label='Counts with Poisson errors')
    
    # Perform double Gaussian fit ------------------------------------------------------------------------

    # Define function for double gauss (normalised)
    def double_gauss(x, N, frac, mu1, sigma1, mu2, sigma2):
        # Notice we introduce a common parameter frac that determines the proportion of the two functions
        return N * binwidth * ( frac * stats.norm.pdf(x, mu1, sigma1) + (1-frac)*stats.norm.pdf(x, mu2, sigma2) )


    # Perform fit by minimizing chi2
    chi2_double = Chi2Regression(double_gauss, x, y, sy)

    if p0 is not None:
        minuit_double = Minuit(chi2_double, pedantic=False, N=len(data), frac = p0[0], mu1=p0[1], sigma1=p0[2], mu2 = p0[3], sigma2=p0[4] ) 
    else:
        minuit_double = Minuit(chi2_double, pedantic=False, N=len(data), frac = 0.5, mu1=np.mean(data), mu2 = np.mean(data),
                                sigma1=np.std(data, ddof=1), sigma2=np.std(data, ddof=1))  
    minuit_double.migrad(); 
    
    # Extract chi2 values
    chi2_double = minuit_double.fval
    Ndof_double = len(x) - minuit_double.narg
    Prob_double = stats.chi2.sf(chi2_double, Ndof_double)

    # Extract values from fit
    mean1 = [minuit_double.values['mu1'], minuit_double.errors['mu1']]
    sigma1 = [minuit_double.values['sigma1'], minuit_double.errors['sigma1']]
    mean2 = [minuit_double.values['mu2'], minuit_double.errors['mu2']]
    sigma2 = [minuit_double.values['sigma2'], minuit_double.errors['sigma2']]
    frac = [minuit_double.values['frac'], minuit_double.errors['frac']]

    # Plot Double Gaussian fit
    xaxis = np.linspace(min(data), max(data), 100)
    y_double = double_gauss(xaxis, *minuit_double.args)
    ax.plot(xaxis, y_double, linewidth=2, color='r')

    # Plot the sub distributions
    if sub_pdfs:
        
        # Define a regular gaussian
        def gauss(x, N, frac, mu, sigma):
            return N * binwidth * frac * stats.norm.pdf(x, mu, sigma)

        ax.plot(xaxis, gauss(xaxis, len(data), frac[0], mean1[0], sigma1[0]), linestyle='dashed', color='r', label='Gauss 1' )
        ax.plot(xaxis, gauss(xaxis, len(data), 1-frac[0], mean2[0], sigma2[0]), linestyle='dashed', color='r', label='Gauss 2' )

    d = {'N': [minuit_double.values['N'], minuit_double.errors['N']],
         'frac': frac,
         'mu1': mean1,
         'sigma1': sigma1,
         'mu2': mean2,
         'sigma2': sigma2,
         'Ndof': Ndof_double,
         'Chi2': chi2_double,
         'Prob': Prob_double}

    text = nice_string_output(d, extra_spacing=3, decimals=3)
    add_text_to_ax(0.02, 0.9, text, ax, fontsize=12, color='red')
    ax.text(0.02, 0.97, 'Double Gaussian Fit', weight='heavy', fontsize=15, transform=ax.transAxes, 
                    verticalalignment='top', color='r')

    ax.legend()
    ax.set_xlabel(label, fontsize=14)

    return mean1, sigma1, mean2, sigma2, Prob_double

