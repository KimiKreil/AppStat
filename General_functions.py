#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:22:28 2020

@author: Kimi Kreilgaard & Jonathan Steensgaard

"""

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
    
    OUTPUT: Probability of r
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
    
    # Find the standard deviation
    sigma = np.std(data, ddof=1)

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
    INPUT:
    data = empirical data to test against a fit
    fit_function = Scipy stats function that is fitted with in a string, fx 'norm'
    args = arguments to give fit function in default order
    
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
    cdf_axis = len(cut3_new_data) * cdf(xaxis, *args)
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
        index = np.argmax(resi)
        
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
    ks_stat, pval = stats.kstest(cut3_new_data, fit_function, args=args )

    # Adding fit results to plot:
    d = {'KS stat':     ks_stat, 'Prob':     pval,}
    
    text = nice_string_output(d, extra_spacing=2, decimals=3)
    add_text_to_ax(0.02, 0.82, text, ax[0], fontsize=14)

    plt.show()
    
    return ks_stat, pval




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








