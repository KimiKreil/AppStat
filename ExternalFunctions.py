# coding: utf-8

#Author: Christian Michelsen, NBI, 2018
#Autor1: Troels Petersen, NBI, 2020

#Edited by Jonathan Baltazar, wmh858, 2020 

import numpy as np                                    
import matplotlib.pyplot as plt

from IPython.core import display
from IPython.core.display import Latex

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


def format_value(value, decimals):
    """ 
    Checks the type of a variable and formats it accordingly.
    Floats has 'decimals' number of decimals.
    """
    
    if isinstance(value, (float, np.float)):
        return f'{value:.{decimals}f}'
    elif isinstance(value, (int, np.integer)):
        return f'{value:d}'
    else:
        return f'{value}'


def values_to_string(values, decimals):
    """ 
    Loops over all elements of 'values' and returns list of strings
    with proper formating according to the function 'format_value'. 
    """
    
    res = []
    for value in values:
        if isinstance(value, list):
            tmp = [format_value(val, decimals) for val in value]
            res.append(f'{tmp[0]} +/- {tmp[1]}')
        else:
            res.append(format_value(value, decimals))
    return res


def len_of_longest_string(s):
    """ Returns the length of the longest string in a list of strings """
    return len(max(s, key=len))


def nice_string_output(d, extra_spacing=5, decimals=5):
    """ 
    Takes a dictionary d consisting of names and values to be properly formatted.
    Makes sure that the distance between the names and the values in the printed
    output has a minimum distance of 'extra_spacing'. One can change the number
    of decimals using the 'decimals' keyword.  
    """
    
    names = d.keys()
    max_names = len_of_longest_string(names)
    
    values = values_to_string(d.values(), decimals=decimals)
    max_values = len_of_longest_string(values)
    
    string = ""
    for name, value in zip(names, values):
        spacing = extra_spacing + max_values + max_names - len(name) - 1 
        string += "{name:s} {value:>{spacing}} \n".format(name=name, value=value, spacing=spacing)
    return string[:-2]


def add_text_to_ax(x_coord, y_coord, string, ax, fontsize=12, color='k'):
    """ Shortcut to add text to an ax with proper font. Relative coords."""
    ax.text(x_coord, y_coord, string, family='monospace', fontsize=fontsize,
            transform=ax.transAxes, verticalalignment='top', color=color)
    return None




# =============================================================================
#  Probfit replacement
# =============================================================================

from iminuit.util import make_func_code
from iminuit import describe #, Minuit,

def set_var_if_None(var, x):
    if var is not None:
        return np.array(var)
    else: 
        return np.ones_like(x)
    
def compute_f(f, x, *par):
    
    try:
        return f(x, *par)
    except ValueError:
        return np.array([f(xi, *par) for xi in x])


class Chi2Regression:  # override the class with a better one
        
    def __init__(self, f, x, y, sy=None, weights=None, bound=None):
        
        if bound is not None:
            x = np.array(x)
            y = np.array(y)
            sy = np.array(sy)
            mask = (x >= bound[0]) & (x <= bound[1])
            x  = x[mask]
            y  = y[mask]
            sy = sy[mask]

        self.f = f  # model predicts y for given x
        self.x = np.array(x)
        self.y = np.array(y)
        
        self.sy = set_var_if_None(sy, self.x)
        self.weights = set_var_if_None(weights, self.x)
        self.func_code = make_func_code(describe(self.f)[1:])

    def __call__(self, *par):  # par are a variable number of model parameters
        
        # compute the function value
        f = compute_f(self.f, self.x, *par)
        
        # compute the chi2-value
        chi2 = np.sum(self.weights*(self.y - f)**2/self.sy**2)
        
        return chi2






def simpson38(f, edges, bw, *arg):
    
    yedges = f(edges, *arg)
    left38 = f((2.*edges[1:]+edges[:-1]) / 3., *arg)
    right38 = f((edges[1:]+2.*edges[:-1]) / 3., *arg)
    
    return bw / 8.*( np.sum(yedges)*2.+np.sum(left38+right38)*3. - (yedges[0]+yedges[-1]) ) #simpson3/8


def integrate1d(f, bound, nint, *arg):
    """
    compute 1d integral
    """
    edges = np.linspace(bound[0], bound[1], nint+1)
    bw = edges[1] - edges[0]
    
    return simpson38(f, edges, bw, *arg)



class UnbinnedLH:  # override the class with a better one
    
    def __init__(self, f, data, weights=None, bound=None, badvalue=-100000, extended=False, extended_bound=None, extended_nint=100):
        
        if bound is not None:
            data = np.array(data)
            mask = (data >= bound[0]) & (data <= bound[1])
            data = data[mask]
            if (weights is not None) :
                weights = weights[mask]

        self.f = f  # model predicts PDF for given x
        self.data = np.array(data)
        self.weights = set_var_if_None(weights, self.data)
        self.bad_value = badvalue
        
        self.extended = extended
        self.extended_bound = extended_bound
        self.extended_nint = extended_nint
        if extended and extended_bound is None:
            self.extended_bound = (np.min(data), np.max(data))

        
        self.func_code = make_func_code(describe(self.f)[1:])

    def __call__(self, *par):  # par are a variable number of model parameters
        
        logf = np.zeros_like(self.data)
        
        # compute the function value
        f = compute_f(self.f, self.data, *par)
    
        # find where the PDF is 0 or negative (unphysical)        
        mask_f_positive = (f>0)

        # calculate the log of f everyhere where f is positive
        logf[mask_f_positive] = np.log(f[mask_f_positive]) * self.weights[mask_f_positive] 
        
        # set everywhere else to badvalue
        logf[~mask_f_positive] = self.bad_value
        
        # compute the sum of the log values: the LLH
        llh = -np.sum(logf)
        
        if self.extended:
            extended_term = integrate1d(self.f, self.extended_bound, self.extended_nint, *par)
            llh += extended_term
        
        return llh
    
    def default_errordef(self):
        return 0.5


class BinnedLH:  # override the class with a better one
    
    def __init__(self, f, data, bins=40, weights=None, weighterrors=None, bound=None, badvalue=1000000, extended=False, use_w2=False, nint_subdiv=1):
        
        if bound is not None:
            data = np.array(data)
            mask = (data >= bound[0]) & (data <= bound[1])
            data = data[mask]
            if (weights is not None) :
                weights = weights[mask]
            if (weighterrors is not None) :
                weighterrors = weighterrors[mask]

        self.weights = set_var_if_None(weights, data)

        self.f = f
        self.use_w2 = use_w2
        self.extended = extended

        if bound is None: 
            bound = (np.min(data), np.max(data))

        self.mymin, self.mymax = bound

        h, self.edges = np.histogram(data, bins, range=bound, weights=weights)
        
        self.bins = bins
        self.h = h
        self.N = np.sum(self.h)

        if weights is not None:
            if weighterrors is None:
                self.w2, _ = np.histogram(data, bins, range=bound, weights=weights**2)
            else:
                self.w2, _ = np.histogram(data, bins, range=bound, weights=weighterrors**2)
        else:
            self.w2, _ = np.histogram(data, bins, range=bound, weights=None)


        
        self.badvalue = badvalue
        self.nint_subdiv = nint_subdiv
        
        
        self.func_code = make_func_code(describe(self.f)[1:])
        self.ndof = np.sum(self.h > 0) - (self.func_code.co_argcount - 1)
        

    def __call__(self, *par):  # par are a variable number of model parameters

        # ret = compute_bin_lh_f(self.f, self.edges, self.h, self.w2, self.extended, self.use_w2, self.badvalue, *par)
        ret = compute_bin_lh_f2(self.f, self.edges, self.h, self.w2, self.extended, self.use_w2, self.nint_subdiv, *par)
        
        return ret


    def default_errordef(self):
        return 0.5




import warnings


def xlogyx(x, y):
    
    #compute x*log(y/x) to a good precision especially when y~x
    
    if x<1e-100:
        warnings.warn('x is really small return 0')
        return 0.
    
    if x<y:
        return x*np.log1p( (y-x) / x )
    else:
        return -x*np.log1p( (x-y) / y )


#compute w*log(y/x) where w < x and goes to zero faster than x
def wlogyx(w, y, x):
    if x<1e-100:
        warnings.warn('x is really small return 0')
        return 0.
    if x<y:
        return w*np.log1p( (y-x) / x )
    else:
        return -w*np.log1p( (x-y) / y )


def compute_bin_lh_f2(f, edges, h, w2, extended, use_sumw2, nint_subdiv, *par):
    
    N = np.sum(h)
    n = len(edges)

    ret = 0.
    
    for i in range(n-1):
        th = h[i]
        tm = integrate1d(f, (edges[i], edges[i+1]), nint_subdiv, *par)
        
        if not extended:
            if not use_sumw2:
                ret -= xlogyx(th, tm*N) + (th-tm*N)

            else:
                if w2[i]<1e-200: 
                    continue
                tw = w2[i]
                factor = th/tw
                ret -= factor*(wlogyx(th,tm*N,th)+(th-tm*N))
        else:
            if not use_sumw2:
                ret -= xlogyx(th,tm)+(th-tm)
            else:
                if w2[i]<1e-200: 
                    continue
                tw = w2[i]
                factor = th/tw
                ret -= factor*(wlogyx(th,tm,th)+(th-tm))

    return ret



def compute_bin_lh_f(f, edges, h, w2, extended, use_sumw2, badvalue, *par):
    
    mask_positive = (h>0)
    
    N = np.sum(h)
    midpoints = (edges[:-1] + edges[1:]) / 2
    b = np.diff(edges)
    
    midpoints_pos = midpoints[mask_positive]
    b_pos = b[mask_positive]
    h_pos = h[mask_positive]
    
    if use_sumw2:
        warnings.warn('use_sumw2 = True: is not yet implemented, assume False ')
        s = np.ones_like(midpoints_pos)
        pass
    else: 
        s = np.ones_like(midpoints_pos)

    
    E_pos = f(midpoints_pos, *par) * b_pos
    if not extended:
        E_pos = E_pos * N
        
    E_pos[E_pos<0] = badvalue
    
    ans = -np.sum( s*( h_pos*np.log( E_pos/h_pos ) + (h_pos-E_pos) ) )

    return ans




def profile_x(x, y, bins=(50, 50), xyrange=[(0, 50), (-1,1)]):
    
    H, xedges, yedges = np.histogram2d(x, y, bins=bins, range=xyrange)
    x_center = 0.5*(xedges[1:] + xedges[:-1])
    y_center = 0.5*(yedges[1:] + yedges[:-1])
    
    wsums = H.sum(1)
    
    mask = wsums > 0
    
    mean = (H*y_center).sum(1)[mask] / wsums[mask]
    mean_squared = (H*y_center**2).sum(1)[mask] / wsums[mask]
    std = np.sqrt( mean_squared - mean**2 ) / np.sqrt(wsums[mask]) 

    return x_center[mask], mean, std

def FDA(spec_A, spec_B, labels=None):

    mu_A = np.mean(spec_A, 0)
    mu_B = np.mean(spec_B, 0)

    cov_A = np.cov(spec_A.T)
    cov_B = np.cov(spec_B.T)
    cov_sum = cov_A + cov_B

    # inverts cov_sum
    cov_sum_inv = np.linalg.inv(cov_sum)

    wf = np.dot(cov_sum_inv.T, (mu_A - mu_B))

    fisher_data_A = np.dot(wf, spec_A.T)
    fisher_data_B = np.dot(wf, spec_B.T)

    return fisher_data_A, fisher_data_B, wf

def calc_ROC(hist1, hist2):

    # hist1 is signal, hist2 is background

    # first we extract the entries (y values) and the edges of the histograms
    y_sig, x_sig_edges, __ = hist1
    y_bkg, x_bkg_edges, __ = hist2

    # Check that the two histograms have the same x edges:
    if np.array_equal(x_sig_edges, x_bkg_edges):

        # extract the center positions (x values) of the bins (doesn't matter if we use signal or background because they are equal)
        x_centers = 0.5*(x_sig_edges[1:] + x_sig_edges[:-1])

        # calculate the integral (sum) of the signal and background
        integral_sig = y_sig.sum()
        integral_bkg = y_bkg.sum()

        # initialize empty arrays for the True Positive Rate (TPR) and the False Positive Rate (FPR).
        TPR = np.zeros_like(y_sig) # True positive rate (sensitivity)
        FPR = np.zeros_like(y_sig) # False positive rate ()

        # loop over all bins (x_centers) of the histograms and calculate TN, FP, FN, TP, FPR, and TPR for each bin
        for i, x in enumerate(x_centers):

            # the cut mask
            cut = (x_centers < x)

            # true positive
            TP = np.sum(y_sig[~cut]) / integral_sig    # True positives
            FN = np.sum(y_sig[cut]) / integral_sig     # False negatives
            TPR[i] = TP / float(TP + FN)                    # True positive rate

            # true negative
            TN = np.sum(y_bkg[cut]) / integral_bkg      # True negatives (background)
            FP = np.sum(y_bkg[~cut]) / integral_bkg     # False positives
            FPR[i] = FP / float(FP + TN)                     # False positive rate

        return FPR, TPR

    else:
        AssertionError("Signal and Background histograms have different bins and ranges")


def plot_ROC(FPR, TPR, labels=None, colors=None, figsize=(8,6), ax=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)


    for i in np.arange(len(FPR)):

        kwargs = {}
        if colors is not None:
            kwargs['color'] = colors[i]
        if labels is not None:
            kwargs['label'] = labels[i]

        ax.plot(FPR[i], TPR[i], **kwargs)


    ax.plot([0,1], [0,1], 'k--')
    if labels is not None:
        ax.legend()
    ax.set(xlabel='False Positive Rate (Background Efficiency)', ylabel='True Positive Rate (Signal Efficiency)', 
           xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    return ax

# Defining a function that returns random correlated numbers: 
def Corr(mu1, sig1, mu2, sig2, rho12): 
    
    # 
    r = np.random
    r.seed(42)
    
    theta = 0.5 * np.arctan( 2.0 * rho12 * sig1 * sig2 / ( sig1**2 - sig2**2 ) )
    sigu = np.sqrt( np.abs( ((sig1*np.cos(theta))**2 - (sig2*np.sin(theta))**2 ) / ( np.cos(theta)**2 - np.sin(theta)**2) ) )
    sigv = np.sqrt( np.abs( ((sig2*np.cos(theta))**2 - (sig1*np.sin(theta))**2 ) / ( np.cos(theta)**2 - np.sin(theta)**2) ) )

    u = r.normal(0.0, sigu)
    v = r.normal(0.0, sigv)

    x = mu1 + np.cos(theta)*u - np.sin(theta)*v
    y = mu2 + np.sin(theta)*u + np.cos(theta)*v

    return x, y

def calc_separation(x, y):
    d = np.abs((np.mean(x) - np.mean(y))) / np.sqrt(np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2)
    return d
