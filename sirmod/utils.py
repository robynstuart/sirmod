'''
Utilities
'''

# Imports
import numpy as np
import pandas as pd
import sciris as sc
from scipy.stats import truncnorm, norm


def sample(dist=None, par1=None, par2=None, size=None, lower_clip=None, upper_clip=None, **kwargs):
    '''
    Draw a sample from the distribution specified by the input. The available
    distributions are:

    - 'trunc_norm'    : truncated normal distribution with

    Args:
        dist (str):     the distribution to sample from
        par1 (float):   the "main" distribution parameter (e.g. mean)
        par2 (float):   the "secondary" distribution parameter (e.g. std)
        size (int):     the number of samples (default=1)
        lower (float):  the lower bound (used when truncating)
        upper (float):  the upper bound (used when truncating)
        kwargs (dict): passed to individual sampling functions

    Returns:
        A length "size" array of samples

    **Examples**::
        cv.sample(dist='trunc_norm', par1=5, par2=3, lower=2, upper=6) # returns a truncated normal distributed set of values with mean 5 and std 3

    Notes:
        Scipy's truncated normal distributions interprets the truncation limits with reference to a standard
        normal distribution (see:https://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.stats.truncnorm.html),
         but this function assumes the user wants to specify them with reference to the normal distribution being used.
    '''

    # Some of these have aliases, but these are the "official" names
    choices = [
        'uniform',
        'trunc_norm',
    ]

    # Ensure it's an integer
    if size is not None:
        size = int(size)

    # Compute distribution parameters and draw samples
    if dist in ['unif', 'uniform']:
        samples = np.random.uniform(low=par1, high=par2, size=size, **kwargs)
    elif dist in ['norm', 'normal']: samples = np.random.normal(loc=par1, scale=par2, size=size, **kwargs)
    elif dist == 'normal_int': samples = np.round(
        np.abs(np.random.normal(loc=par1, scale=par2, size=size, **kwargs)))
    elif dist == 'poisson': samples = n_poisson(rate=par1, n=size, **kwargs)  # Use Numba version below for speed
    elif dist == 'neg_binomial': samples = n_neg_binomial(rate=par1, dispersion=par2, n=size,
                                                          **kwargs)  # Use custom version below
    elif dist in ['lognorm', 'lognormal', 'lognorm_int', 'lognormal_int']:
        if par1 > 0:
            mean = np.log(
                par1 ** 2 / np.sqrt(par2 ** 2 + par1 ** 2))  # Computes the mean of the underlying normal distribution
            sigma = np.sqrt(np.log(par2 ** 2 / par1 ** 2 + 1))  # Computes sigma for the underlying normal distribution
            samples = np.random.lognormal(mean=mean, sigma=sigma, size=size, **kwargs)
        else:
            samples = np.zeros(size)
        if '_int' in dist:
            samples = np.round(samples)
    elif dist in ['trunc_norm', 'truncnorm', 'truncated_normal', 'truncatednormal', 'tnorm']:
        if lower_clip is None: lower_clip = 0 # Arbitrary limits
        if upper_clip is None: upper_clip = 10  # Arbitrary limits
        a, b = (lower_clip - par1) / par2, (upper_clip - par1) / par2
        samples = truncnorm.rvs(a, b, loc=par1, scale=par2, size=size)
    else:
        errormsg = f'The selected distribution "{dist}" is not implemented; choices are: {sc.newlinejoin(choices)}'
        raise NotImplementedError(errormsg)

    return samples


def get_q(x, dist, par1, par2, lower_clip=None, upper_clip=None):
    ''' Function to get quantiles of a distribution'''
    if dist=='trunc_norm':
        a, b = (lower_clip - par1) / par2, (upper_clip - par1) / par2
        if x<=par1:
            q = truncnorm.cdf(x, a, b, loc=par1, scale=par2)
        else:
            q = 1-truncnorm.cdf(x, a, b, loc=par1, scale=par2)
    elif dist in ['norm', 'normal']:
        if x<=par1:
            q = norm.cdf(x, loc=par1, scale=par2)
        else:
            q = 1-norm.cdf(x, loc=par1, scale=par2)
    else:
        errormsg = f'The selected distribution "{dist}" is not implemented.'
        raise NotImplementedError(errormsg)

    return q


def extractdata(xdata, ydata):
    ''' Return the x and y data values for non-nan y data '''
    nonnanx = np.array(xdata)[~np.isnan(np.array(ydata))]
    nonnany = np.array(ydata)[~np.isnan(np.array(ydata))]
    return nonnanx, nonnany


class Result(object):
    ''' Class to hold overall and by-population results '''

    def __init__(self, name=None, ispercentage=False, pops=None, tot=None, datapops=None, datatot=None, estimate=False,
                 defaultplot=None):
        self.name = name  # Name of this parameter
        self.ispercentage = ispercentage  # Whether or not the result is a percentage
        self.pops = pops  # The model result by population, if available
        self.tot = tot  # The model result total, if available
        self.datapops = datapops  # The input data by population, if available
        self.datatot = datatot  # The input data total, if available
        self.estimate = estimate  # Whether or not the "data" is actually a model-based output
        self.defaultplot = defaultplot if defaultplot is not None else 'stacked'

    def __repr__(self):
        ''' Print out useful information when called '''
        output = sc.prepr(self)
        return output


def gettvecdt(tvec=None, dt=None, justdt=False):
    '''
    Function to encapsulate the logic of returning sensible tvec and dt based on flexible input.
    If tvec and dt are both supplied, do nothing.
    Will always work if tvec is not None, but will use default value for dt if dt==None and len(tvec)==1.
    Usage:
        tvec,dt = gettvecdt(tvec, dt)
    '''
    defaultdt = 0.2
    if tvec is None:
        if justdt:
            return defaultdt  # If it's a constant, maybe don' need a time vector, and just return dt
        else:
            raise Exception('No time vector supplied, and unable to figure it out')  # Usual case, crash
    elif sc.isnumber(tvec):
        tvec = np.array([tvec])  # Convert to 1-element array
    elif np.shape(tvec):  # Make sure it has a length -- if so, overwrite dt
        if len(tvec) >= 2:
            dt = tvec[1] - tvec[0]  # Even if dt supplied, recalculate it from the time vector
        else:
            dt = dt  # Use input
    else:
        raise Exception('Could not understand tvec of type "%s"' % type(tvec))
    if dt is None: dt = defaultdt  # Or give up and use default
    return tvec, dt


####################################################################################
# Distance functions
####################################################################################
def prep_distance(model, fitto=None):
    ''' Pre-process model to get out model results and data for calculating distance'''

    # Handle inputs
    results = sc.dcp(model.results)
    if fitto is None:
        fitto = ['pop_size', 'prev']

    modely, datay = np.array([]), np.array([])

    for key in fitto:
        thisres = results[key]
        tmpdata = thisres.datapops
        modelrows = thisres.pops

        # First extract datay
        if tmpdata is not None:
            datarows = tmpdata[0]  # Pull out best data
            nrows = len(datarows)
            for row in range(nrows):  # Loop over each available row (pops if it's a by-pop result, or single row if it's a total result)
                datarow = datarows[row]
                thisdatay = datarow[~np.isnan(datarow)]
                datay = np.append(datay, thisdatay)

                # Now extract modely, if available
                if modelrows is not None:
                    if len(modelrows.shape) > 1:    modelrow = modelrows[row]
                    else:                           modelrow = modelrows
                    thismodely = modelrow[~np.isnan(datarow)]
                    modely = np.append(modely, thismodely)

    return modely, datay


def distance(model=None, data=None, method='wape', eps=1e-3):
    ''' Evaluate how well the model fits to the data '''

    if method == 'wape':
        mismatch = np.sum(abs(model["Y"] - data["Y"]) / np.mean(data["Y"] + eps))
    elif method == 'mape':
        mismatch = np.sum(abs(model["Y"] - data["Y"]) / (data["Y"] + eps))
    elif method == 'mad':
        mismatch = np.sum(abs(model["Y"] - data["Y"]))
    elif method == 'mse':
        mismatch = np.sum((model["Y"] - data["Y"]) ** 2)
    else:
        errormsg = f'"method" not known; you entered {method}, but must be one of:\n'
        errormsg += '"wape" = weighted absolute percentage error (default)\n'
        errormsg += '"mape" = mean absolute percentage error\n'
        errormsg += '"mad"  = mean absolute difference\n'
        errormsg += '"mse"  = mean squared error'
        raise Exception(errormsg)

    return mismatch
