'''
Function to create and manipulate parameters
'''

# Imports
import numpy as np
import pandas as pd
import sciris as sc
import sciris as sc
from . import utils as smu
from xlrd import open_workbook, colname


# # Helper functions
# def make_pars(format='dists', size=None):
#     '''
#     Create scale factors over each of the parameters in the model
#     Args:
#         format (str): 'sample' returns a sample of size "size" of each parameter, 'dists' returns the distributions, 'point_estimate' returns the mean (best) estimates
#     '''
#     pardists = sc.objdict()
#     pardists.init_prev      = dict(dist='trunc_norm', par1=1, lower_95=0.5, upper_95=2, lower_clip=0, upper_clip=5)  # Initial prevalence
#     pardists.death_other    = dict(dist='trunc_norm', par1=1, lower_95=0.5, upper_95=2, lower_clip=0, upper_clip=5)  # Annual probability of death from other causes
#     pardists.death          = dict(dist='trunc_norm', par1=1, lower_95=0.5, upper_95=2, lower_clip=0, upper_clip=5)  # Annual probability of death from modelled disease
#     pardists.birth          = dict(dist='trunc_norm', par1=1, lower_95=0.5, upper_95=2, lower_clip=0, upper_clip=5)  # Birth rate
#     pardists.foi            = dict(dist='trunc_norm', par1=1, lower_95=0.5, upper_95=2, lower_clip=0, upper_clip=5)  # Force of infection
#
#     # Derive stdev from 95% confidence intervals
#     for k,v in pardists.items():
#         if v.get('par2') is None:
#             v['par2'] = (v['upper_95']-v['lower_95'])/1.96
#
#     # Return parameters in the desired format
#     if format=='point_estimate':
#         pars = {k:v['par1'] for k,v in pardists.items()}
#     elif format=='sample':
#         pars = {}
#         for k,v in pardists.items():
#             v.update(size=size)
#             pars = {k:hmu.sample(**v) for k,v in pardists.items()}
#     elif format=='dists':
#         pars = pardists
#     else:
#         errormsg = f'The selected format "{format}" is not implemented.'
#         raise NotImplementedError(errormsg)
#
#     return pars
#

def load_par_specs(file_name='model-inputs.xlsx', folder='../sirmod', sheet_name='Model parameters'):
    '''  Function to parse the parameter definitions from the spreadsheet and return a structure that can be used to generate the parameters '''
    sheet_name = 'Model parameters'
    full_path = sc.makefilepath(filename=file_name, folder=folder)
    workbook = open_workbook(full_path)
    sheet = workbook.sheet_by_name(sheet_name)

    rawpars = []
    for rownum in range(sheet.nrows - 1):
        rawpars.append({})
        for colnum in range(sheet.ncols):
            attr = sheet.cell_value(0, colnum)
            rawpars[rownum][attr] = sheet.cell_value(rownum + 1, colnum) if sheet.cell_value(rownum + 1, colnum) != 'None' else None
            if sheet.cell_value(0, colnum) in ['limits']:
                rawpars[rownum][attr] = eval(sheet.cell_value(rownum + 1, colnum))  # Turn into actual values
    return rawpars


def load_trans_specs(file_name='model-inputs.xlsx', folder='../sirmod', sheet_name='Transitions', npops=None):
    ''' Function to load the allowable transitions from the spreadsheet '''
    if npops is None: npops = 1  # Use just one population if not told otherwise
    fullpath = sc.makefilepath(filename=file_name, folder=folder)
    workbook = open_workbook(fullpath)
    sheet = workbook.sheet_by_name(sheet_name)

    if sheet.nrows != sheet.ncols:
        errormsg = f'Transition matrix should have the same number of rows and columns ({sheet.nrows} vs. {sheet.ncols})'
        raise Exception(errormsg)
    nstates = sheet.nrows - 1  # First row is header

    fromto = []
    transmatrix = np.zeros((nstates, nstates, npops))
    for rownum in range(nstates):  # Loop over each health state: the from state
        fromto.append([])  # Append two lists: the to state and the probability
        for colnum in range(nstates):  # ...and again
            if sheet.cell_value(rownum + 1, colnum + 1):
                fromto[rownum].append(colnum)  # Append the to states
                transmatrix[rownum, colnum, :] = np.ones(npops)  # Append the probabilities

    return fromto, transmatrix


#################################################################################################################################
### Define the other classes
#################################################################################################################################

class Par(object):
    '''
    The base class for epidemiological model parameters.
    '''

    def __init__(self, short=None, name=None, limits=(0., 1.), by=None, manual='', fromdata=None, **defaultargs):
        ''' To initialize with a prior, prior should be a dict with keys 'dist' and 'pars' '''
        self.short      = short     # The short name
        self.name       = name      # The full name
        self.limits     = limits    # Parameter limits
        self.by         = by        # Whether it's by population, partnership, or total
        self.fromdata   = fromdata  # Whether or not the parameter is made from data

    def __repr__(self):
        ''' Print out useful information when called'''
        output = sc.prepr(self)
        return output



class Metapar(Par):
    ''' The definition of a single metaparameter '''

    def __init__(self, y=None, prior=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        self.y = y
        self.ysample = None
        if not self.limits or self.limits is None:
            self.limits = (0.05, 50)  # Arbitrary for metapars
        if isinstance(prior, dict):
            self.prior = prior
        elif prior is None:
            self.prior = sc.odict()
            for key in self.keys():
                self.prior[key] = dict(dist='trunc_norm', par1=1, par2=1, lower_clip=self.limits[0], upper_clip=self.limits[1])
        else:
            errormsg = f'Prior for metaparameters must be an odict, not {type(prior)}'
            raise Exception(errormsg)

    def keys(self):
        ''' Return the valid keys for using with this parameter '''
        return self.y.keys()

    def sample(self, randseed=None, size=1):
        ''' Recalculate ysample '''
        self.ysample = sc.odict()
        for key in self.keys():
            args = sc.mergedicts(self.prior[key],dict(randseed=randseed, size=size))
            self.ysample[key] = smu.sample(**args)
        return None

    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=None, randseed=None,
               popkeys=None, size=1):  # Keyword arguments are for consistency but not actually used
        # Figure out sample
        if not sample:
            y = self.y
        else:
            if sample == 'new' or self.ysample is None: self.sample(randseed=randseed, size=size)  # msample doesn't exist, make it
            y = self.ysample

        dt = smu.gettvecdt(tvec=tvec, dt=dt, justdt=True)  # Method for getting dt
        outkeys = getoutkeys(self, popkeys)  # Get the list of keys for the output
        if asarray:
            output = np.zeros(len(outkeys))
        else:
            output = sc.odict()

        for pop, key in enumerate(outkeys):  # Loop over each population, always returning an [npops x npts] array
            if key in self.keys():
                yval = y[key]
            else:
                yval = 0.  # Population not present, set to zero
            yinterp = apply_limits(par=self, y=yval, limits=self.limits, dt=dt)
            if asarray:
                output[pop] = yinterp
            else:
                output[key] = yinterp
        return output


class Timepar(Par):
    ''' The definition of a single time-varying parameter, which may or may not vary by population '''

    def __init__(self, t=None, y=None, **defaultargs):
        Par.__init__(self, **defaultargs)
        if t is None: t = odict()
        if y is None: y = odict()
        self.t = t  # Time data, e.g. [2002, 2008]
        self.y = y  # Value data, e.g. [0.3, 0.7]

    def keys(self):
        ''' Return the valid keys for using with this parameter '''
        return self.y.keys()

    def df(self, key=None, data=None):
        '''
        Return t,y data as a data frame; or if data is supplied, replace current t,y values.
        Example: use df() to export data, work with it as a dataframe, and then import it back in:
        screen = P.pars()['screen'].df()
        screen.addrow([2015, 0.3])
        P.pars()['screen'].df(data=screen)
        '''
        if key is None: key = self.keys()[0]  # Pull out first key if not specified -- e.g., 'tot'
        output = dataframe(['t', 'y'], [self.t[key], self.y[key]])
        if data is not None:
            if isinstance(data, dataframe):
                self.t[key] = np.array(data['t'], dtype=float)
                self.y[key] = np.array(data['y'], dtype=float)
            else:
                errormsg = f'Data argument must be a dataframe, not {type(data)}'
                raise Exception(errormsg)
            return None
        else:
            return output

    def interp(self, tvec=None, dt=None, smoothness=None, asarray=True, sample=None, randseed=None,
               popkeys=None):
        """ Take parameters and turn them into model parameters """

        # Validate input
        if tvec is None:
            errormsg = f'Cannot interpolate parameter {self.name} with no time vector specified'
            raise Exception(errormsg)
        tvec, dt = smu.gettvecdt(tvec=tvec, dt=dt)  # Method for getting these as best possible
        if smoothness is None: smoothness = int(1 / dt)  # Handle smoothness
        outkeys = getoutkeys(self, popkeys)  # Get the list of keys for the output

        # Set things up and do the interpolation
        npops = len(outkeys)
        if self.by == 'pship': asarray = False  # Force odict since too dangerous otherwise
        if asarray:
            output = np.zeros((npops, len(tvec)))
        else:
            output = sc.odict()

        for pop, key in enumerate(outkeys):  # Loop over each population, always returning an [npops x npts] array
            if key in self.keys():
                yinterp = sc.smoothinterp(tvec, self.t[key], self.y[key], smoothness=smoothness)  # Use interpolation
                yinterp = apply_limits(par=self, y=yinterp, limits=self.limits, dt=dt)
            else:
                yinterp = np.zeros(len(tvec))  # Population not present, just set to zero
            if asarray:
                output[pop, :] = yinterp
            else:
                output[key] = yinterp
        if npops == 1 and self.by == 'tot' and asarray:
            return output[0, :]  # npops should always be 1 if by==tot, but just be doubly sure
        else:
            return output



#############################################################################################################################
### Functions for handling the parameters
#############################################################################################################################

def getoutkeys(par=None, popkeys=None):
    ''' Small method to decide whether to return 'tot', a subset of population keys, or all population keys '''
    if par.by in ['mpop', 'fpop'] and popkeys is not None:
        return popkeys  # Expand male or female only keys to all
    else:
        return par.keys()  # Or just return the default


def getvalidyears(years, validdata, defaultind=0):
    ''' Return the years that are valid based on the validity of the input data '''
    if sum(validdata):  # There's at least one data point entered
        if len(years) == len(validdata):  # They're the same length: use for logical indexing
            validyears = array(array(years)[validdata])  # Store each year
        elif len(validdata) == 1:  # They're different lengths and it has length 1: it's an assumption
            validyears = array(
                [array(years)[defaultind]])  # Use the default index; usually either 0 (start) or -1 (end)
    else:
        validyears = array([0.0])  # No valid years, return 0 -- NOT an empty array, as you might expect!
    return validyears


def data2init(data=None, parname='prev', keys=None, index=0, **defaultargs):
    ''' Take an array of data return either the first or last (...or some other) non-NaN entry '''

    par = Metapar(y=sc.odict([(key, None) for key in keys]), **defaultargs)  # Create structure -- need key:None for prior
    for row, key in enumerate(keys):
        par.y[key] = sc.sanitize(data[parname][0][row])[index]  # Return the specified index
        parrange = sc.sanitize(data[parname][2][row])[index]-sc.sanitize(data[parname][1][row])[index]
        par.prior[key]['par1'] = par.y[key] # Mean
        par.prior[key]['par2'] = parrange/4 # Assuming normal distribution
    return par


def data2timepar(data=None, years=None, keys=None, defaultind=0, verbose=2, **defaultargs):
    ''' Take an array of data and turn it into default parameters, including priors '''
    # Check that at minimum, name and short were specified, since can't proceed otherwise
    try:
        name, short = defaultargs['name'], defaultargs['short']
    except:
        errormsg = f'Cannot create a time parameter without keyword arguments "name" and "short"! \n\nArguments:\n {defaultargs.items()}'
        raise Exception(errormsg)

    # Process data
    if isinstance(data, dict):  # The entire structure has been passed
        thisdata = data[short]
        years = data['years']
    elif isinstance(data, list):  # Just the relevant entry has been passed
        thisdata = data

    par = Timepar(m=1.0, y=sc.odict(), t=sc.odict(), **defaultargs)  # Create structure
    for row, key in enumerate(keys):
        try:
            validdata = ~np.isnan(thisdata[row])
            par.t[key] = sc.getvaliddata(years, validdata, defaultind=defaultind)
            if sum(validdata):
                par.y[key] = sc.sanitize(thisdata[row])
            else:
                sc.printv(f'data2timepar(): no data for parameter {name}, key {key}', 3, verbose)
                par.y[key] = np.array([0.0])  # Blank, assume zero
                par.t[key] = np.array([0.0])
        except:
            errormsg = f'Error converting time parameter {name}, key {key}'
            sc.printv(errormsg, 1, verbose)
            raise

    return par


def convert_limits(limits=None, tvec=None, dt=None, safetymargin=None, settings=None, verbose=None):
    '''
    Method to calculate numerical limits from qualitative strings.
    '''

    if verbose is None:
        if settings is not None:
            verbose = settings.verbose
        else:
            verbose = 2

    sc.printv('Converting to numerical limits...', 4, verbose)
    if dt is None:
        if settings is not None:
            dt = settings.dt
        else:
            raise Exception('convert_limits() must be given either a timestep or a settings object')
    if safetymargin is None:
        if settings is not None:
            safetymargin = settings.safetymargin
        else:
            sc.printv('Note, using default safetymargin since could not find it', 4, verbose)
            safetymargin = 0.8  # Not that important, so just set safety margin

    # Update dt
    dt = smu.gettvecdt(tvec=tvec, dt=dt, justdt=True)

    # Actually define the rates
    maxrate = safetymargin / dt
    maxpopsize = 1e9
    maxduration = 1000.
    maxmeta = 1000.0
    maxacts = 5000.0
    maxyear = settings.end if settings is not None else 2050.  # Set to a default maximum year

    # It's a single number: just return it
    if sc.isnumber(limits): return limits

    # Just return the limits themselves as a dict if no input argument
    if limits is None:
        return {'maxrate': maxrate, 'maxpopsize': maxpopsize, 'maxduration': maxduration, 'maxmeta': maxmeta, 'maxacts': maxacts}

    # If it's a string, convert to list, but remember this
    isstring = (type(limits) == str)
    if isstring: limits = [limits]  # Convert to list

    # If it's a tuple, convert to a list before converting back at the end
    istuple = (type(limits) == tuple)
    if istuple: limits = list(limits)

    # If list argument is given, replace text labels with numeric limits
    for i, m in enumerate(limits):
        if m == 'maxrate':
            limits[i] = maxrate
        elif m == 'maxpopsize':
            limits[i] = maxpopsize
        elif m == 'maxduration':
            limits[i] = maxduration
        elif m == 'maxmeta':
            limits[i] = maxmeta
        elif m == 'maxacts':
            limits[i] = maxacts
        else:
            limits[i] = limits[i]  # This leaves limits[i] untouched if it's a number

    # Wrap up
    if isstring: return limits[0]  # Return just a scalar
    if istuple:
        return tuple(limits)  # Convert back to a tuple
    else:
        return limits  # Or return the whole list


def apply_limits(y, par=None, limits=None, dt=None, warn=True, verbose=2):
    '''
    A function to intelligently apply limits (supplied as [low, high] list or tuple) to an output.
    Needs dt as input since that determines maxrate.
    '''

    # If parameter object is supplied, use it directly
    parname = ''
    if par is not None:
        if limits is None: limits = par.limits
        parname = par.name

    # If no limits supplied, don't do anything
    if limits is None:
        printv('No limits supplied for parameter "%s"' % parname, 4, verbose)
        return y

    if dt is None:
        raise Exception('No timestep specified: required for convert_limits()')

    # Convert any text in limits to a numerical value
    limits = convert_limits(limits=limits, dt=dt, verbose=verbose)

    # Apply limits, preserving original class -- WARNING, need to handle nans
    if sc.isnumber(y):
        if ~np.isfinite(y): return y  # Give up
        newy = np.median([limits[0], y, limits[1]])
        if warn and newy != y: sc.printv(f'Note, parameter value {parname} reset from {y} to {newy}', 3, verbose)
    elif np.shape(y):
        newy = np.array(y)  # Make sure it's an array and not a list
        infiniteinds = sc.findinds(~np.isfinite(newy))
        infinitevals = newy[infiniteinds]  # Store these for safe keeping
        if len(infiniteinds): newy[infiniteinds] = limits[0]  # Temporarily reset -- value shouldn't matter
        newy[newy < limits[0]] = limits[0]
        newy[newy > limits[1]] = limits[1]
        newy[infiniteinds] = infinitevals  # And stick them back in
        if warn and any(newy != np.array(y)):
            sc.printv(f'Note, parameter {parname} value reset from:\n{y}\nto:\n{newy}', 3, verbose)
    else:
        if warn:
            raise Exception(
                f'Data type {type(y)} not understood for applying limits for parameter {parname}')
        else:
            newy = np.array(y)

    if np.shape(newy) != np.shape(y):
        errormsg = f'Something went wrong with applying limits for parameter {parname}:\ninput and output do not have the same shape:\n{np.shape(y)} vs. {np.shape(newy)}'
        raise Exception(errormsg)

    return newy


## Acts
def balance(which=None, data=None, popkeys=None, limits=None, popsizepar=None, eps=None):
    '''
    Combine the different estimates for the number of acts or protection use and return the "average" value.
    Set which='numacts' to compute for number of acts, which='protection' to compute for protection.
    '''
    if eps is None: eps = 1e-3

    if which not in ['numacts', 'protection']: raise Exception('Can only balance numacts or protection, not "%s"' % which)
    mixmatrix = np.array(data['part'])  # Get the partnerships matrix
    npops = len(popkeys)  # Figure out the number of populations
    symmetricmatrix = np.zeros((npops, npops));
    for pop1 in range(npops):
        for pop2 in range(npops):
            if which == 'numacts': symmetricmatrix[pop1, pop2] = symmetricmatrix[pop1, pop2] + (
                        mixmatrix[pop1, pop2] + mixmatrix[pop2, pop1]) / float(
                eps + ((mixmatrix[pop1, pop2] > 0) + (mixmatrix[pop2, pop1] > 0)))
            if which == 'protection': symmetricmatrix[pop1, pop2] = bool(
                symmetricmatrix[pop1, pop2] + mixmatrix[pop1, pop2] + mixmatrix[pop2, pop1])

    # Decide which years to use -- use the earliest year, the latest year, and the most time points available
    yearstouse = []
    for row in range(npops): yearstouse.append(sc.getvaliddata(data['years'], data[which][row]))
    minyear = np.Inf
    maxyear = -np.Inf
    npts = 1  # Don't use fewer than 1 point
    for row in range(npops):
        minyear = np.minimum(minyear, min(yearstouse[row]))
        maxyear = np.maximum(maxyear, max(yearstouse[row]))
        npts = np.maximum(npts, len(yearstouse[row]))
    if minyear == np.Inf:  minyear = data['years'][0]  # If not set, reset to beginning
    if maxyear == -np.Inf: maxyear = data['years'][-1]  # If not set, reset to end
    ctrlpts = np.linspace(minyear, maxyear, npts).round()  # Force to be integer

    # Interpolate over population acts data for each year
    tmppar = data2timepar(name='tmp', short=which, limits=(0, 'maxacts'), data=data[which],
                          years=data['years'], keys=popkeys, by='pop',
                          verbose=0)  # Temporary parameter for storing acts
    tmpsim = tmppar.interp(tvec=ctrlpts)
    if which == 'numacts': popsize = popsizepar.interp(tvec=ctrlpts)
    npts = len(ctrlpts)

    # Compute the balanced acts
    output = np.zeros((npops, npops, npts))
    for t in range(npts):
        if which == 'numacts':
            smatrix = sc.dcp(symmetricmatrix)  # Initialize
            psize = popsize
            popacts = tmpsim[:, t]
            for pop1 in range(npops): smatrix[pop1, :] = smatrix[pop1, :] * psize[pop1]  # Yes, this needs to be separate! Don't try to put in the next for loop, the indices are opposite!
            for pop1 in range(npops): smatrix[:, pop1] = psize[pop1] * popacts[pop1] * smatrix[:, pop1] / float(eps + sum(smatrix[:,pop1]))  # Divide by the sum of the column to normalize the probability, then multiply by the number of acts and population size to get total number of acts

        # Reconcile different estimates of number of acts, which must balance
        thispoint = np.zeros((npops, npops));
        for pop1 in range(npops):
            for pop2 in range(npops):
                if which == 'numacts':
                    balanced = (smatrix[pop1, pop2] * psize[pop1] + smatrix[pop2, pop1] * psize[pop2]) / (psize[pop1] + psize[pop2])  # There are two estimates for each interaction; reconcile them here
                    thispoint[pop2, pop1] = balanced / psize[pop2]  # Divide by population size to get per-person estimate
                    thispoint[pop1, pop2] = balanced / psize[pop1]  # ...and for the other population
                if which == 'protection':
                    thispoint[pop1, pop2] = (tmpsim[pop1, t] + tmpsim[pop2, t]) / 2.0
                    thispoint[pop2, pop1] = thispoint[pop1, pop2]

        output[:, :, t] = thispoint

    return output, ctrlpts


def make_pars(data=None, verbose=2, die=0):
    """
    Translates the raw data (which were read from the spreadsheet) into
    parameters that can be used in the model. These data are then used to update
    the corresponding model (project). This method should be called before a
    simulation is run.
    """

    sc.printv('Converting data to parameters...', 1, verbose)

    ###############################################################################
    ## Loop over quantities
    ###############################################################################

    pars = sc.odict()

    # Shorten information on which populations are male, which are female, which inject, which provide commercial sex
    pars['male'] = np.array(data['pops']['male']).astype(bool)  # Male populations
    pars['female'] = np.array(data['pops']['female']).astype(bool)  # Female populations

    # Set up keys
    totkey = ['tot']  # Define a key for when not separated by population
    popkeys = data['pops']['short']  # Convert to a normal string and to lower case...maybe not necessary
    fpopkeys = [popkey for popno, popkey in enumerate(popkeys) if data['pops']['female'][popno]]
    mpopkeys = [popkey for popno, popkey in enumerate(popkeys) if data['pops']['male'][popno]]
    pars['popkeys'] = sc.dcp(popkeys)
    pars['age'] = np.array(data['pops']['age'])

    # Read in parameters automatically
    try:
        rawpars = load_par_specs()  # Read the parameters structure
    except Exception as E:
        errormsg = 'Could not load parameter table: "%s"' % repr(E)
        raise Exception(errormsg)

    pars['fromto'], pars['transmatrix'] = load_trans_specs(npops=len(popkeys))  # Read the transitions

    for rawpar in rawpars:  # Iterate over all automatically read in parameters
        sc.printv(f'Converting data parameter {rawpar["short"]}...', 3, verbose)

        try:  # Optionally keep going if some parameters fail

            # Shorten key variables
            partype = rawpar.pop('partype')
            parname = rawpar['short']
            by = rawpar['by']
            fromdata = rawpar['fromdata']
            rawpar['verbose'] = verbose  # Easiest way to pass it in

            # Decide what the keys are
            if by == 'tot':
                keys = totkey
            elif by == 'pop':
                keys = popkeys
            elif by == 'fpop':
                keys = fpopkeys
            elif by == 'mpop':
                keys = mpopkeys
            else:
                keys = []

            # Decide how to handle it based on parameter type
            if partype == 'init':  # Used to determine initial conditions
                if parname == 'init_prev':
                    pars[parname] = data2init(data=data, parname='prev', keys=keys, **rawpar)  # Pull out first available popsize and prevalence point
                elif parname == 'init_pop':
                    pars[parname] = data2init(data=data, parname='pop_size', keys=keys, **rawpar)  # Pull out first available popsize and prevalence point

            elif partype == 'timepar':  # Otherwise it's a regular time par, made from data
                domake = False  # By default, don't make the parameter
                if by != 'pship' and fromdata: domake = True  # If it's not a partnership parameter and it's made from data, then make it
                if domake:
                    pars[parname] = data2timepar(data=data, keys=keys, **rawpar)
                else:
                    pars[parname] = Timepar(y=sc.odict([(key, np.array([np.nan])) for key in keys]),
                                            t=sc.odict([(key, np.array([0.0])) for key in keys]),
                                            **rawpar)  # Create structure

            elif partype == 'meta':  # Scaling parameters over the transitions used during calibration
                pars[parname] = Metapar(y=sc.odict([(key, 1.) for key in keys]), **rawpar)
                # if parname == 'birth_meta':


        except Exception as E:
            errormsg = f'Failed to convert parameter {parname}:\n{repr(E)}'
            if die:
                raise Exception(errormsg)
            else:
                sc.printv(errormsg, 1, verbose)

    ###############################################################################
    ## Tidy up -- things that can't be converted automatically
    ###############################################################################

    # Birth transitions
    npopkeys = len(popkeys)
    birthtransit = np.zeros((npopkeys, npopkeys))
    c = 0
    for pkno, popkey in enumerate(popkeys):
        if data['pops']['female'][pkno]:
            for colno, col in enumerate(data['birthtransit'][c]):
                if sum(data['birthtransit'][c]):
                    birthtransit[pkno, colno] = col / sum(data['birthtransit'][c])
            c += 1
    pars['birthtransit'] = birthtransit

    # Handle acts
    tmpacts, tmpactspts= balance(which='numacts', data=data, popkeys=popkeys, popsizepar=pars['init_pop'])
    tmpcond, tmpcondpts = balance(which='protection', data=data, popkeys=popkeys)

    # Convert matrices to lists of of population-pair keys
    for i, key1 in enumerate(popkeys):
        for j, key2 in enumerate(popkeys):
            if sum(np.array(tmpacts)[i, j, :]) > 0:
                pars['numacts'].y[(key1, key2)] = np.array(tmpacts)[i, j, :]
                pars['numacts'].t[(key1, key2)] = np.array(tmpactspts)
                if key1 in mpopkeys or key1 not in fpopkeys:  # For condom use, only store one of the pair -- and store male first
                    pars['protection'].y[(key1, key2)] = np.array(tmpcond)[i, j, :]
                    pars['protection'].t[(key1, key2)] = np.array(tmpcondpts)

    return pars

