'''
Core code for defining and running the model
'''

# Imports
import numpy as np
import pandas as pd
import sciris as sc
from scipy.stats import truncnorm
from . import utils as smu
from . import parameters as smpars
from . import data as smd
# from . import plotting as smplt


# Define the model
class BaseModel():
    def __init__(self, label=None, verbose=1):
        # Set attributes
        self.label          = label    # The label/name of the simulation
        self.created        = None     # The datetime the sim was created
        self.t              = None     # The current time in the simulation (during execution); outside of sim.step(), its value corresponds to next timestep to be computed
        self.dt             = 0.2  # Timestep
        self.start          = 2015.0  # Default start year
        self.now            = 2022.0  # Default current year
        self.dataend        = 2025.0  # Default end year for data entry
        self.end            = 2035.0  # Default end year for projections
        self.datayears      = sc.inclusiverange(self.start,self.dataend)
        self.projyears      = sc.inclusiverange(self.dataend+1,self.end)
        self.years          = sc.inclusiverange(self.start,self.end)
        self.nyears         = len(self.years)
        self.tvec           = sc.inclusiverange(start=self.start, stop=self.end, step=self.dt)
        self.npts           = len(self.tvec)
        self.ntvec          = np.arange(self.npts)
        self.states         = ['sus', 'inf', 'rec', 'dead', 'dead_other']
        self.nstates        = len(self.states)
        self.healthstatesfull = ['Susceptible', 'Infected', 'Recovered', 'Dead', 'Dead other']
        self.verbose        = verbose
        self.eps            = 1e-3  # A small number used to avoid /0 errors
        return



class Model(BaseModel):

    def __init__(self, data_sheet=None, init_results=True, make_pars=True, make_simpars=False, **kwargs):
        super().__init__(**kwargs)
        # Load datasheet, if provided
        if data_sheet is not None:
            self.load_data(data_sheet=data_sheet, make_pars=make_pars, make_simpars=make_simpars)
        if init_results:    self.init_results()
        else:               self.results = {}
        return

    def load_data(self, data_sheet=None, make_pars=True, make_simpars=False):
        ''' Load in a datasheet '''
        self.data = smd.load_data(data_sheet=data_sheet)
        self.npops = self.data['npops']
        self.popkeys = self.data['pops']['short']
        if make_pars:
            self.pars       = smpars.make_pars(data=self.data)
            if make_simpars:
                self.simpars    = self.interpolate_pars(self.pars)
        return

    def interpolate_pars(self, pars, smoothness=1):
        ''' Interpolate parameters '''

        # Handle inputs and initialization
        simpars = sc.odict()
        keys = list(pars.keys())  # Just get all keys
        if smoothness is None: smoothness = int(1 / dt)

        generalkeys = ['male', 'female', 'popkeys', 'fromto', 'transmatrix']
        staticmatrixkeys = ['birthtransit']

        # Copy default keys by default
        for key in generalkeys: simpars[key] = sc.dcp(pars[key])
        for key in staticmatrixkeys: simpars[key] = sc.dcp(np.array(pars[key]))

        # Loop over requested keys
        for key in keys:  # Loop over all keys
            if isinstance(pars[key], smpars.Par):
                try:
                    simpars[key] = pars[key].interp(tvec=self.tvec, dt=self.dt, popkeys=self.popkeys, smoothness=smoothness)
                except Exception as E:
                    errormsg = f'Could not figure out how to interpolate parameter {key}'
                    errormsg += f'Error: {repr(E)}'
                    raise Exception(errormsg)

        return simpars

    def init_compartments(self):
        ''' Initialize people into compartments '''

        simpars = self.simpars
        init_people = np.zeros((self.nstates, self.npops))  # Initialize
        init_people[0,:] = simpars['init_pop'][:] * (1 - simpars['init_prev'][:])
        init_people[1,:] = simpars['init_pop'][:] * simpars['init_prev'][:]  # Set initial infected population
        return init_people

    def process_acts(self):
        ''' Compute the effective numbers of acts '''

        simpars = self.simpars
        sexactslist = []
        for popkey in simpars['numacts']:
            this = sc.odict()
            this['wholeacts'] = np.floor(self.dt * simpars['numacts'][popkey])
            this['fracacts'] = self.dt * simpars['numacts'][popkey] - this['wholeacts']  # Probability of an additional act

            if simpars['protection'].get(popkey) is not None:
                condkey = simpars['protection'][popkey]
            elif simpars['protection'].get((popkey[1], popkey[0])) is not None:
                condkey = simpars['protection'][(popkey[1], popkey[0])]
            else:
                errormsg = label + f'Cannot find condom use between {popkey[0]} and {popkey[1]}, assuming there is none.'
                sc.printv(errormsg, 1, verbose)
                condkey = 0.0

            this['protection'] = 1.0 - condkey * simpars['effprotection']
            this['pop1'] = self.popkeys.index(popkey[0])
            this['pop2'] = self.popkeys.index(popkey[1])
            if simpars['male'][this['pop1']] and simpars['female'][this['pop2']]:
                this['trans'] = simpars['transmfi']
            elif simpars['female'][this['pop1']] and simpars['male'][this['pop2']]:
                this['trans'] = simpars['transmfr']
            else:
                errormsg = label + f'Not able to figure out the sex of {popkey[0]} and {popkey[0]}'
                sc.printv(errormsg, 3, verbose)
                this['trans'] = (simpars['transmfi'] + simpars['transmfr']) / 2.0

            sexactslist.append(this)

            # Error checking
            for key in ['wholeacts', 'fracacts', 'protection']:
                if not (all(this[key] >= 0)):
                    errormsg = label + f'Invalid sexual behavior parameter {popkey}: values are:\n{this[key]}'
                    sc.printv(errormsg, 1, verbose)
                    this[key][this[key] < 0] = 0.0  # Reset values

        for i, this in enumerate(sexactslist): sexactslist[i] = tuple(
            [this['pop1'], this['pop2'], this['wholeacts'], this['fracacts'], this['protection'], this['trans']])

        return sexactslist

    def process_births(self):
        ''' Allocate births '''
        birth = self.simpars['birth'] * self.dt
        birthtransit = self.simpars['birthtransit']
        birthslist = []
        for p1 in range(self.npops):
            for p2 in range(self.npops):
                birthrates = birthtransit[p1, p2] * birth[p1, :]
                if birthrates.any():
                    birthslist.append(tuple([p1, p2, birthrates]))
        motherpops = set([thisbirth[0] for thisbirth in birthslist])
        return birthslist, motherpops

    def init_results(self):
        ''' Store the results '''
        self.results = sc.objdict()

        # Store some general useful things
        self.results['popkeys'] = self.popkeys
        self.results['inci']            = smu.Result('New infections')
        self.results['births']          = smu.Result('Total births)')
        self.results['deaths']          = smu.Result('Disease-related deaths')
        self.results['deaths_other']    = smu.Result('Non-disease-related deaths')
        self.results['newrec']          = smu.Result('New recoveries')
        self.results['prev']            = smu.Result('Prevalence (%)', ispercentage=True)
        self.results['force']           = smu.Result('Incidence (per 100 p.y.)', ispercentage=True)
        self.results['pop_size']        = smu.Result('Population size')
        self.results['numsus']          = smu.Result('Number currently susceptible')
        self.results['numinf']          = smu.Result('Number currently infected')
        self.results['numrec']          = smu.Result('Number currently recovered')

        # Add the results that have data associated with them
        trailingnans = np.zeros((3,self.npops,len(self.projyears))) + np.nan # 3 represents dimensions of [best, low, high]
        self.results['prev'].datapops       = np.append(np.array(self.data['prev']), trailingnans, axis=2)
        self.results['pop_size'].datapops   = np.append(np.array(self.data['pop_size']), trailingnans, axis=2) # 3 represents dimensions of [best, low, high]
        self.results['pop_size'].datatot    = self.results['pop_size'].datapops.sum(axis=1)
        self.results['numinf'].datapops     = self.results['prev'].datapops*self.results['pop_size'].datapops
        self.results['numinf'].datatot      = self.results['numinf'].datapops.sum(axis=1)
        self.results['prev'].datatot        = self.results['numinf'].datatot / self.results['pop_size'].datatot

        # Raw arrays -- reporting annual quantities (so need to divide by dt!)
        raw = sc.odict()
        raw['inci']         = np.zeros((self.npops, self.npts))  # Total incidence acquired by each population
        raw['newrec']       = np.zeros((self.npops, self.npts)) # Recoveries by population
        raw['births']       = np.zeros((self.npops, self.npts))  # Total number of births to each population
        raw['deaths']       = np.zeros((self.npops, self.npts))  # Number of deaths per timestep
        raw['deaths_other'] = np.zeros((self.npops, self.npts))  # Number of other deaths per timestep

        return raw

    def process_results(self, raw, allpeople, data=None, annual=True, doround=True):
        ''' Process results for plots and other outputs'''

        # Initialize
        if annual is False: # Decide what to do with the time vector
            indices = np.arange(len(self.ntvec)) # Use all indices
            ntvec   = sc.dcp(self.ntvec)
        else:
            indices = np.arange(0, len(self.ntvec), int(round(1.0/(self.tvec[1]-self.tvec[0]))))
            ntvec   = self.ntvec[indices]

        self.results['tvec'] = self.tvec[indices]
        self.results['ntvec'] = ntvec
        self.results['datayears'] = self.datayears

        # Actually do calculations
        for key in ['inci', 'births', 'deaths', 'deaths_other', 'newrec']:
            self.results[key].pops = raw[key][:, indices]
            self.results[key].tot = raw[key][:, indices].sum(axis=0)

        self.results['force'].pops      = raw['inci'][:, indices] / (self.eps + allpeople[0, :, indices].transpose())
        self.results['force'].tot       = raw['inci'][:, indices].sum(axis=0) / (self.eps + allpeople[0, :, indices].sum(axis=1))
        self.results['numinf'].pops     = allpeople[1, :, indices].transpose()
        self.results['numinf'].tot      = allpeople[1, :, indices].sum(axis=1)
        self.results['numsus'].pops     = allpeople[0, :, indices].transpose()
        self.results['numsus'].tot      = allpeople[0, :, indices].sum(axis=1)
        self.results['numrec'].pops     = allpeople[2, :, indices].transpose()
        self.results['numrec'].tot      = allpeople[2, :, indices].sum(axis=1)

        # Add the results that have data associated with them
        self.results['prev'].pops           = allpeople[1, :, indices].transpose() / (self.eps + allpeople[:, :, indices].sum(axis=0))
        self.results['prev'].tot            = allpeople[1, :, indices].sum(axis=1) / (self.eps + allpeople[:, :, indices].sum(axis=(0, 1)))
        self.results['pop_size'].pops       = allpeople[:, :, indices].sum(axis=0)
        self.results['pop_size'].tot        = allpeople[:, :, indices].sum(axis=(0, 1))

        # Now add the data distributions
        datadists = {'prev': 'trunc_norm', 'pop_size': 'norm'}
        for datatype, dist in datadists.items():
            datadist = sc.objdict()
            for pn, popkey in enumerate(self.pars['popkeys']):
                years = sc.getvaliddata(self.years, ~np.isnan(self.results[datatype].datapops[0, pn, :]), defaultind=0)
                best = sc.sanitize(self.results[datatype].datapops[0, pn, :])
                low = sc.sanitize(self.results[datatype].datapops[1, pn, :])
                high = sc.sanitize(self.results[datatype].datapops[2, pn, :])
                par2 = abs(high - low) / 4 # This is because low and high are reversed for popsize
                datadist[popkey] = {years[i]: dict(dist=dist, par1=best[i], par2=par2[i], lower_clip=0, upper_clip=1) for i in range(len(best))}
            self.results[datatype].datadist = datadist

    def step(self, t, people, raw, sexactslist, birthslist, motherpops, metapars=None):
        ''' Calculate probability of getting infected '''

        # Shorten pars
        death               = self.simpars['death'][:, t] * self.dt
        death_other         = self.simpars['death_other'][:, t] * self.dt
        recrate             = 1. - np.exp(-self.dt / (np.maximum(self.eps, self.simpars['duration'][t])))
        thistransit         = sc.dcp(self.simpars['transmatrix'])
        new_people          = np.zeros(people.shape)

        # Deal with metapars
        if metapars is None:
            foi_meta            = np.ones(self.npops)
            death_other_meta    = np.ones(self.npops)
            death_meta          = np.ones(self.npops)
            birth_meta          = np.ones(self.npops)
            rec_meta            = np.ones(self.npops)
        else:
            foi_meta            = np.array([metapars['foi_meta_M'],metapars['foi_meta_F']])
            death_other_meta    = np.array([metapars['death_other_meta_M'],metapars['death_other_meta_F']])
            death_meta          = np.array([metapars['death_meta_M'],metapars['death_meta_F']])
            birth_meta          = np.array([metapars['birth_meta_M'],metapars['birth_meta_F']])
            rec_meta            = np.array([metapars['rec_meta_M'],metapars['rec_meta_F']])

        # Calculate the FOI using behavioral data
        forceinffull = np.ones((self.npops, self.npops)) # First dimension: infection acquired by (pop). Second dimension: infection caused by (pop)

        # Loop over all partnerships to get the probability of pop1 getting infected by pop2
        for pop1, pop2, wholeacts, fracacts, protection, thistrans in sexactslist:
            allprev = people[1, :] / people[:, :].sum(axis=0)
            thisforceinfsex = 1 - fracacts[t] * thistrans[t] * protection[t] * allprev[pop2]
            if wholeacts[t]>0: thisforceinfsex *= np.power(1 - thistrans[t] * protection[t] * allprev[pop2], int(wholeacts[t]))
            forceinffull[pop1, pop2] *= thisforceinfsex

            if not ((forceinffull[pop1, pop2] >= 0).all()):
                errormsg = label + f'Force-of-infection is invalid between populations {popkeys[pop1]} and {popkeys[pop2]}, time {tvec[t]}, FOI:\n{forceinffull[pop1, :, pop2]})'
                raise Exception(errormsg)

        prob_acquire_inf    = (1-forceinffull).sum(axis=1)  # Probability of each population acquiring an infection
        prob_cause_inf      = (1-forceinffull).sum(axis=0)  # Probability of each population causing an infection

        # Safety checks: don't let too many people move - warning, not tested
        if (rec_meta * recrate + death_other_meta * death_other + death_meta * death>1).any():
            scalefactor = rec_meta * recrate + death_other_meta * death_other + death_meta * death
            rec_meta /= scalefactor
            death_other_meta /= scalefactor
            death_meta /= scalefactor
        if (foi_meta*prob_acquire_inf + death_other_meta*death_other>1).any():
            scalefactor = foi_meta*prob_acquire_inf + death_other_meta*death_other
            foi_meta /= scalefactor
            death_other_meta /= scalefactor

        # Define transitions
        thistransit[0, 0, :] = 1. - foi_meta*prob_acquire_inf - death_other_meta*death_other            # sus to sus
        thistransit[0, 1, :] = foi_meta*prob_acquire_inf                                                # sus to inf
        thistransit[0, 4, :] = death_other_meta*death_other                                             # sus to dead_other
        thistransit[1, 1, :] = 1. - rec_meta*recrate - death_other_meta*death_other - death_meta*death  # inf to inf
        thistransit[1, 2, :] = rec_meta*recrate                                                         # inf to rec
        thistransit[1, 3, :] = death_meta*death                                                         # inf to dead
        thistransit[1, 4, :] = death_other_meta*death_other                                             # inf to dead_other
        thistransit[2, 2, :] = 1. - death_other_meta*death_other                                        # rec to rec
        thistransit[2, 4, :] = death_other_meta*death_other                                             # rec to dead

        # Add births
        fsums = dict()
        for p1 in motherpops:
            fsumpop     = people[:, p1]
            fsums[p1]   = fsumpop[:].sum()
        for p1, p2, birthrates in birthslist:
            thisbirthrate = birthrates[t] * birth_meta[p1]
            popbirths = thisbirthrate * fsums[p1]
            raw['births'][p2, t] += popbirths / self.dt

        # Store results
        raw['inci'][:, t]           = (people[0, :] * thistransit[0, 1, :]).sum()/self.dt
        raw['newrec'][:, t]         = (people[1, :] * thistransit[1, 2, :]).sum()/self.dt
        raw['deaths'][:, t]         = (people[:, :] * thistransit[:, 3, :]).sum()/self.dt
        raw['deaths_other'][:, t]   = (people[:, :] * thistransit[:, 4, :]).sum()/self.dt

        # Shift people
        for s in range(self.nstates):
            new_people += people[s, :] * thistransit[s, :, :]
        new_people[0,:] += raw['births'][:, t]

        return new_people, raw

    def run(self, metapars=None, verbose=1, keep_raw=True):
        ''' Run the model once '''

        # Initialize results, acts, births, and people
        self.simpars            = self.interpolate_pars(self.pars)
        raw                     = self.init_results()
        sexactslist             = self.process_acts()
        birthslist, motherpops  = self.process_births()
        allpeople = np.zeros((self.nstates, self.npops, self.npts))

        for t in self.ntvec:  # Loop over time
            sc.printv(f'Timestep {t+1} of {self.npts}', 4, verbose)
            if t==0:
                allpeople[:,:,0]    = self.init_compartments()
            else:
                new_people, raw     = self.step(t, allpeople[:,:,t-1], raw, sexactslist, birthslist, motherpops, metapars=metapars)
                allpeople[:,:,t]    = new_people

        # Finalize and process
        if keep_raw: self.raw = raw
        self.people = allpeople
        self.process_results(raw, allpeople)

        return

    def abc_run(self, metapars=None, verbose=1):
        ''' A slightly different API for running the model with ABC '''
        self.run(metapars=metapars,verbose=verbose)
        modely, _ = smu.prep_distance(self)
        return modely

    def plot(self, toplot=None):
        ''' Plot the outputs of the model '''
        pass


##############################################################################
# Running functions
##############################################################################
def single_run(model, ind=0, run_args=None, verbose=1, **kwargs):
    '''
    Convenience function to perform a single model run.
    '''
    run_args = sc.mergedicts({'verbose':verbose}, run_args, kwargs)
    if verbose is None:
        verbose = model.verbose
    if not model.label:
        model.label = f'Model {ind}'
    if verbose>=1:
        print(f'Running {model.label}')
    model.run(**run_args)
    return model


def multi_run(model, n_runs=4, metapars=None, run_args=None, par_args=None, parallel=True, n_cpus=None, verbose=1, **kwargs):
    ''' For running multiple runs in parallel '''

    # Handle inputs
    par_args = sc.mergedicts({'ncpus':n_cpus}, par_args) # Handle blank

    # Handle metapars
    if metapars is None:
        metapars = {}
    else:
        n_runs = len(metapars) # Reset and get from length of list instead

    # Run the sims
    if isinstance(model, Model): # One model
        iterkwargs = dict(ind=np.arange(n_runs))
        iterkwargs.update(dict(metapars=metapars))
        kwargs = dict(model=model, verbose=verbose, run_args=run_args)
    elif isinstance(sim, list): # List of sims
        iterkwargs = dict(model=model, ind=np.arange(len(model)))
        kwargs = dict(verbose=verbose, run_args=run_args)
    else:
        errormsg = f'Must be Model object or list, not {type(model)}'
        raise TypeError(errormsg)

    # Actually run
    if parallel:
        models = sc.parallelize(single_run, iterkwargs=iterkwargs, kwargs=kwargs, **par_args) # Run in parallel
    else: # Run in serial, not in parallel
        models = []
        n_models = len(list(iterkwargs.values())[0]) # Must have length >=1 and all entries must be the same length
        for s in range(n_models):
            this_iter = {k:v[s] for k,v in iterkwargs.items()} # Pull out items specific to this iteration
            this_iter.update(kwargs) # Merge with the kwargs
            this_iter['model'] = sc.dcp(this_iter['model']) # Ensure we have a fresh sim; this happens implicitly on pickling with multiprocessing
            model = single_run(**this_iter) # Run in series
            models.append(model)

    return models



class MultiModel(BaseModel):
    def __init__(self, models=None, base_model=None, **kwargs):
        super().__init__(**kwargs)
        self.models     = models
        self.base_model = base_model
        self.run_args  = sc.mergedicts(kwargs)
        return

    def __len__(self):
        try:    return len(self.models)
        except: return 0

    def run(self, **kwargs):
        ''' Set up multiple runs of the model with different metapars'''

        if self.models is None:
            models = self.base_model
        else:
            models = self.models

            # Handle missing labels
            for m,model in enumerate(models):
                if model.label is None:
                    model.label = f'Model {m}'

        # Run
        kwargs = sc.mergedicts(self.run_args, kwargs)
        self.models = multi_run(models, **kwargs)
        self.process_results()

        return self


    def process_results(self, reduce=True, quantiles=None, doround=True):
        ''' Process results of multiruns for plots and other outputs'''

        # Initialize
        if quantiles is None: quantiles = [0.5, 0.25, 0.75]
        npops = self.base_model.npops
        nyears = self.base_model.nyears

        # Compile results
        allresults = sc.objdict()
        ref_model = self.models[0]
        for key,res in ref_model.results.items():
            if type(res) == smu.Result: # Results need to be aggregated
                allresults[key]         = sc.dcp(res) # Initially just copy the whole thing
                allresults[key].pops    = np.zeros((len(self),npops,nyears))
                allresults[key].tot     = np.zeros((len(self),nyears))
                for m, model in enumerate(self.models): # Now overwrite simulation results
                    allresults[key].pops[m, :, :] = model.results[key].pops
                    allresults[key].tot[m, :]     = model.results[key].tot
            else: # Just copy over other things, like tvecs
                allresults[key] = sc.dcp(res)

        # Optionally reduce
        if reduce:
            self.results = sc.objdict()
            for key,res in allresults.items():
                if type(res) == smu.Result:  # Results need to be aggregated
                    self.results[key]       = sc.dcp(res)  # Initially just copy the whole thing
                    self.results[key].pops  = np.quantile(res.pops,quantiles,axis=0)
                    self.results[key].tot   = np.quantile(res.tot,quantiles,axis=0)
                else: # Just copy over other things, like tvecs
                    self.results[key] = sc.dcp(res)
            self.reduced=True

        else:
            self.reduced=False
            self.results = allresults

        self.results.type = 'multi'
        self.results.nruns = len(self)

    def plot(self):
        ''' Plot the outputs of the model '''

        pass
