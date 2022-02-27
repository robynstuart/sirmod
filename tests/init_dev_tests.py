''' Dev tests '''

import numpy as np
import sciris as sc
import pylab as pl
import sys
import os

# Add module to paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

# Helper function for getting lognormal parameters
def get_ln_pars(desired_mu, desired_sigma):
    s = np.sqrt(np.log(desired_sigma ** 2 / desired_mu ** 2 + 1))
    scale = desired_mu ** 2 / np.sqrt(desired_sigma ** 2 + desired_mu ** 2)
    return s, scale

def test_pars():
    from sirmod.data import load_data
    from sirmod.parameters import make_pars

    data = load_data(data_sheet='simple.xlsx')
    pars = make_pars(data=data)
    return pars

def test_model():
    from sirmod.model import Model
    model = Model(data_sheet='simple.xlsx')
    model.run()
    return model


def test_fit():
    # Imports
    from sirmod.model import Model
    from sirmod.utils import prep_distance, distance_popsize, distance_numinf, distance_prev
    import pyabc
    import tempfile
    from pyabc.populationstrategy import AdaptivePopulationSize

    # Generate model and data
    model = Model(data_sheet='simple.xlsx')
    _, datay = prep_distance(model)

    # Set up priors
    metaparkeys = ['foi_meta', 'death_other_meta', 'death_meta', 'birth_meta', 'rec_meta']
    metaparkeysfull = [mpk+'_M' for mpk in metaparkeys]+[mpk+'_F' for mpk in metaparkeys]

    s, scale = get_ln_pars(1, 1)
    hyperpars = dict(s=s, scale=scale)
    prior_pars = {k:hyperpars for k in metaparkeysfull}
    prior = pyabc.Distribution(
        **{
            key: pyabc.RV("lognorm", **val)
            for key, val in prior_pars.items()
        },
    )

    distance = pyabc.AdaptiveAggregatedDistance([distance_popsize, distance_prev])

    # Define actual ABC call
    abc = pyabc.ABCSMC(model.abc_run,
                       prior,
                       distance,
                       population_size=50,
                       # population_size=AdaptivePopulationSize(500, 0.15)
                       )

    # Set up database
    db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test.db")
    abc.new(db_path, datay);

    # Run ABC
    history = abc.run(max_nr_populations=500)
    return history, model


def test_multi(sample=False, vals=2):
    from sirmod.model import Model, MultiModel
    model = Model(data_sheet='simple.xlsx')
    mmod = MultiModel(base_model=model)

    # Define metapars
    from scipy.stats import lognorm
    np.random.seed(2022)
    metaparkeys = ['foi_meta', 'death_other_meta', 'death_meta', 'birth_meta', 'rec_meta']
    metaparkeysfull = [mpk+'_M' for mpk in metaparkeys]+[mpk+'_F' for mpk in metaparkeys]
    if sample:
        s, scale = get_ln_pars(1, 1)
        metapars = [{key: lognorm.rvs(s=s, scale=scale) for key in metaparkeysfull},
                    {key: lognorm.rvs(s=s, scale=scale) for key in metaparkeysfull},
                    {key: lognorm.rvs(s=s, scale=scale) for key in metaparkeysfull},
                    {key: lognorm.rvs(s=s, scale=scale) for key in metaparkeysfull}]
    else:
        metapars = [{key: 1 for key in metaparkeysfull}]
        if vals == 2:
            metapars[0]['birth_meta_F'] = .3
            metapars[0]['foi_meta_F'] = 5
            metapars[0]['foi_meta_M'] = 5
            metapars[0]['death_other_meta_F'] = 2
            metapars[0]['death_other_meta_M'] = 2

    mmod.run(metapars=metapars)
    return mmod

def mmod_plot(do_save=True, fname=None):
    fig = pl.figure(figsize=(12, 8))
    colors = sc.gridcolors(2)
    fillzorder = 0  # The order in which to plot things -- fill at the back
    datazorder = 100  # Then data
    linezorder = 200  # Finally, lines
    plotyears = mmod.years

    # Plot prevalence and pop sizes
    for kn, key in enumerate(['pop_size', 'prev']):
        for pn, popkey in enumerate(mmod.results.popkeys):
            n = int(f'{kn}{pn}', base=2) + 1
            pl.subplot(2, 2, n)
            factor = 1. if not key == 'prev' else 100.
            best = factor * mmod.results[key].pops[0, pn, :]
            pl.plot(plotyears, best, color=colors[kn], label=popkey)
            bottom = factor * mmod.results[key].pops[1, pn, :]
            top = factor * mmod.results[key].pops[2, pn, :]
            pl.fill_between(plotyears, bottom, top, facecolor=colors[kn], alpha=.2, label=popkey)

            for y in range(len(mmod.results.datayears)):
                ydata = factor * np.array([mmod.results[key].datapops[1, pn, y], mmod.results[key].datapops[2, pn, y]])
                pl.plot(plotyears[y] * np.array([1, 1]), ydata, c='k', lw=1)
            pl.scatter(plotyears, factor * mmod.results[key].datapops[0, pn, :], c='k', s=30, lw=0,
                       zorder=datazorder)  # Without zorder, renders behind the graph

            pl.title(mmod.results[key].name + ' - ' + popkey)
            sc.boxoff()
            sc.commaticks()
            # pl.gca().set_xticklabels(plotyears.astype(int))
    pl.show()
    if fname is None: fname = 'mmod-plot.png'
    if do_save: sc.savefig(fname+'.png', fig)

    fig = pl.figure(figsize=(12, 8))
    colors = sc.gridcolors(4)
    fillzorder = 0  # The order in which to plot things -- fill at the back
    datazorder = 100  # Then data
    linezorder = 200  # Finally, lines
    plotyears = mmod.years

    # Plot SIRD
    for kn, key in enumerate(['numsus', 'numinf', 'numrec', 'deaths']):
        pl.subplot(2, 2, kn + 1)
        best = factor * mmod.results[key].tot[0, :]
        pl.plot(plotyears, best, color=colors[kn])
        bottom = factor * mmod.results[key].tot[1, :]
        top = factor * mmod.results[key].tot[2, :]
        pl.fill_between(plotyears, bottom, top, facecolor=colors[kn], alpha=.2)

        if mmod.results[key].datatot is not None:
            for y in range(len(mmod.results.datayears)):
                ydata = factor * np.array([mmod.results[key].datatot[1, y], mmod.results[key].datatot[2, y]])
                pl.plot(plotyears[y] * np.array([1, 1]), ydata, c='k', lw=1)
            pl.scatter(plotyears, factor * mmod.results[key].datatot[0, :], c='k', s=30, lw=0,
                       zorder=datazorder)  # Without zorder, renders behind the graph

        pl.title(mmod.results[key].name)
        sc.boxoff()
        sc.commaticks()
        # pl.gca().set_xticklabels(plotyears.astype(int))
    pl.show()
    if fname is None: fname = 'mmod-plot-sir.png'
    if do_save: sc.savefig(fname+'-sir.png', fig)


if __name__ == '__main__':

    # Plotting
    import matplotlib.pyplot as pl
    pl.rcParams['font.size'] = 22
    pl.rcParams['font.family'] = 'Libertinus Sans'

    # Define what to run
    torun = [
        # 'handcal',
        'abccal'
    ]
    run_abc=True # Whether to run ABC or just load a previously run version
    makeplots = False

    # ABC-SMC runs
    if 'abccal' in torun:

        # Either run or load ABC-SMC
        if run_abc:
            history, model = test_fit()
            sc.saveobj('trial_run.obj', history)
        else:
            from sirmod.model import Model
            model = Model(data_sheet='simple.xlsx') # Initalize model again
            history = sc.loadobj('trial_run.obj')

        if makeplots:
            # Plot the simulated trajectories
            df, w = history.get_distribution()
            from sirmod.model import MultiModel
            mmod = MultiModel(base_model=model)
            metapars = df.to_dict(orient='records')
            mmod.run(metapars=metapars)
            mmod_plot(fname='fig5-abc_calibration')

            # Plot the posteriors
            from pyabc.visualization import plot_kde_1d, plot_kde_matrix

            metaparkeys = ['foi_meta', 'death_other_meta', 'death_meta', 'birth_meta', 'rec_meta']
            metaparkeysfull = [mpk + '_M' for mpk in metaparkeys] + [mpk + '_F' for mpk in metaparkeys]
            titles = ['Force of infection', 'All-cause mortality', 'Disease-related mortality', 'Birth rate', 'Recovery rate']
            titlesfull = [t + ' - M' for t in titles] + [t + ' - F' for t in titles]

            # Plot posterior distributions of the FOI metaparameters over generation steps
            pl.rcParams['font.size'] = 22
            pl.rcParams['font.family'] = 'Libertinus Sans'
            fig = pl.figure(figsize=(12, 16))
            colors = pl.cm.GnBu(np.linspace(0.1, 0.9, history.n_populations))
            for npar,parkey in enumerate(metaparkeysfull):
                pl.subplot(5, 2, npar+1)
                ax = pl.gca()
                ax.axvline(1, color="black", linestyle="dotted")
                for t in range(0, history.n_populations, 10):
                    df, w = history.get_distribution(t=t)
                    if len(w) > 0:  # Particles in a model might die out
                        plot_kde_1d(
                            df,
                            w,
                            parkey,
                            ax=ax,
                            label=f"t={t}" if npar==0 else None,
                            c=colors[t],
                            xmin=0,
                            xmax=2,
                            numx=200,
                        )
                ax.set_title(titlesfull[npar])
                ax.set_ylabel('')
                ax.set_xlabel('')
                sc.boxoff()
            fig.legend(title="Generation", fontsize=18, ncol=1, loc="upper right", frameon=False)
            fig.tight_layout(pad=1.0)
            fig.subplots_adjust(right=0.8)
            sc.savefig('fig4-posteriors-all.png', fig)

    # Hand calibration runs
    if 'handcal' in torun:
        for sample in [False, True]:
            mmod = test_multi(sample=sample)
            if makeplots:
                if not sample:  fname = 'fig1-init_calibration'
                else:           fname = 'fig2-2nd_calibration'
                mmod_plot(fname)


