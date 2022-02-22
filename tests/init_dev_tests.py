''' Dev tests '''

import numpy as np
import sciris as sc
import pylab as pl
import sys
import os

# Add module to paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

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

    from sirmod.model import Model
    from sirmod.utils import prep_distance, distance

    model = Model(data_sheet='simple.xlsx')
    _, datay = prep_distance(model)

    return model, datay, distance

if __name__ == '__main__':

    model, datay, distance = test_fit()

    import pyabc
    import tempfile
    from pyabc.populationstrategy import AdaptivePopulationSize

    # Set up prior
    metaparkeys = ['foi_meta', 'death_other_meta', 'death_meta', 'birth_meta', 'rec_meta']
    metaparkeysfull = [mpk+'_M' for mpk in metaparkeys]+[mpk+'_F' for mpk in metaparkeys]
    hyperpars = dict(s=1, loc=0, scale=1)
    prior_pars = {k:hyperpars for k in metaparkeysfull}
    prior = pyabc.Distribution(
        **{
            key: pyabc.RV("lognorm", **val)
            for key, val in prior_pars.items()
        },
    )

    # Define actual ABC call
    abc = pyabc.ABCSMC(model.abc_run,
                       prior,
                       distance,
                       population_size=20,
                       # population_size=AdaptivePopulationSize(500, 0.15)
                       )

    # Set up database
    db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test.db")
    abc.new(db_path, {"Y": datay});

    # Run ABC
    history = abc.run(max_nr_populations=15)

    # Plotting
    from pyabc.visualization import plot_kde_1d, plot_kde_matrix, plot_data_callback
    import matplotlib.pyplot as plt

    # # Plot posterior distributions of the FOI metaparameters over generation steps
    # fig, axes = plt.subplots(2)
    # fig.set_size_inches((6, 6))
    # axes = axes.flatten()
    # axes[0].axvline(1, color="black", linestyle="dotted")
    # axes[1].axvline(1, color="black", linestyle="dotted")
    # for m, ax in enumerate(axes):
    #     for t in range(0, history.n_populations, 2):
    #         df, w = history.get_distribution(m=0, t=t)
    #         if len(w) > 0:  # Particles in a model might die out
    #             plot_kde_1d(
    #                 df,
    #                 w,
    #                 ["foi_meta_M", "foi_meta_F"][m],
    #                 ax=ax,
    #                 label=f"t={t}",
    #                 xmin=0,
    #                 xmax=4,
    #                 numx=200,
    #             )
    #     ax.set_title(f"Theta {m + 1}")
    # axes[0].legend(title="Generation", loc="upper left", bbox_to_anchor=(1, 1))
    #
    # fig.tight_layout()
    # fig.show()

    # # Pairwise plot of joint posterior distributions
    # df, w = history.get_distribution(m=0)
    # _, fig = plot_kde_matrix(df, w)
    # fig.show()


    # Plot the simulated trajectories
    _, ax = plt.subplots()

    def plot_data(sum_stat, weight, ax, **kwargs):
        """Plot a single trajectory"""
        ax.plot(measurement_times, sum_stat['X_2'], color='grey', alpha=0.1)

    def plot_mean(sum_stats, weights, ax, **kwargs):
        """Plot mean over all samples"""
        weights = np.array(weights)
        weights /= weights.sum()
        data = np.array([sum_stat['X_2'] for sum_stat in sum_stats])
        mean = (data * weights.reshape((-1, 1))).sum(axis=0)
        ax.plot(measurement_times, mean, color='C2', label='Sample mean')


    ax = plot_data_callback(h, plot_data, plot_mean, ax=ax)

    plt.plot(true_trajectory, color="C0", label='Simulation')
    plt.scatter(measurement_times, measurement_data, color="C1", label='Data')
    plt.xlabel('Time $t$')
    plt.ylabel('Measurement $Y$')
    plt.title('Conversion reaction: Simulated data fit')
    plt.legend()
    plt.show()


