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
    model = Model(data_sheet='simple.xlsx')

    import pyabc

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

    abc = pyabc.ABCSMC(model.run, prior, model.get_fit, population_size=10)

    return abc

if __name__ == '__main__':

    # pars = test_pars()
    # model = test_model()
    abc = test_fit()
    import tempfile
    db_path = os.path.join(tempfile.gettempdir(), "test.db")
    abc.new("sqlite:///" + db_path)
    history = abc.run(minimum_epsilon=0.1, max_nr_populations=10)



    # import tempfile

    # def model(parameter):
    #     return {"data": parameter["mu"] + 0.5 * np.random.randn()}
    # prior = pyabc.Distribution(mu=pyabc.RV("uniform", 0, 5))
    # def distance(x, x0):
    #     return abs(x["data"] - x0["data"])
    # abc = pyabc.ABCSMC(model, prior, distance, population_size=1000)
    # db_path = os.path.join(tempfile.gettempdir(), "test.db")
    # observation = 2.5
    # abc.new("sqlite:///" + db_path, {"data": observation})
    # history = abc.run(minimum_epsilon=0.1, max_nr_populations=10)
