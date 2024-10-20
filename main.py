from model import Model
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle as pkl
import seaborn as sns
import yaml

# Define plotting functions --------------------------------------------------------------------------------------------

def plot_relerror(model):
    """
    This function plots the relative errors of the New/Total RNA ratios.
    :param model: model of RNA metabolism
    :return: None
    """

    pal = sns.color_palette('Set2')
    pal = pal.as_hex()

    times = model.times
    nuc_ratios = model.meas_nuc
    cyt_ratios = model.meas_cyt
    nuc_predicted = model.nuc_newtotal_mlest
    cyt_predicted = model.cyt_newtotal_mlest

    nuc_relerr = (np.abs(nuc_ratios - nuc_predicted) / nuc_ratios) * 100
    cyt_relerr = (np.abs(cyt_ratios - cyt_predicted) / cyt_ratios) * 100

    fig, ax = plt.subplots()
    sns.lineplot(x=times, y=nuc_relerr, label='Nucleus', color=pal[1], ax=ax)
    sns.lineplot(x=times, y=cyt_relerr, label='Cytosol', color=pal[2], ax=ax)
    ax.set_ylabel(r'Relative Error of $\frac{New}{Total}$ RNA ratios [%]')
    ax.set_xlabel('Time [min]')
    ax.axhline(y=5, linestyle='--', color='grey')
    ax.grid(which='major', axis='y')
    ax.set_yscale('log')
    ax.set_ylim(0.01, 100)
    ax.set_yticks(ticks = [0.01, 0.01, 0.1, 1, 5, 10, 100, 1000], labels=[0.01, 0.01, 0.1, 1, 5, 10, 100, 1000])
    fig.tight_layout()
    fig.savefig('data/relativeerror_lineplot.pdf')
    fig.savefig('data/relativeerror_lineplot.png')
    plt.close()

    return None

def plot_ratios(model):
    """
    This function plots the the New/Total RNA ratios over time.
    :param model: model of RNA metabolism
    :return: None
    """

    pal = sns.color_palette('Set2')
    pal = pal.as_hex()
    times = model.times
    nuc_ratios = model.meas_nuc * 100
    cyt_ratios = model.meas_cyt * 100
    nuc_predicted = model.nuc_newtotal_mlest * 100
    cyt_predicted = model.cyt_newtotal_mlest * 100

    fig, ax = plt.subplots()
    sns.lineplot(x=times, y=nuc_ratios, label='Nucleus (measured)', color=pal[1], ax=ax)
    sns.lineplot(x=times, y=nuc_predicted, linestyle="--",  label='Nucleus (predicted)', color=pal[1], ax=ax)
    sns.lineplot(x=times, y=cyt_ratios, color=pal[2], label='Cytosol (measured)', ax=ax)
    sns.lineplot(x=times, y=cyt_predicted, color=pal[2], linestyle='--',  label='Cytosol (predicted)', ax=ax)
    ax.set_ylabel(r'$\frac{New}{Total}$ RNA ratios [%]')
    ax.set_xlabel('Time [min]')
    ax.grid(which='major', axis='y')
    fig.tight_layout()
    fig.savefig('data/ratios_lineplot.pdf')
    fig.savefig('data/ratios_lineplot.png')
    plt.close()

    return None

# Define Wrapper function ----------------------------------------------------------------------------------------------

def modelmetabol():
    """
    Wrapper of the RNA metabolism model.
    :return: None
    """

    print("Modeling RNA metabolism ...")

    with open("docs/parameters.yml", "r") as parameterfile:
        params = yaml.load(parameterfile, Loader=yaml.SafeLoader)

    res = Model(
        parameters=[params['mu'], params['nu'], params['tau'], params['lam']],
        readlength=params['readlength'],
        coverage= params['coverage'],
        labelingefficiency=params['labelingefficiency'],
        times=np.arange(params['times'][0], params['times'][1], params['times'][2])
    )

    if not os.path.exists('data'):
        os.makedirs('data')

    plot_relerror(res)
    plot_ratios(res)
    with open(r"data/model.pkl", "wb") as output_file:
        pkl.dump(res, output_file)

    print("Done.")

    return None

# Run the model and grab a coffee --------------------------------------------------------------------------------------

if __name__ == "__main__":
    modelmetabol()
