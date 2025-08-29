# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 13:47:02 2025

@author: trist
"""

import os, json, re, shutil
from os.path import join, exists, isdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bolthausen import get_expected_Tp
from tqdm import tqdm
from tree_analysis import get_Tp_from_T2

def find_Ns(path, suffix, verbose=True):
    """
    Get all the population sizes simulated for a given model.

    Parameters
    ----------
    path : str
        Path where T2 trees are recorded.
    suffix : str
        Small string with only letters and numbers specifying which model is
        simulated.
    verbose : bool, optional
        Choose if the graph of the number of simulations for each population 
        size should be ploted. The default is True.

    Returns
    -------
    Ns_list : list of ints
        List of simulated population sizes for this model.

    """
    dirs = os.listdir(join(path, suffix))
    Ns = {}
    for dir in dirs:
        num_files = len(os.listdir(join(path, suffix, dir)))
        if num_files > 0:
            Ns[int(dir)] = num_files
        else:
            print(f'Removed folder for {dir}{suffix} because empty.')
            os.remove(join(path, suffix, dir))
    Ns_list = sorted(Ns.keys())
    if verbose:
        plt.semilogx(Ns_list, [Ns[N] for N in Ns_list], '+')
        plt.xlabel('$N$')
        plt.ylabel('Number of\nsimulations')
        plt.grid(which='both', alpha=.4)
        plt.title('Number of simulations for each population size $N$.')
        plt.tight_layout()
        plt.show()
    return Ns_list


def get_or_write_avg(T2_path, data_path, suffix, N, p, timestamp=''):
    """
    Compute and write the average Tp of a given model for a given population
    size, or just read it if it was already computed.

    Parameters
    ----------
    T2_path : str
        Path of the folder where the T2 dictionnaries are stored.
    data_path : str
        Path of the folder where the average Tp will be writted or read.
    suffix : str
        Small string with only letters and numbers specifying which model is
        simulated.
    N : int
        Size of the population.
    p : int
        Size of the p-sample.

    Raises
    ------
    ValueError
        Folder is empty.

    Returns
    -------
    avg : float
        Average Tp for the model and population size.

    """
    #If file already exists, read just read data
    filename = f'Brunet2 {N}{suffix}_T{p}_{timestamp}_avg'
    if exists(join(data_path, filename)):
        with open(join(data_path, filename), 'r') as f:
            avg = float(f.read())
    else:
        if not isdir(join(T2_path, suffix, str(N))):
            print(f'Folder not found for {suffix}{N}.')
            return np.nan
        #Otherwise, loop over each corresponding T2 file and add the data
        print(f'Computing average T{p} for {N}{suffix}...   {timestamp}')
        if timestamp == '':
            values = [get_or_write_avg(T2_path, data_path, suffix, N, p, 
                                       timestamp=timestamp) for timestamp\
                      in os.listdir(join(T2_path, suffix, str(N)))]
        else:
            values = []
            for T2file in os.listdir(join(T2_path, suffix, str(N), timestamp)):
                with open(join(T2_path, suffix, str(N), timestamp, T2file), 'r') as f:
                    try:
                        T2 = json.load(f)
                    except json.JSONDecodeError:
                        print(f'JSONDecodeError with file {join(T2_path, suffix, str(N), timestamp, T2file)}')
                        continue
                if p == 2:
                    values.extend(T2.values())
                else:
                    Tp = get_Tp_from_T2(T2, p)
                    values.extend(Tp.values())
        if len(values) == 0:
            raise ValueError('values is empty. No files detected.')
            
        #Compute the average of all the data.
        avg = np.mean(values)
               
        #If p>2, we rescale directly here.
        if p != 2 and timestamp != '':
            avg /= get_or_write_avg(T2_path, data_path, suffix, N, 2,
                                    timestamp=timestamp)
        
        #Write down average to not have to recompute it next time.
        with open(join(data_path, filename), 'w') as f:
            f.write(str(avg))
    return avg

def get_or_write_SE(T2_path, data_path, suffix, N, p):
    """
    Compute and write the standard error for the average Tp of a given model 
    for a given population size, or just read it if it was already computed.
    
    Parameters
    ----------
    T2_path : str
        Path of the folder where the T2 dictionnaries are stored.
    data_path : str
        Path of the folder where the average Tp will be writted or read.
    suffix : str
        Small string with only letters and numbers specifying which model is
        simulated.
    N : int
        Size of the population.
    p : int
        Size of the p-sample.
        
    Returns
    -------
    SE : float
        Standard error for the average Tp for the model and population size.
    """
    filename = f'Brunet2 {N}{suffix}_T{p}_SE'
    if exists(join(data_path, filename)):
        with open(join(data_path, filename), 'r') as f:
            SE = float(f.read())
    else:
        print(f'Computing SE for T{p} for {N}{suffix}...')
        values = np.array([get_or_write_avg(T2_path, data_path, suffix, N, p, 
                                   timestamp=datetime) for datetime in \
                  os.listdir(join(T2_path, suffix, str(N)))])
        sigma_X = np.sqrt(np.mean((values - np.mean(values))**2))
        SE = sigma_X/len(values)**.5
        with open(join(data_path, filename), 'w') as f:
            f.write(str(SE))
    return SE
            
            
def curve_T2(Ns, T2_path, data_path, suffix):
    """
    Curve showing the time-scale of the model.
    """
    avgs = []
    for N in Ns:
        avgs.append(get_or_write_avg(T2_path, data_path, suffix, N, 2))
    plt.loglog(Ns, avgs, 'k-', alpha=0.2)
    plt.loglog(Ns, avgs, 'k+', label=r'Empirical values')
    plt.loglog(Ns, Ns, '-.', label='Kingman')
    plt.loglog(Ns, [np.log(N)**2 for N in Ns], '--', label='B-S')
    plt.loglog(Ns, [np.log(N)**3 for N in Ns], '--', label='B-S')
    plt.grid(which='both', alpha=.4)
    plt.title('Time-scale of the model')
    plt.xlabel('N')
    plt.ylabel(r'$\langle T_2\rangle$')
    plt.legend()
    plt.show()
    
    
    
def curves_Tp(Ns, T2_path, data_path, suffix, p_max=5, display_title=True):
    """
    Curve analyzing the model and comparing it to the Kingman coalescent
    and to the Bolthausen-Sznitman coalescent.
    """
    colors = ['b', 'r', 'y', 'k', 'm']
    for p in range(3, p_max+1):
        avgs = []
        q1, q2 = [], []
        for N in Ns:
            avg = get_or_write_avg(T2_path, data_path, suffix, N, p)
            avgs.append(avg)
            SE = get_or_write_SE(T2_path, data_path, suffix, N, p)
            k = 2*5**.5
            q1.append(avg + k*SE)
            q2.append(avg - k*SE)
        if p_max == 3:
            plt.semilogx(Ns, [2-2/p]*len(Ns), 'b--', alpha=1, label='Kingman')
            plt.semilogx(Ns, [get_expected_Tp(p),]*len(Ns), 'r-.', label='Bolthausen-Sznitman')
            plt.semilogx(Ns, avgs, 'g+', label='Empirical')
            plt.semilogx(Ns, avgs, 'g-', alpha=.2)
            plt.fill_between(Ns, q1, q2, color='g', alpha=.05, label='95% confidence interval')
        else:
            plt.fill_between(Ns, q1, q2, color=colors[p-3], alpha=.05, label='95% confid. interval')
            plt.semilogx(Ns, [2-2/p]*len(Ns), f'{colors[p-3]}--', alpha=1, label='Kingman')
            plt.semilogx(Ns, [get_expected_Tp(p),]*len(Ns), f'{colors[p-3]}-.', label='Bolthausen-Sznitman')
            plt.semilogx(Ns, avgs, f'{colors[p-3]}+', label='Empirical values')
            plt.semilogx(Ns, avgs, f'{colors[p-3]}-', alpha=.2)
    plt.grid(which='both', alpha=.4)
    plt.xlabel('N')
    if display_title:
        plt.title('Comparison with coalescent models')
    if p_max == 3:
        plt.ylabel(r'$\langle T_3\rangle/\langle T_2\rangle$')
    else:
        plt.ylabel(r'$\langle T_p\rangle/\langle T_2\rangle$ for several values of $p$')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show()
    
    
def several_curves_Tp(Ns, T2_path, data_path, suffixes, labels, 
                      title=None):
    """
    Plot the Tp curves for 
    """
    colors = ['b', 'c', 'g', 'y', 'r']
    plt.semilogx(Ns, [2-2/3]*len(Ns), 'b--', alpha=1, label='Kingman')
    plt.semilogx(Ns, [get_expected_Tp(3),]*len(Ns), 'r-.', label='Bolthausen-Sznitman')
    for i, suffix in tqdm(enumerate(suffixes)):
        avgs = [get_or_write_avg(T2_path, data_path, suffix, N, 3) \
               for N in Ns]
        plt.semilogx(Ns, avgs, f'{colors[i]}+', label=labels[i])
        plt.semilogx(Ns, avgs, f'{colors[i]}-', alpha=.4)
    plt.grid(which='both', alpha=.4)
    plt.xlabel('N')
    if title is not None:
        plt.title(title)
    plt.ylabel(r'$\langle T_3\rangle/\langle T_2\rangle$')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show()
    
    
    
def plot_random_fluctuations(T2_path, data_path, suffix, N, p=3, lim=99999):
    """
    Curve showing all simulations independantly, showing the importance of
    doing several sims.
    """
    datetimes = os.listdir(join(T2_path, suffix, str(N)))
    if len(datetimes) >= lim:
        datetimes = datetimes[:lim]
    gens = range(3000, 8000+1)
    regex_gens = re.compile(r'gen(\d+)_')
    for datetime in tqdm(datetimes):
        filename = f'{N}{suffix}_{datetime}_{str(gens)}'
        if exists(join(data_path, filename)):
            with open(join(data_path, filename), 'r') as f:
                to_plot = eval(f.read())
        else:
            to_plot = [None for _ in gens]
            available_gens = set([regex_gens.search(file).group(1) for\
                              file in os.listdir(join(T2_path, suffix, str(N), datetime))\
                                  if int(regex_gens.search(file).group(1)) in gens])
            for gen in tqdm(available_gens):
                T2vals, Tpvals = [], []
                for T2file in os.listdir(join(T2_path, suffix, str(N), datetime)):
                    if str(gen) in T2file:
                        with open(join(T2_path, suffix, str(N), datetime, T2file), 'r') as f:
                            T2 = json.load(f)
                        T2vals.extend(T2.values())
                        Tpvals.extend(get_Tp_from_T2(T2, p).values())
                if len(T2vals) == 0:
                    raise ValueError (f'No file detected for T2 for {N}{suffix}gen{gen} in {datetime}')
                else:
                    to_plot[int(gen) - gens[0]] = (np.mean(Tpvals)/np.mean(T2vals))
            with  open(join(data_path, filename), 'w') as f:
                f.write(str(to_plot))
        plt.plot(gens, to_plot)
    plt.plot(gens, [4/3,]*len(gens), 'g--', label='Kingman')
    plt.plot(gens, [5/4,]*len(gens), 'g-.', label='Bolthausen')
    plt.xlabel('generations')
    plt.ylabel(r'$\langle T_3\rangle/\langle T_2\rangle$')
    plt.title(f'Influence of random flucutations for model A (haploid) with N={N}')
    plt.grid(which='both', alpha=.4)
    plt.legend()
    plt.show()
    
    
def plot_correlation(T2_path, data_path, suffix, N,
                     gaps=range(1, 1001),#[1, 2, 5, 10, 20, 50, 80, 100, 200, 500, 700, 1000], 
                     lim=999999):
    """
    Plot the autocorrelation values of the average T3 divided by the average T2
    for different gaps.
    """
    datetimes = os.listdir(join(T2_path, suffix, str(N)))
    regex_gens = re.compile(r'gen(\d+)_')
    plt.semilogx(gaps, [0 for _ in gaps], 'r--', alpha=.2)
    if len(datetimes) >= lim:
        datetimes = datetimes[:lim]    
    for datetime in datetimes:
        gens = set([regex_gens.search(file).group(1) for \
                    file in os.listdir(join(T2_path, suffix, str(N), 
                                            datetime))])
        gens = sorted(list([int(gen) for gen in gens]))
        filename = f'{N}{suffix}_{datetime}_range({min(gens)}, {max(gens)})'
        if exists(join(data_path, filename)):
            with open(join(data_path, filename), 'r') as f:
                avgs = eval(f.read())
        else:
            avgs = [None for _ in gens]
            for gen in tqdm(gens):
                T2vals, Tpvals = [], []
                for T2file in os.listdir(join(T2_path, suffix, str(N), datetime)):
                    if str(gen) in T2file:
                        with open(join(T2_path, suffix, str(N), datetime, T2file), 'r') as f:
                            T2 = json.load(f)
                        T2vals.extend(T2.values())
                        Tpvals.extend(get_Tp_from_T2(T2, 3).values())
                if len(T2vals) == 0:
                    raise ValueError (f'No file detected for T2 for {N}{suffix}gen{gen} in {datetime}')
                else:
                    avgs[int(gen) - gens[0]] = (np.mean(Tpvals)/np.mean(T2vals))
            with open(join(data_path, filename), 'w') as f:
                f.write(str(avgs))
        to_plot = [np.corrcoef(avgs[:-gap], avgs[gap:])[0, 1] for gap in gaps]
        plt.semilogx(gaps, to_plot)
    plt.xlabel('lag (number of generations between recorded trees)')
    plt.ylabel(r'Autocorrelation between values of $\langle T_3\rangle/\langle T_2\rangle$')
    plt.title(f'Autocorrelation between subsequent trees for model A with\nindividual fitnessfor several simulations with N={N}')
    plt.grid(which='both', alpha=.4)
    plt.legend()
    plt.show()
    
    
def plot_SE(T2_path, data_path, suffix, pmax):
    """
    Curve of the standard error of T3/T2 with respect to N.
    """
    colors = ['b', 'r', 'y', 'k', 'm']
    Ns = find_Ns(T2_path, suffix, verbose=False)
    for p in range(3, pmax+1):
        to_plot = [get_or_write_SE(T2_path, data_path, suffix, N, p) for N in Ns]
        plt.loglog(Ns, to_plot, f'{colors[p-3]}+', label=fr'$\langle T_{p}\rangle/\langle T_2\rangle$')
        plt.loglog(Ns, to_plot, f'{colors[p-3]}-', alpha=.2)
    plt.grid(which='both', alpha=.4)
    plt.xlabel('N')
    plt.ylabel('SEM')
    plt.legend()
    plt.title('Standard error of the mean')
    plt.tight_layout()
    plt.show()
    
    
def plot_diff_bolt(T2_path, data_path, models, Ns, p=3):
    """
    Curve of the squared difference between empirical values and theoretical 
    values for Bolthausen-Sznitman, for all models.
    """
    dict_names = {'H':'Model A (haploid)',
                  '_':'Model A_1 (diploid)',
                  'A2':'Model A_2 (diploid)',
                  'EH':'Exponential model (haploid)',
                  'E':'Exponential model (diploid)',
                  'M':'Mutations model (diploid)',
                  'MH':'Mutations model (haploid)',
                  '_MIN':'Model A_1 with min',
                  '_MAX':'Model A_1 with max',
                  'A2MIN':'Model A_2 with min',
                  'A2MAX':"Model A_2 with max",
                  'X':'Model A_1 but hermaphroditic',
                  'EMIN':'Exponential model with min',
                  'EMAX':'Exponential model with max',
                  'A2M':'Model A_2 with males selection',
                  'A2MX':'Model A_2 with herm. selec.'}
    dict_colors = {'H':'c-',
                   '_':'y-',
                   'A2':'y--',
                   'EH':'c+-',
                   'E':'y+-',
                   'M':'yx-',
                   'MH':'cx-',
                   '_MIN':'b-',
                   'A2MIN':'b--',
                   '_MAX':'r-',
                   'A2MAX':'r--',
                   'X':'y:',
                   'EMAX':'r+-',
                   'EMIN':'b+-',
                   'A2M':'k--',
                   'A2MX':'k-.',
                   }
    diff_w_bolt = lambda x: (x - get_expected_Tp(p))**2
    for suffix in models: 
        to_plot = [diff_w_bolt(get_or_write_avg(T2_path, data_path, suffix, N, p))\
                   for N in Ns]
        plt.loglog(Ns, to_plot, dict_colors[suffix], label=dict_names[suffix])
    plt.xlabel('N')
    plt.ylabel('Squared difference between empirical $\\langle T_3\\rangle/\\langle T_2\\rangle$\n and theoretical value of the Bolthausen-Sznitman')
    plt.title('Convergence speed of several models towards the Bolthausen')
    plt.legend()
    plt.grid(which='both', alpha=.4)
    plt.tight_layout()
    plt.show()
    
    
def plot_unimodal(T2_path, data_path, suffix, N, d=.1):
    """
    Check if T3/T2 is unimodal (doesn't seem to be)
    """
    avgs = [get_or_write_avg(T2_path, data_path, suffix, N, 3, 
                             timestamp=timestamp) for timestamp in \
            tqdm(os.listdir(join(T2_path, suffix, str(N))))]
    bins = np.arange(min(avgs), max(avgs)+d, d)
    plt.hist(avgs, bins=bins, edgecolor='black', align='left', rwidth=0.8)
    plt.xlabel(f'Intervals (length = {d})')
    plt.ylabel(r'Number of $\langle T_3\rangle/\langle T_2\rangle$ in the interval')
    plt.title('Density estimate of simulation results')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.gca().set_xticklabels([])
    plt.xticks(bins)
    plt.show()


def main(trees_path, data_path, p_max=4, plot_type=1): #_ H E EH M MH
    positions = [(0,-50), (650,-50), (1300,-50), (0,500), (650,500), 
                 (1300,500), (250, 250), (1000, 250),
                 (0,-50), (650,-50), (1300,-50), (0,500), (650,500), 
                              (1300,500), (250, 250), (1000, 250),
                (0,-50), (650,-50), (1300,-50), (0,500), (650,500), 
                             (1300,500), (250, 250), (1000, 250)]
    i = 0
    suffixes = {'A, haploid':'H', 
                '$A_1$ (fitness on individuals) diploid':'_',
                '$A_2$ (fitness on chromosomes), diploid':'A2',
                'Exponential, haploid':'EH',
                'Exponential, diploid':'E',
                'A_2max (fitness on chromosomes, dominant fitness), diploid':'A2MAX',
                'Mutations model (haploid)':'MH',
                'Mutations model (diploid)':'M',
              }
    ('_', 'E', 'A2', 'A2MAX', '_MAX', 'EMAX', 'A3MAX')
    suffixes2 = {'A with max':'_MAX',
               'A with chromosomes and max':'A2MAX',
               'Exponential with max':'EMAX',
               'A with chromosomes and superfit':'A3',
               'A with chromosomes and selection of male based on fitness':'A2M',
               'A with chromosomes, max and selection of males based on fitness':'A2MMAX'}
    suffixes3 = {'Model A1 with min, diploid':'_MIN',
                 'Model A2 with min, diploid':'A2MIN',
                 'Model A but asexual':'X',
                 'Model A but asexual and with max':'XMAX',
                 'Model A but with chromosomes, asexual and selection based on fitness':'A2MX',}
    suffixes4 = {'Model A, haploid':'H',
                 r'Model A_1 (diploid)':'_',
                 r'Model A_2 (diploid)':'A2',
                 'Exponential model, haploid':'EH',
                 'Exponential model, diploid':'E',
                 'Mutations model, haploid':'MH',
                 'Mutations model, diploid':'M',
                 r'Model A_1 with max, diploid':'_MAX',
                 r'Model A_2 with max, diploid':'A2MAX',
                 r'Model A_1 with min, diploid':'_MIN',
                 r'Model A_2 with min, diploid':'A2MIN',
                 r'Model A_1 but hermaphroditic':'X',
                 #r'Model $\rm A_2$ but hermaphroditic':'A2X',
                 #r'Model $\rm A_2$ with a small chance of huge fitness gain':'A3',
                 'Model A_2 with selection of males based on fitness':'A2M',
                 'Model A_2 but hermaphrotitic and only fittest can reproduce':'A2MX',
                 'Exponential model with max':'EMAX',
                 'Exponential model with min':'EMIN',
                 }
    if plot_type == 1:
        for name, suffix in tqdm(suffixes.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=suffix)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=4, colspan=3)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                plt.subplot2grid((4, 6), (1, 3), rowspan=4, colspan=3)
                plot_SE(T2_path, data_path, suffix, 4)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
    elif plot_type == 2:
        #N = 1000
        #for name, suffix in tqdm(suffixes.items()):
        #    fig = plt.figure(num=name)
        #    fig.canvas.manager.window.move(*positions[i])
        #    try:
        plot_random_fluctuations(T2_path, data_path, '_', 1000, lim=4)
        #    except FileNotFoundError:
        #        print(f'Skipped {suffix} because no files yet for N={N}.')
        #    i += 1
    elif plot_type == 3:
        for name, suffix in tqdm(suffixes2.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=suffix)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=4, colspan=3)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                plt.subplot2grid((4, 6), (1, 3), rowspan=4, colspan=3)
                plot_SE(T2_path, data_path, suffix, 4)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
    elif plot_type == 6:
        for name, suffix in tqdm(suffixes3.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=suffix)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=4, colspan=3)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                plt.subplot2grid((4, 6), (1, 3), rowspan=4, colspan=3)
                plot_SE(T2_path, data_path, suffix, 4)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
    if plot_type == 7:
        for name, suffix in tqdm(suffixes4.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=suffix)
            #fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=4, colspan=6)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                #plt.subplot2grid((4, 6), (1, 3), rowspan=4, colspan=3)
                #plot_SE(T2_path, data_path, suffix, 4)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
        plt.figure(num=f'Table for p={p_max}')
        plot_diff_bolt(T2_path, data_path, ['H', '_', 'A2', 'EH', 'E', 
                                            '_MAX', 'A2MAX', '_MIN',
                                            'A2MIN', 'X', 'EMIN', 'EMAX',
                                            'A2M', 'A2MX'], 
                       [100, 316, 1000, 3160, 10000])
    elif plot_type == 8:
        plot_unimodal(T2_path, data_path, '_', 100, d=.005)
    elif plot_type == 9:
        #Dominance coeff for Model A1
        several_curves_Tp([100, 316, 1000, 3160, 10000], 
                          T2_path, data_path, 
                          ['_MIN', '_ARGS0.3_4', '_', '_ARGS0.7_4', '_MAX'],
                          ['s=0', 's=0.3', 's=0.5', 's=0.7', 's=1'],
                          'Comparison of Model A_1 with coalescent models, \n for different values of the dominance coefficient s of mutations')
    elif plot_type == 10:
        #Nb children for A1
        several_curves_Tp([100, 316, 1000, 3160, 10000], 
                          T2_path, data_path, 
                          ['_', '_ARGS0.5_8', '_ARGS0.5_16', '_ARGS0.5_32'],
                          ['4 children per female', '8 children per female',
                           '16 children per female', '32 children per female'],
                          'Comparison of Model A_1 with coalescent models,\nfor different numbers of potential children per female')
    elif plot_type == 11:
        plot_correlation(T2_path, data_path, '_', 1000, lim=99999)
    elif plot_type == 12:
        several_curves_Tp([100, 316, 1000, 3160, 10000, 31600], 
                          T2_path, data_path, 
                          ['H', '_', 'A2'],
                          ['Haploïde', 'Modèle 1 (diploïde)', 'Modèle 2 (diploïde)'],
                          'Premiers résultats de convergence')
    elif plot_type == 13:
        several_curves_Tp([100, 316, 1000, 3160, 10000], 
                          T2_path, data_path, 
                          ['_', 'X'], ['Modèle 1 (diploïde)',
                                            'Modèle 1, mais hermaphrodite'])
    elif plot_type == 0:
        regex = re.compile(r"Brunet2 \d+([A-Z0-9_\.]+)_T\d_(SE|_avg)")
        for file in os.listdir(data_path):
            m = regex.search(file)
            if m is not None:
                os.remove(join(data_path, file))


if __name__ == '__main__':
    T2_path = r"C:\Users\trist\Desktop\Trees5"
    data_path = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard\Stage M1\Code\plots_data"
    plt.close('all')
    main(T2_path, data_path, p_max=4, plot_type=7) #1,3,6,7
