# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 13:47:02 2025

@author: trist
"""

import os, json
from os.path import join, exists
import numpy as np
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
        plt.show()
    return Ns_list


def get_or_write_avg(T2_path, data_path, suffix, N, p):
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
        DESCRIPTION.

    Returns
    -------
    avg : TYPE
        DESCRIPTION.

    """
    filename = f'Brunet2 {N}{suffix}_T{p}'
    if exists(join(data_path, filename)):
        with open(join(data_path, filename), 'r') as f:
            avg = float(f.read())
    else:
        print(f'Computing average T{p} for {N}{suffix}...')
        values = []
        for datetime in os.listdir(join(T2_path, suffix, str(N))):
            for T2file in os.listdir(join(T2_path, suffix, str(N), datetime)):
                with open(join(T2_path, suffix, str(N), datetime, T2file), 'r') as f:
                    T2 = json.load(f)
                if p == 2:
                    values.extend(T2.values())
                else:
                    Tp = get_Tp_from_T2(T2, p)
                    values.extend(Tp.values())
        if len(values) == 0:
            raise ValueError('values is empty. No files detected.')
        avg = np.mean(values)
        if p != 2:
            avg /= get_or_write_avg(T2_path, data_path, suffix, N, 2)
        with open(join(data_path, filename), 'w') as f:
            f.write(str(avg))
    return avg
            

def curve_T2(Ns, trees_path, data_path, suffix):
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
    
    
    
def curves_Tp(Ns, T2_path, data_path, suffix, p_max=5):
    colors = ['b', 'r', 'y', 'k', 'm']
    for p in range(3, p_max+1):
        avgs = []
        for N in Ns:
            avgs.append(get_or_write_avg(T2_path, data_path, suffix, N, p))
        plt.semilogx(Ns, [2-2/p]*len(Ns), f'{colors[p-3]}--',
                     label='Kingman' if p==3 else '', alpha=1)
        plt.semilogx(Ns, [get_expected_Tp(p),]*len(Ns), color=colors[p-3],
                     label='Bolthausen' if p==3 else '')
        plt.semilogx(Ns, avgs, f'{colors[p-3]}x', 
                     label='Empirical values')
        plt.semilogx(Ns, avgs, f'{colors[p-3]}-', alpha=.2)
    plt.grid(which='both', alpha=.4)
    plt.xlabel('N')
    plt.title('Analysis of the model')
    plt.ylabel(r'$\langle T_p\rangle/\langle T_2\rangle$ for several values of $p$')
    plt.legend()
    plt.show()
    
    
def curve_dominance(Ns, T2_path, data_path):
    colors = {'MA':'b', 'MB':'r', 'MC':'y', 'MD':'k', 'ME':'m'}
    coeffs = {'MA':.1, 'MB':.3, 'MC':.5, 'MD':.8, 'ME':1.}
    for suffix, color in tqdm(colors.items()):
        avgs = []
        for N in Ns:
            avgs.append(get_or_write_avg(T2_path, data_path, suffix, N, 3))
        plt.semilogx(Ns, avgs, f'{colors[suffix]}x', 
                     label=f'dom. coeff. = {coeffs[suffix]}')
        plt.semilogx(Ns, avgs, f'{colors[suffix]}-', alpha=.2)
    plt.semilogx(Ns, [4/3,]*len(Ns), 'g--', label='Kingman', alpha=.2)
    plt.semilogx(Ns, [5/4,]*len(Ns), 'g-.', label='Bolthausen')
    plt.grid(which='both', alpha=.4)
    plt.xlabel('N')
    plt.title('Importance of the dominance coefficicient of mutations, for N=100')
    plt.ylabel(r'$\langle T_3\rangle/\langle T_2\rangle$')
    plt.legend()
    plt.show()
    
    
def plot_random_fluctuations(T2_path, data_path, suffix, N, p=3):
    datetimes = os.listdir(join(T2_path, suffix, str(N)))
    gens = range(10*N, 10*N+10+1)
    for datetime in datetimes:
        if exists(join(data_path, f'{N}{suffix}_{datetime}')):
            with open(join(data_path, f'{N}{suffix}_{datetime}'), 'r') as f:
                to_plot = eval(f.read())
        else:
            to_plot = []
            for gen in gens:
                T2vals, Tpvals = [], []
                for T2file in os.listdir(join(T2_path, suffix, str(N), datetime)):
                    if str(gen) in T2file:
                        with open(join(T2_path, suffix, str(N), datetime, T2file), 'r') as f:
                            T2 = json.load(f)
                        T2vals.extend(T2.values())
                        Tpvals.extend(get_Tp_from_T2(T2, p).values())
                if len(T2vals) == 0:
                    to_plot.append(None)
                    #raise ValueError (f'No file detected for T2 for {N}{suffix}gen{gen} in {datetime}')
                else:
                    to_plot.append(np.mean(Tpvals)/np.mean(T2vals))
            with  open(join(data_path, f'{N}{suffix}_{datetime}'), 'w') as f:
                f.write(str(to_plot))
        plt.plot(gens, to_plot)
    plt.plot(gens, [4/3,]*len(gens), 'g--', label='Kingman')
    plt.plot(gens, [5/4,]*len(gens), 'g-.', label='Bolthausen')
    plt.xlabel('generations')
    plt.ylabel(r'$\langle T_3\rangle/\langle T_2\rangle$')
    plt.title(f'Influence of random flucutations for model {suffix} with N={N}')
    plt.grid(which='both', alpha=.4)
    plt.legend()
    plt.show()
    
    
def plot_nb_children(T2_path, data_path, Ns, p=3):
    suffixes = ('A2', 'AA2', 'AAA2', 'AAAA2', 'AAAAA2', 'AAAAAA2')
    nb_childrenF = (2, 4, 8, 16, 32, 64)
    colors = ['b', 'r', 'y', 'k', 'm']
    for i, N in enumerate(Ns):
        to_plot = [get_or_write_avg(T2_path, data_path, suffix, N, p) \
                   for suffix in suffixes]
        plt.semilogx(nb_childrenF, to_plot, f'{colors[i]}+', label=f'N={N}')
        plt.semilogx(nb_childrenF, to_plot, f'{colors[i]}-', alpha=.2)
    plt.plot(nb_childrenF, [4/3,]*len(suffixes), 'g--', label='Kingman')
    plt.plot(nb_childrenF, [5/4,]*len(suffixes), 'g-.', label='Bolthausen')
    plt.xlabel('Number $k$ of children per female')
    plt.ylabel(rf'$\langle T_{p}\rangle/\langle T_2\rangle$')
    plt.title('Influence of the number of children on Model A')
    plt.grid(which='both', alpha=.4)
    plt.legend()
    plt.show()
    
        
    
def main(trees_path, data_path, p_max=4, plot_type=1): #_ H E EH M MH
    positions = [(0,-50), (650,-50), (1300,-50), (0,500), (650,500), 
                 (1300,500), (250, 250), (1000, 250)]
    i = 0
    suffixes = {'A, haploid':'H', 'Exponential, haploid':'EH', 
              'Mutations, haploid':'MH', 
              'A_1 (fitness on individuals) diploid':'_', 
              'Exponential, diploid':'E', 'Mutations, diploid':'M',
              'A_2 (fitness on chromosomes), diploid':'A2',
              'A_2max (fitness on chromosomes, dominant fitness), diploid':'A2MAX',
              }
    ('_', 'E', 'A2', 'A2MAX', '_MAX', 'EMAX', 'A3MAX')
    suffixes2 = {'A with chromosomes':'A2',
               'A with max':'_MAX',
               'A with chromosomes and max':'A2MAX',
               'Exponential with max':'EMAX',
               'A with chromosomes, max and superfit':'A3MAX',
               'A with chromosomes and selection of male based on fitness':'A2M',
               'A with chromosomes, max and selection of males based on fitness':'A2MMAX'}
    if plot_type == 1:
        for name, suffix in tqdm(suffixes.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=name)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=3, colspan=3)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                plt.subplot2grid((4, 6), (1, 3), rowspan=3, colspan=3)
                curve_T2(Ns, trees_path, data_path, suffix)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
    elif plot_type == 2:
        N = 316
        for name, suffix in tqdm(suffixes.items()):
            fig = plt.figure(num=name)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plot_random_fluctuations(T2_path, data_path, suffix, N)
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet for N={N}.')
            i += 1
    elif plot_type == 3:
        for name, suffix in tqdm(suffixes2.items()):
            #plt.figure(figsize=(2*plt.rcParams["figure.figsize"][0], 3*plt.rcParams["figure.figsize"][1]))
            fig = plt.figure(num=name)
            fig.canvas.manager.window.move(*positions[i])
            try:
                plt.subplot2grid((4, 6), (0, 0), rowspan=1, colspan=6)
                Ns = find_Ns(trees_path, suffix=suffix, verbose=True)
                print(suffix, ':', Ns)
                plt.subplot2grid((4, 6), (1, 0), rowspan=3, colspan=3)
                curves_Tp(Ns, trees_path, data_path, suffix, p_max=p_max)
                plt.subplot2grid((4, 6), (1, 3), rowspan=3, colspan=3)
                curve_T2(Ns, trees_path, data_path, suffix)
                plt.suptitle(f'model : {name}')
                plt.tight_layout()
            except FileNotFoundError:
                print(f'Skipped {suffix} because no files yet.')
            i += 1
    elif plot_type == 4:
        Ns = [100, 316]
        curve_dominance(Ns, trees_path, data_path)
    elif plot_type == 5:
        Ns = [100, 316]
        plot_nb_children(T2_path, data_path, Ns)
 


if __name__ == '__main__':
    T2_path = r"C:\Users\trist\Desktop\Trees4"
    data_path = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard\Stage M1\Code\plots_data"
    plt.close('all')
    main(T2_path, data_path, p_max=4, plot_type=3)