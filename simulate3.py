# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 11:03:10 2025

@author: trist
"""


import re, subprocess, os, json, tskit
from datetime import datetime
from os.path import join, isdir
from tqdm import tqdm, trange
from time import sleep
import random as rd

from tree_analysis import get_T2, sample_tree


def modify_pop_size(slim_path, N):
    """
    Modify the population size constant of a SLiM file.

    Parameters
    ----------
    slim_path : str
        Path of the SLiM file.
    N : int
        New population size.

    Returns
    -------
    None.

    """
    regex_N = re.compile(r"defineConstant\('N', \d+\);")
    with open(slim_path, 'r') as f:
        text = f.read()
    text = regex_N.sub(f"defineConstant('N', {N});", text)
    with open(slim_path, 'w') as f:
        f.write(text)


def run_sims(slim_path, cwd, T2_dest, suffix, N, nb_sims=20, n_samples=10,
             sample_size=10):
    """
    Run SLiM simulations.

    Parameters
    ----------
    slim_path : str
        Path of the SLiM file.
    cwd : str
        Path of the folder where the trees will be recorded by SLiM.
        (should be different from T2_dest to prevent bugs)
    T2_dest : str
        Path of the folder where the T2 dict will be stored.
    suffix : str
        Small string with only letters and numbers specifying which model is
        simulated.
    N : int
        Size of the population in the simulations.
    nb_sims : int, optional
        Number of different simulations that will be executed.
        Basically this function is a big for loop. The default is 20.
    n_samples : int, optional
        Number of T2 dictionnaries that will be created at each generations. 
        The samples will be disjointed and chosen randomly. The default is 10.
    sample_size : int, optional
        Number of individuals in each sample. Allows for more precise averages
        but increase exponentially the number of computations necessary to get
        the Tp dictionnaries. The default is 10.

    Returns
    -------
    None.

    """
    #Set up regex to get the generation number of a tree file and change 
    #pop size in the SLiM file
    regex = re.compile(r'gen(\d+)')
    modify_pop_size(slim_path, N)
    
    for _ in trange(nb_sims):
        #Run a SLiM simulation and record trees
        subprocess.run([r"C:\msys64\mingw64\bin\slim.exe", slim_path], cwd=cwd)
        
        #Record date and time for naming T2 files later
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S').replace(':', '_')
        
        #Iterate over every generated tree
        files = os.listdir(cwd)
        for file in files:
            
            #Check if tree has been created by this instance of run_sims()
            #(not always the case with parallel computing like in ALICE)
            if f'Brunet2 {N}{suffix}gen' in file and str(N) in file:
                
                #Get generation number
                m = regex.search(file)
                if m is not None:
                    gen = int(m.group(1))
                    
                #Load tree
                ts = tskit.load(join(cwd, file))
                ts = ts.simplify()
                
                #If model is haploid, remove all disconnected leaves
                if 'H' in suffix:
                    tree = ts.first()
                    root = max(tree.roots) if tree.roots else None
                    ts = sample_tree(ts, sample=list(tree.samples(root)))
                
                #Get the disjointed samples and iterate over them
                pre_sample = rd.sample(range(N), n_samples*sample_size)
                samples = [pre_sample[i*sample_size:(i+1)*sample_size] for i in range(n_samples)]
                for i, sample in enumerate(samples):
                    
                    #Sample the tree and get T2 dict of sampled tree
                    if sample_size is not None:
                        ts2 = sample_tree(ts, sample=sample)
                    T2 = get_T2(ts2, disable=True)
                    if T2 is None:
                        print('Error : more than one root exists')
                        continue
                    
                    #Create all necessary folders (to not raise an error later)
                    for dir in (join(T2_dest, suffix), 
                                join(T2_dest, suffix, str(N)),
                                join(T2_dest, suffix, str(N), timestamp)):
                        if not isdir(dir):
                            os.mkdir(dir)
                            
                    #Write the T2 dict as a json file
                    with open(join(T2_dest, suffix, str(N), timestamp, 
                                   f'T2_gen{gen}_sampl{i}.json'), 'w') as f:
                        json.dump(T2, f)
                
                #Delete the tree
                os.remove(join(cwd, file))
                
                #Wait one second to be sure that the timestamp will change 
                #for the next simulation
                sleep(1)
        

def main(suffix, N, nb_sims):
    #slim_path = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard\Stage M1\Code\SLiM\Brunet2_"
    slim_path = r"C:\Users\trist\Desktop\Brunet2_"
    cwd = r"C:\Users\trist\Desktop\temp"
    T2_dest = r"C:\Users\trist\Desktop\Trees4"
    slim_path += suffix
    run_sims(slim_path, cwd, T2_dest, suffix, N, nb_sims=nb_sims, n_samples=10,
             sample_size=10)


if __name__ == '__main__':
    for suffix in tqdm(('AA2', 'AAA2', 'AAAA2', 'AAAAA2', 'AAAAAA2')):
    #for suffix in tqdm(('A2MMAX', 'A2MAX', '_MAX', 'A3MAX', 'EMAX')): #A2M
    #for suffix in tqdm(('_', 'H', 'E', 'EH', 'M', 'MH')):
        print(suffix)
        main(suffix, 100, 15)
        main(suffix, 316, 15)
    
    
    