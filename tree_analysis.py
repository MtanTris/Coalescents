# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 11:00:41 2025

@author: trist
"""

import tskit
from IPython.display import display
import random as rd
from scipy.special import binom
from itertools import combinations
from tqdm import tqdm

def sample_tree(tree, p=None, sample=None):
    """
    Sample p leaves of a tree randomly or with a given sample.

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.
    p : int, optional
        Number of leaves to sample. The default is None.
    sample : tuple of ints, optional
        Contains the id of leaves that will be sampled. The default is None.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    N = tree.num_samples
    if sample is None:
        if p is not None:
            try:
                sample = rd.sample(range(N), p)
            except ValueError:
                #plot whole tree and raise error
                svg_size = (800, 400)
                svg_string = tree.draw_svg(size=svg_size)
                display(svg_string)
                raise ValueError (f'Sample larger than population size or is negative. ({p}, {N})')
        else:
            raise ValueError ('p or samples must be specified')
    return tree.simplify(samples=sample)


def display_sample(tree, p=None, sample=None, cut_haploid=False):
    sampled_tree = sample_tree(tree, p=p, sample=sample, cut_haploid=cut_haploid)
    svg_size = (800, 400)
    svg_string = sampled_tree.draw_svg(size=svg_size)
    display(svg_string)
    
    
def get_MRCA_time_4sampled(tree, sample):
    """
    Get the MRCA time of the sample of a tree.

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.
    sample : tuple of ints
        Set of leaves to sample.

    Returns
    -------
    float
        MRCA time of sampled tree.

    """
    sampled = sample_tree(tree, sample=sample).first()
    return sampled.time(sampled.mrca(*sampled.samples()))



def get_T2(tree, disable=True):
    """
    Extract the T2 data (MRCA of every couple of leaves) of a tree.

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.
    disable : bool, optional
        Choose if the tqdm.tqdm loading bar should be displayed. 
        The default is True.

    Returns
    -------
    T2 : dict (tuple of ints):float
        MRCA of every ordered couple of leaves.

    """
    if not disable:
        print('Computing T2 matrix...')
    #Get leaves (if more than one root exists, the error will be picked up
    #by simulate3.run_sims()).
    try:
        leaves = tuple(tree.first().leaves(tree.first().root))
    except ValueError: #more than one root exists
        return None
    N = len(leaves)
    
    #Iterate over every ordered couple of leaves
    couples = combinations(leaves, 2)
    T2 = {}
    for couple in tqdm(couples, total=binom(N, 2), disable=disable):
        CA_time = get_MRCA_time_4sampled(tree, couple)
        T2[str(couple)] = CA_time
    return T2


def get_Tp_from_T2(T2, p, disable=True):
    """
    Extract the Tp data (MRCA time of p-sample) of a tree, based on T2.

    Parameters
    ----------
    T2 : dict (tuple of ints):float
        Contains MRCA of every couple of individuals in a sampled tree.
    p : int
        Size of p-sample.
    disable : bool, optional
        Choose if the tqdm.tqdm loading bar should be displayed. 
        The default is True.

    Raises
    ------
    ValueError
        Tp is empty at the end of the function.
        Likely due to bad T2.

    Returns
    -------
    Tp : dict
        Contains MRCA of every p-sample of individuals in a sampled tree.

    """
    if not disable:
        print(f'computing T{p}...')
    #Sort T2 keys (necessary to not raise KeyErrors later)
    T2 = {tuple(sorted(eval(k))):v for k,v in T2.items()}
    
    Tp = {}
    #Get leaves of tree from T2 keys
    leaves = tuple(set().union(*[X for X in T2.keys()]))
    
    #Iterate over every ordered p-sample of leaves
    p_tuples = combinations(leaves, p)
    for p_tuple in tqdm(p_tuples, disable=disable):
        
        #Tp value for a given p-sample key obtained by a max of T2 values
        #(process used in "Noisy travelling waves" by Brunet et al.)
        Tp[str(p_tuple)] = max([T2[(p_tuple[0], leaf)] for leaf in p_tuple[1:]])
    
    #Raise error if Tp is empty (to not write an empty file later)
    if Tp == {}:
        raise ValueError ('Tp is empty.')
    return Tp
    
if __name__ == '__main__':
    semipath = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard" +\
                r"\Stage M1\Code\Trees\Brunet"
    path = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard" +\
            r"\Stage M1\Code\Trees\Brunet10000H.trees"
    #ts = tskit.load(path)
    #ts = ts.simplify()
    #tree_path = r"C:\Users\trist\OneDrive\Notes Cours\🔢 MATHÉMATIQUES\M1 Hadamard\Stage M1\Code\Trees\Brunet10H.trees"
    tree_path = r"C:\Users\trist\Desktop\Trees3\Brunet2 100Hgen1141 2025-05-28 17_26_37.trees"
    
    ts = tskit.load(tree_path)
    ts = ts.simplify()
    #tree = ts.first(); connected = set(); [connected.update(tree.nodes(root)) for root in tree.roots]; ts = ts.simplify([u for u in ts.samples() if u in connected])

    tree = ts.first()
    root = max(tree.roots) if tree.roots else None
    #ts = ts.samples(tree.children(root))
    display_sample(ts, sample=list(tree.samples(root)))
    print(root)
    
    #compute_fractions_T(ts, pre_sample_size=200)
    #get_curves_p(ts, 5, pre_sample_size=50)
    #Ns = [ 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 
    #      1400, 1500, 1600, 1800, 2500]
    #get_curves_N(semipath, Ns, sample_size=200, curve=1)
    #get_curves_p(ts, 5)
    #ts = prune_tree_from_root(ts)
    #display_sample(ts, 100)
    
    #L'ARBRE EST EN DOUBLE PARCE QUE HAPLOIDE !!!!!!!!!
    #DANS LE CAS DIPLOIDE, C'EST PARCE QUE CHAQUE CHROMOSOME EST PISTE