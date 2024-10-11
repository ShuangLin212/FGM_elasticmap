"""Convienence functions for building dilute structures"""
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from itertools import permutations, product, chain
from collections import Counter
import copy

def dilute_substitution(endmembers, sublattice_dict, supercell_matrix=None):
    """
    Return dilute structures base on endmembers and elements

    Parameters
    ----------
    endmembers: list of endmember structures from endmember.py

    sublattice_dict: dict of elements for substitutions for each sublattice, e.g., {'a': ['Hf', 'Mo'], 'f': ['Hf', 'Mo'], 'h': ['Hf', 'Mo']}
    
    supercell_matrix: list of matrix to make supercell, e.g., [2,2,2], [[2,1,0],[0,3,0],[0,0,1]], default is none

    Returns
    -------
    dilute structures: list
    
    dilute sublattice configurations
    """
    dilute_strs=[]
    dilute_conf=[]
    structures=copy.deepcopy(endmembers)
    if supercell_matrix is not None:
        for i in structures:
            i.make_supercell(supercell_matrix)
    for strs in structures:
        dictstr=strs.as_dict()
        sub_lattice_list=[]
        for i in dictstr['sites']:
            j=i['properties']['sublattice_sites']
            sub_lattice_list.append(j)
        num=[]
        for i in range(0, len(sub_lattice_list)):
            if i != len(sub_lattice_list):
                if sub_lattice_list[i] != sub_lattice_list[i-1]:
                    num.append(i)          
        for i in num:
            dictstr=strs.as_dict()
            old_ele=dictstr['sites'][i]['species'][0]['element']
            site_sub=dictstr['sites'][i]['properties']['sublattice_sites']
            sublattice_num=sub_lattice_list[i]
            element_dict=sublattice_dict[sublattice_num]
            for j in element_dict:
                if j != 'fix':
                    if old_ele != j:
                        dictstr['sites'][i]['species'][0]['element']=j
                        dilute_strs.append(Structure.from_dict(dictstr))
                                      
    for dilute in dilute_strs:
        comb={}
        dictdilute=dilute.as_dict()
        for i in num:
            sfirst_ele=dictdilute['sites'][i]['species'][0]['element']
            ssecond_ele=dictdilute['sites'][i+1]['species'][0]['element']
            site_sub=dictdilute['sites'][i]['properties']['sublattice_sites']
            if sfirst_ele == ssecond_ele:
                comb[site_sub]=sfirst_ele
            else:
                conflist=list()
                conflist.append(sfirst_ele)
                conflist.append(ssecond_ele)
                comb[site_sub]=conflist
        dilute_conf.append(comb)
    return dilute_strs, dilute_conf
