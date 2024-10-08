
"""
Convenience functions for building dilute structures
"""

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from itertools import permutations, product, chain
from collections import Counter
import copy

def dilute_substitution(endmembers, sublattice_dict, supercell_matrix=None):
    """
    Return dilute structures based on endmembers and elements
    """
    dilute_strs = []
    dilute_conf = []
    structures = copy.deepcopy(endmembers)
    if supercell_matrix is not None:
        for i in structures:
            i.make_supercell(supercell_matrix)
    for strs in structures:
        dictstr = strs.as_dict()
        sub_lattice_list = [i['properties']['sublattice_sites'] for i in dictstr['sites']]
        num = [i for i in range(len(sub_lattice_list)) if i != 0 and sub_lattice_list[i] != sub_lattice_list[i - 1]]
        for i in num:
            dictstr = strs.as_dict()
            old_ele = dictstr['sites'][i]['species'][0]['element']
            sublattice_num = sub_lattice_list[i]
            element_dict = sublattice_dict[sublattice_num]
            for j in element_dict:
                if j != 'fix' and old_ele != j:
                    dictstr['sites'][i]['species'][0]['element'] = j
                    dilute_strs.append(Structure.from_dict(dictstr))

    for dilute in dilute_strs:
        comb = {}
        dictdilute = dilute.as_dict()
        for i in num:
            sfirst_ele = dictdilute['sites'][i]['species'][0]['element']
            ssecond_ele = dictdilute['sites'][i + 1]['species'][0]['element']
            site_sub = dictdilute['sites'][i]['properties']['sublattice_sites']
            if sfirst_ele == ssecond_ele:
                comb[site_sub] = sfirst_ele
            else:
                comb[site_sub] = [sfirst_ele, ssecond_ele]
        dilute_conf.append(comb)
    return dilute_strs, dilute_conf
