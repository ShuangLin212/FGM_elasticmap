"""Convienence functions for building endmembers"""
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from itertools import permutations, product, chain
from collections import Counter

def get_sublattice_information(structure, use_equivalent_atom=False):
    """
    Return defined sublattice_name based on wyckoff positions and other sublattice information

    Parameters
    ----------
    structure : pymatgen.Structure

    use_equivalent_atom: 
        From pymatgen function, it shows the wyckoff letter for each site. Some sites with the same
        wyckoff letters may still be symmetric inequivalent due to their coordinates. Here equivalent atom
        setting will distinguish symmetry inequivalent sites and provide more information for build up sublattice.
        The default setting is False. If you would like to get information based on equivalent atoms, use True.

    Returns
    -------
    true_sublattice_sites: list
        List of wyckoff letters for each site

    subl_model_name: list
        Name for each sublattice based on the name of wyckoff sites
    
    sublattice_site_ratio: list
        Site ratio for each sublattice
        
    Note:
    Need to consider order-disorder case
    """
    structure.replace_species({sp.name: "H" for sp in structure.species})
    sga = SpacegroupAnalyzer(structure)
    wyckoff_sites = sga.get_symmetry_dataset()['wyckoffs']
    equal_atom = sga.get_symmetry_dataset()['equivalent_atoms']
    num_wyckoff_sites = sorted(set(wyckoff_sites))
    num_eq_atom=sorted(set(equal_atom))
    if len(num_wyckoff_sites)!=len(num_eq_atom):
        print('Wyckoff sites may not be consistent with equivalent atoms in the structure, please check symmerty information and choose the sublattice model.')
    if use_equivalent_atom==True:
        sub_wyckoff_name=[]
        for i in num_eq_atom:
            sub_wyckoff_name.append(wyckoff_sites[i])
        replace_list=[]
        original_list=[]
        for i in sub_wyckoff_name:
            original_list.append(i)
            if i in replace_list:
                count_subl=dict(Counter(original_list))
                j=i+str(count_subl[i])
            else:
                j=i
            replace_list.append(j)
        replace_dict=dict(zip(num_eq_atom, replace_list))
        true_sublattice_sites=[replace_dict[i] for i in equal_atom]
    else:
        true_sublattice_sites = wyckoff_sites
    subl_model_name = sorted(set(true_sublattice_sites))
    sublattice_site_ratio=[]
    for i in subl_model_name:
        j=true_sublattice_sites.count(i)
        sublattice_site_ratio.append(j)
    return true_sublattice_sites, subl_model_name, sublattice_site_ratio

def get_templates(structure, wyckoff_site_list, subl_model_name, equivalent_sites=None):
    """
    Return templates of structure and configuration for substitution of endmembers

    Parameters
    ----------
    structure : pymatgen.Structure

    wyckoff_site_list: list
        List of wyckoff letter for each site, can be get from function get_sublattice_information. 
        If input manually, please make sure the order matches with the order of postions in pymatgen.Structrue

    subl_model_name: list
        Name for each sublattice based on the name of wyckoff sites.

    equivalent_sites: dict
        Set the equivalent sites when needed. For example, {'b': 'a'} means merge site 'a' and 'b' to the same sublattice 'a'.
    
    Returns
    -------
    template_structure: pymatgen.Structure

    template_configuration: list
        Template configuration from the template structure and sublattice model.
    """
    if equivalent_sites is not None:
        true_sublattice_sites=[equivalent_sites[i] if i in equivalent_sites else i for i in wyckoff_site_list]
        rep_sublattices=[equivalent_sites[i] if i in equivalent_sites else i for i in subl_model_name]
    else:
        true_sublattice_sites=wyckoff_site_list
        rep_sublattices=subl_model_name
    true_sublattices=sorted(set(rep_sublattices))
    site_list=[]
    dict_struct=structure.as_dict()
    for i in true_sublattice_sites:
        site_dict={'sublattice_sites':i}
        site_list.append(site_dict)
    for i in range(0,len(dict_struct['sites'])):
        dict_struct['sites'][i]['properties'].update(site_list[i])
        dict_struct['sites'][i]['species'][0]['element']=site_list[i]['sublattice_sites']
    template_configuration=[]
    for i in range(0,len(true_sublattices)):
        for element in dict_struct['sites']:
            if element['species'][0]['element']==true_sublattices[i]:
                element['species'][0]['element'] = str(Element.from_Z(i+1))
        template_configuration.append(str(Element.from_Z(i+1)))
    template_structure=Structure.from_dict(dict_struct)
    return template_structure, template_configuration


def get_endmembers_with_templates(template_structure, template_configuration, sublattice_configuration):
    """
    Return endmembers of the structure

    Parameters
    ----------
    template_structure: pymatgen.Structure

    template_configuration: list
        Configuration in the template structure, should be consitent with sublattice model, e.g., one unique
        element for each sublattice. If input manually, make sure be consistent with the sublattice_model_name 

    sublattice_configuration: list
        List of configurations for each sublattice

    Returns
    -------
    endmembers: list 
        List of structures in pymatgen.Structure
    """
    element_dict=dict(map(lambda x, y: [x, y], template_configuration, sublattice_configuration))
    subl_combination=[]
    for comb in product(*map(element_dict.get, sorted(element_dict))):
        subl_combination.append(list(comb))
    endmembers=[]
    for i in subl_combination:
        endmembers.append(substitute_configuration(template_structure, [template_configuration], [i]))
    return subl_combination, endmembers
