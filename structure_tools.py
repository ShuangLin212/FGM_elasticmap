
"""
Tools for substituting structures and generating metadata
"""

from copy import deepcopy

def sort_x_by_y(x, y):
    """Sort a list of x in the order of sorting y"""
    return [xx for _, xx in sorted(zip(y, x), key=lambda pair: pair[0])]

def canonicalize_config(configuration, occupancies):
    """
    Return canonicalized (sorted) configurations and occupancies.

    Parameters
    ----------
    configuration : list of lists
        DFTTK-style configuration
    occupancies :
        DFFTK-style occupancies

    Returns
    -------
    tuple
        Tuple of canonical (configuration, occupancies)
    """
    new_occupancies = [sort_x_by_y(occ, config) for occ, config in zip(occupancies, configuration)]
    new_configuration = [sorted(config) for config in configuration]
    return (new_configuration, new_occupancies)

def get_density_from_pt(ele_list):
    """
    Get density(g/cm^3) from periodictable package

    Parameters
    ----------
        ele_list : list/dict
            The list of elements, e.g. ['Nb', 'Ti']/{'Nb': 3, 'Ti': 1}
    Returns
    -------
        density_dict : dict
            Dictionary of {element: density}, e.g. {'Nb': 8.57, 'Ti': 4.507}. 
    """
    from pymatgen.core.periodic_table import Element
    density_dict = {}
    for ele in ele_list:
        density_dict[ele] = float(Element(ele).density_of_solid)/1000.
    return density_dict

def get_ele_list_from_struct(struct):
    """
    Get elements list from pymatgen structure objective

    Parameters
    ----------
        struct : pymatgen objective
            The structure
    Returns
    -------
        ele_list : [str]
            The list of elements
    """
    ele_list = []
    for ele in struct.species:
        ele_list.append(str(ele))
    return ele_list

def scale_struct(struct):
    """Scale the structure according to the weighted average density of each element."""
    species_amnt_dict = struct.composition.get_el_amt_dict()
    density_dict = get_density_from_pt(species_amnt_dict)
    expected_density = float(sum([density_dict[species]*amnt for species, amnt in species_amnt_dict.items()]))/struct.composition.num_atoms
    current_density = struct.density
    current_volume = struct.volume
    expected_volume = current_volume/expected_density*current_density
    struct.scale_lattice(float(expected_volume))
    return struct

def gen_replacement_dict(old_config, new_config):
    """Create a pymatgen replacement dict based on old and new sublattice configurations."""
    replacement_dict = {}
    for new_subl, old_subl in zip(new_config, old_config):
        for new_atom, old_atom in zip(new_subl, old_subl):
            replacement_dict[old_atom] = new_atom
    return replacement_dict

def substitute_configuration(template_structure, template_config, config, check_sorting=True):
    """
    Replace the species in the template structure by switching the template_config elements for the config elements.
    """
    struct = deepcopy(template_structure)
    struct.replace_species(gen_replacement_dict(template_config, config))
    scale_struct(struct)
    return struct

def substitute_configuration_with_metadata(template_structure, template_config, config, occupation, phase_name, site_ratios):
    """
    Replace the species in the template structure by switching the template_config elements for the config elements.
    Returns a tuple of the structure and metadata.
    """
    struct = substitute_configuration(template_structure, template_config, config, check_sorting=False)
    config, occupation = canonicalize_config(config, occupation)
    metadata = {'phase_name': phase_name, 'sublattice': {'configuration': config, 'occupancies': occupation, 'site_ratios': site_ratios}}
    return struct, metadata
