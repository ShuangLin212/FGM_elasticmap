
"""
The sqs module handles converting abstract SQS Structure objects to concrete structures.
The SQS are regular pymatgen Structures with the species named according to sublattice and species type.
These species in pymatgen Structures are named to `Xab`, which corresponds to atom `B` in sublattice `a`.
"""

from __future__ import division
import copy
import itertools
import pymatgen as pmg
from pymatgen.core import Structure

class AbstractSQS(Structure):
    """A pymatgen Structure with special features for SQS."""

    def __init__(self, *args, **kwargs):
        self.sublattice_model = kwargs.pop('sublattice_model', None)
        self._sublattice_names = kwargs.pop('sublattice_names', None)
        super(AbstractSQS, self).__init__(*args, **kwargs)

    @property
    def normalized_sublattice_site_ratios(self):
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = self.composition.as_dict()
        site_ratios = [[comp_dict['X'+name+e+'0+']/self.num_sites for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

    @property
    def sublattice_site_ratios(self):
        subl_model = self.sublattice_model
        subl_names = self._sublattice_names
        comp_dict = {k: int(v) for k, v in self.composition.reduced_composition.as_dict().items()}
        site_ratios = [[comp_dict['X'+name+e+'0+'] for e in subl] for subl, name in zip(subl_model, subl_names)]
        return site_ratios

    def get_concrete_sqs(self, subl_model, scale_volume=True):
        if len(subl_model) != len(self.sublattice_model):
            raise ValueError('Concrete sublattice model {} does not match size of abstract sublattice model {}'.format(subl_model, self.sublattice_model))

        replacement_dict = {}
        site_occupancies = []
        for abstract_subl, concrete_subl, subl_name, subl_ratios in zip(self.sublattice_model, subl_model, self._sublattice_names, self.sublattice_site_ratios):
            sublattice_ratio_sum = sum(subl_ratios)
            sublattice_occupancy_dict = {}
            for abstract_specie, concrete_specie, site_ratio in zip(abstract_subl, concrete_subl, subl_ratios):
                specie = 'X' + subl_name + abstract_specie
                replacement_dict[specie] = concrete_specie
                sublattice_occupancy_dict[concrete_specie] = sublattice_occupancy_dict.get(concrete_specie, 0) + site_ratio/sublattice_ratio_sum
            site_occupancies.append(sublattice_occupancy_dict)

        self_copy = copy.deepcopy(self)
        self_copy.replace_species(replacement_dict)

        if scale_volume:
            fractional_comp = dict(self_copy.composition.fractional_composition)
            estimated_density = sum((pmg.core.periodic_table.Element(component).density for component in self_copy.composition.elements)) / 1000
            self_copy.scale_lattice(float((self_copy.volume/estimated_density)*self_copy.density))

        sublattice_configuration = [sorted(set(subl)) for subl in subl_model]
        sublattice_occupancies = [[occupancies[specie] for specie in subl] for occupancies, subl in zip(site_occupancies, sublattice_configuration)]
        site_ratios = [sum(ratios) for ratios in self.sublattice_site_ratios]

        concrete_sqs = PRLStructure.from_sites(self_copy.sites)
        concrete_sqs.sublattice_configuration = sublattice_configuration
        concrete_sqs.sublattice_occupancies = sublattice_occupancies
        concrete_sqs.sublattice_site_ratios = site_ratios
        return concrete_sqs

    def get_endmember_space_group_info(self, symprec=1e-2, angle_tolerance=5.0):
        endmember_subl = [['X' + subl_name for _ in subl] for subl, subl_name in zip(self.sublattice_model, self._sublattice_names)]
        endmember_speices = {specie for subl in endmember_subl for specie in subl}
        real_species_dict = {abstract_specie: real_specie for abstract_specie, real_specie in zip(endmember_speices, pmg.core.periodic_table._pt_data.keys())}
        endmember_subl = [[real_species_dict[specie] for specie in subl] for subl in endmember_subl]
        endmember_struct = self.get_concrete_sqs(endmember_subl, scale_volume=False)
        return endmember_struct.get_space_group_info(symprec=symprec, angle_tolerance=angle_tolerance)

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(AbstractSQS, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_model'] = self.sublattice_model
        d['sublattice_names'] = self._sublattice_names
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        d['symmetry'] = {'symbol': self.get_endmember_space_group_info()[0], 'number': self.get_endmember_space_group_info()[1]}
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        sqs = super(AbstractSQS, cls).from_dict(d, fmt=fmt)
        sqs.sublattice_model = d.get('sublattice_model')
        sqs._sublattice_names = d.get('sublattice_names')
        return sqs

def enumerate_sqs(structure, subl_model, scale_volume=True, skip_on_failure=False):
    if len(subl_model) != len(structure.sublattice_model):
        raise ValueError('Passed sublattice model ({}) does not agree with the passed structure ({})'.format(subl_model, structure.sublattice_model))
    possible_subls = [itertools.product(subl, repeat=len(abstract_subl)) for subl, abstract_subl in zip(subl_model, structure.sublattice_model)]
    unique_subl_models = itertools.product(*possible_subls)

    unique_sqs = []
    unique_configurations_occupancies = []
    for model in unique_subl_models:
        proposed_sqs = structure.get_concrete_sqs(model, scale_volume)
        proposed_config_occupancy = (proposed_sqs.sublattice_configuration, proposed_sqs.sublattice_occupancies)
        if proposed_config_occupancy not in unique_configurations_occupancies:
            unique_configurations_occupancies.append(proposed_config_occupancy)
            unique_sqs.append(proposed_sqs)
    return unique_sqs
