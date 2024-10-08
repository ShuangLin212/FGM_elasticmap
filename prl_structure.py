
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from collections import Counter

class PRLStructure(Structure):
    """A pymatgen Structure object, with some customizations for ESPEI.
    """

    def __init__(self, *args, **kwargs):
        """Create a Structure object, with some customizations for ESPEI

        Parameters
        ----------
        args :
            args to pass to Structure
        sublattice_configuration : [[str]]
            Sublattice configuration  e.g. `[['Fe', 'Ni'], ['Fe']]`.
        sublattice_occupancies : [[float]]
            Fraction of the sublattice each element in the configuration has e.g. `[[0.3333, 0.6666], [1]]`.
        sublattice_site_ratios : [float]
            Ratios of sublattice multiplicity  e.g. `[3, 1]`.
        kwargs :
            kwargs to pass to Structure
        """
        self.sublattice_configuration = kwargs.pop('sublattice_configuration', None)
        self.sublattice_occupancies = kwargs.pop('sublattice_occupancies', None)
        self.sublattice_site_ratios = kwargs.pop('sublattice_site_ratios', None)
        self.wyckoff_sites = kwargs.pop('wyckoff_sites', None)
        super(PRLStructure, self).__init__(*args, **kwargs)

    def __eq__(self, other):
        """
        self and other are equivalent if the sublattice models are equal

        Parameters
        ----------
        other : PRLStructure
        """
        if not isinstance(other, PRLStructure):
            return False
        subl_config = self.sublattice_configuration == other.sublattice_configuration
        subl_site_ratios = self.sublattice_site_ratios == other.sublattice_site_ratios
        subl_occupancies = self.sublattice_occupancies == other.sublattice_occupancies
        wyckoff_sites = self.wyckoff_sites == other.wyckoff_sites
        return subl_config and subl_site_ratios and subl_occupancies and wyckoff_sites

    @property
    def espei_sublattice_configuration(self):
        """Return ESPEI-formatted sublattice model [['a', 'b'], 'a'] for the concrete case."""
        canonicalize_sublattice = lambda sl: sl[0] if len(sl) == 1 else sl
        return [canonicalize_sublattice(sl) for sl in self.sublattice_configuration]

    @property
    def espei_sublattice_occupancies(self):
        """Return ESPEI-formatted sublattice occupancies [[0.3333, 0.6666], 1] for the concrete case."""
        canonicalize_sublattice = lambda sl: sl[0] if len(sl) == 1 else sl
        return [canonicalize_sublattice(sl) for sl in self.sublattice_occupancies]

    def as_dict(self, verbosity=1, fmt=None, **kwargs):
        d = super(PRLStructure, self).as_dict(verbosity=verbosity, fmt=fmt, **kwargs)
        d['sublattice_configuration'] = self.sublattice_configuration
        d['sublattice_occupancies'] = self.sublattice_occupancies
        d['sublattice_site_ratios'] = self.sublattice_site_ratios
        return d

    @classmethod
    def from_dict(cls, d, fmt=None):
        struct = super(PRLStructure, cls).from_dict(d, fmt=fmt)
        struct.sublattice_configuration = d.get('sublattice_configuration')
        struct.sublattice_occupancies = d.get('sublattice_occupancies')
        struct.sublattice_site_ratios = d.get('sublattice_site_ratios')
        return struct

    @classmethod
    def from_structure(cls, structure, equivalent_sites=None):
        struct = PRLStructure.from_dict(structure.as_dict())
        structure = Structure.from_dict(structure.as_dict())
        sga = SpacegroupAnalyzer(structure)
        wyckoff_sites = sga.get_symmetry_dataset()['wyckoffs']
        equal_atom = sga.get_symmetry_dataset()['equivalent_atoms']
        num_wyckoff_sites = sorted(set(wyckoff_sites))
        num_eq_atom = sorted(set(equal_atom))

        if num_wyckoff_sites == num_eq_atom:
            true_sites = wyckoff_sites
        else:
            subl_wyckoff_name = []
            for i in num_eq_atom:
                subl_wyckoff_name.append(wyckoff_sites[i])
            replace_list = []
            original_list = []
            for i in subl_wyckoff_name:
                original_list.append(i)
                if i in replace_list:
                    count_subl = dict(Counter(original_list))
                    j = i + str(count_subl[i])
                else:
                    j = i
                replace_list.append(j)
            replace_dict = dict(zip(num_eq_atom, replace_list))
            true_sites = [replace_dict[i] for i in equal_atom]

        true_sublattices = sorted(set(true_sites))
        if equivalent_sites is not None:
            combined_sublattices = ['-'.join(sorted(sites)) for sites in equivalent_sites]

            def match_subl(candidate):
                for subl in combined_sublattices:
                    if candidate in subl:
                        return subl
                return candidate

            new_subl_model = sorted(set([match_subl(subl) for subl in true_sublattices]))
        else:
            new_subl_model = true_sublattices

        config = []
        occ = []
        ratios = []
        for subl in new_subl_model:
            species_frequency_dict = {}
            for site, wyckoff_site in zip(struct.sites, true_sites):
                if '-' in subl:
                    if wyckoff_site in subl:
                        species = site.specie.name.upper()
                        species_frequency_dict[species] = species_frequency_dict.get(species, 0) + 1
                else:
                    if wyckoff_site == subl:
                        species = site.specie.name.upper()
                        species_frequency_dict[species] = species_frequency_dict.get(species, 0) + 1
            total_subl_occupation = sum(species_frequency_dict.values())
            subl_species = sorted(set(species_frequency_dict.keys()))
            subl_occpancy = [species_frequency_dict[sp] / total_subl_occupation for sp in subl_species]
            config.append(subl_species)
            occ.append(subl_occpancy)
            ratios.append(total_subl_occupation)

        struct.sublattice_configuration = config
        struct.sublattice_occupancies = occ
        struct.sublattice_site_ratios = ratios
        struct.wyckoff_sites = sorted(set(true_sites))
        return struct

    @staticmethod
    def reindex_sublattice(new_indices, subl_model, subl_occupancies, subl_site_ratios):
        if sorted(new_indices) != list(range(len(subl_model))):
            raise ValueError('Passed re-indexing indices do not match the sublattice model indices.')
        new_subl_model = [subl_model[i] for i in new_indices]
        new_subl_occupancies = [subl_occupancies[i] for i in new_indices]
        new_subl_site_ratios = [subl_site_ratios[i] for i in new_indices]
        return new_subl_model, new_subl_occupancies, new_subl_site_ratios

    def reindex(self, new_indices):
        self.sublattice_configuration, self.sublattice_occupancies, self.sublattice_site_ratios = \
            PRLStructure.reindex_sublattice(new_indices, self.sublattice_configuration, self.sublattice_occupancies, self.sublattice_site_ratios)
