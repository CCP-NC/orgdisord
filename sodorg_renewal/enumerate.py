'''Module to enumerate ordered structures given a cif object containing disorder information.'''
import numpy as np
import warnings
from ase.spacegroup.spacegroup import SpacegroupValueError
from ase import Atoms, Atom
from ase.geometry import get_duplicate_atoms
from soprano.properties.linkage import Molecules
import spglib
import itertools
from tqdm import tqdm
from scipy.spatial.distance import pdist
import logging
from sodorg_renewal.utils import reload_as_molecular_crystal

logger = logging.getLogger("sodorg.enumerate")



def binary_to_idx(config, maxgroups):
    return sum([j*(maxgroups**i) for i,j in list(enumerate(reversed(config)))])

def select_configs(config, supercell, maxgroups):
    '''
    from the supercell config, return indices of each
    component primitive cell 
    e.g. 
    a 2x2x1 supercell of a Z=4 system
    would have configs of length 16.
    This function chunks them into 4
    sets of 4 and returns the integer representation
    (0-15) of each 'binary' number.
    So for example:
    [1,1,1,1, 0,0,0,0, 1,1,0,0, 0,0,0,1]
    would return
    [15, 0, 12, 1] 

    '''
    na, nb, nc = supercell
    nsuper = na * nb * nc
    chunk = len(config) // nsuper
    # should equal Z
    return [binary_to_idx(c, maxgroups) for c in  chunks(config, chunk)]



def chunks(lst, n, offset=0):
    """Yield successive n-sized chunks from lst."""
    if offset != 0:
    # rotate copy of list by offset
        temp = lst[offset:] + lst[:offset]
    else:
        temp = lst
    for i in np.arange(0, len(lst), n):
        yield temp[i:i + n ]

class OrderedfromDisordered:
    def __init__(self, disordered_structure, symprec=1e-4,quiet = False):
        '''
        Args:
            disordered_structure (DisorderedStructure)
            symprec (float): tolerance for symmetry finding
            quiet (bool): if True, suppresses warnings

        '''
        self.quiet = quiet
        self.disordered_structure = disordered_structure
        self.symprec = symprec
        self.logger = logging.getLogger("sodorg.enumerate")
        self.logger.debug("\n\n")
        self.logger.debug("------------------------------")
        self.logger.debug("--- ENUMERATING STRUCTURES ---")
        self.logger.debug("------------------------------")


    def _get_config_symbols_coords(self, config, asslabel):
        '''
        takes config as an int (often binary) list
        of length Z (number of disordered molecular units in cell)
        and the index of the chosen assembly

        Args:
            config (list): list of integers representing the config
            asslabel (str): Assembly label e.g. 'A'

        Returns:
            [
                symbols (list): list of symbols for the config
                sites (list): list of sites for the config
                tags (list): list of tags for the config
                labels (list): list of labels for the config
            ]

        '''
        assembly = self.disordered_structure.get_assembly(asslabel)
        disorder_groups = assembly.disorder_groups
        disorder_group_labels = sorted(group.label for group in disorder_groups)
        
        # TODO: fix tags by assembly
        # get the assembly index for the specified assembly label
        iassembly = 1000
        
        
        # a few aliases
        # how many disorder groups in this assembly?
        ndisordergroups = len(disorder_groups)
        nassemblies = self.disordered_structure.get_number_of_assemblies()

        sites = []
        symbols = []
        tags = []
        labels = []
        # loop over groups
        for igroup in range(ndisordergroups):
            group = disorder_groups[igroup]
            assert len(group.atoms) > 0
            group_symmops = group.symmetry_operations
            # self.logger.debug(f"{len(group_symmops)} group_symmops: {group_symmops}")
            # if Zprime < 1:
            #     ngroups  = int(1/Zprime)
            #     # chunk the symbols and coords into ngroups
            #     group_symmops = [chunk[c] for chunk, c in zip(chunks(ops, ngroups), config)]
            #     group_symmops = self.cif.get_group_symmops(asslabel, disorder_group_labels[igroup])
            # else:
            # at this point len(group_symmops) = len(config) 
            # if len(group_symmops) != len(config):
            #     raise ValueError(f"nops ({len(group_symmops)}) != len(config) ({len(config)})")

            if len(group_symmops) == len(config):
                # choose the subset of symmops for this group given config
                group_config_symmops = [group_symmops[i] for i, c in enumerate(config) if c == igroup]
            else:
                # Zprime < 1 so config instead needs to index into subset of group symops
                # e.g. config [0,0,1,0] for would select those subsets of symops for this group
                # probably need to select the inverse for the next group...
                
                group_config_symmops = [group_symmops[c][i] for i, c in enumerate(config)]
                
            # tag = igroup + 1 if nassemblies == 1 else 1e2*(iassembly+1) + igroup + 1
            tag = group.tag
            [symbols_g, sites_g, tags_g, labels_g] = self._get_group_config_sites(tag, group, group_config_symmops)
            symbols += symbols_g
            sites += sites_g
            tags += tags_g
            labels += labels_g

        return [symbols, sites, tags, labels]

    def _get_group_config_sites(self, tag, group, group_config_symmops):
        sites = []
        symbols = []
        tags = []
        labels = []
        # loop over selected symmetry operations -- apply them to the group
        for iops, op in enumerate(group_config_symmops):
            # g = config[iops]
            # # Which group should we pick? 
            # # if Zprime <1 or disorder_group_labels[igroup] < 0:
            # #     grouplabel = disorder_group_labels[igroup]
            # # else:
            # grouplabel = disorder_group_labels[igroup]
            # group = assembly.get_disorder_group(grouplabel)
            unique_positions = group.atoms.get_scaled_positions()
            unique_symbols   = group.atoms.get_chemical_symbols()
            unique_labels    = group.atoms.get_array("labels")
            kinds = range(len(unique_symbols))
                # loop over site in each disorder group:
            for kind, pos in enumerate(unique_positions):
                    # apply symmetry operation
                rot, trans = op
                site = np.mod(np.dot(rot, pos) + trans, 1.)                  
                
                if not sites:
                    sites.append(site)
                    symbols.append(unique_symbols[kind])
                    tags.append(tag)
                    labels.append(unique_labels[kind])
                    continue
                t = site - sites
                mask = np.all(
                        (abs(t) < self.symprec) | (abs(abs(t) - 1.0) < self.symprec), axis=1)
                if np.any(mask):
                    inds = np.argwhere(mask).flatten()
                    if len(inds) > 1:
                            raise SpacegroupValueError(
                                    'Found multiple equivalent sites!'.format())
                    else:
                        pass
                else:
                    sites.append(site)
                    symbols.append(unique_symbols[kind])
                    tags.append(tag)
                    labels.append(unique_labels[kind])
        return [symbols, sites, tags, labels]
    
    def get_config(self, config, asslabel):
        '''
        Return a fully ordered atoms object for the specified config

        '''
        # set up base atoms object to which sites will be 
        # added (just an empty box).
        atoms = Atoms(cell=self.disordered_structure.cell, pbc=True)
        
        # add in disordered sites, ordered according to config, (tag = 1e2*(iassembly+1) + igroup + 1)
        symbols, sites, tags, labels = self._get_config_symbols_coords(config, asslabel=asslabel)
        for site, symbol, tag in zip(sites, symbols, tags):
            atom = Atom(symbol=symbol, position = atoms.cell.T.dot(site), tag=tag)
            atoms.append(atom)
        atoms.set_array('labels', np.array(labels))

        return atoms

    def get_all_configs(self, exclude_ordered = False):
        '''Generate all possible configs for the primitive cell.'''
        
        # TODO this could be made more efficient!


        # a few aliases
        cell = self.disordered_structure.cell
        Z = self.disordered_structure.Z
        # Zprime = self.cif.Zprime
        nops_all = self.disordered_structure.spacegroup.nsymop
        nassemblies=  self.disordered_structure.get_number_of_assemblies()
        correlated_assemblies = self.disordered_structure.correlated_assemblies

        
        max_spacegroup_kinds = 0
        # set up base atoms object to which sites will be 
        # added. 
        if exclude_ordered:
        # just an empty box
            atoms = Atoms(cell=cell, pbc=True)
        else: 
            atoms = self.disordered_structure.ordered_atoms.copy()
            max_spacegroup_kinds = atoms.get_array('spacegroup_kinds').max()
            # give ordered sites a tag of 0
            atoms.set_tags(0)


        # what's the final number of configs we will generate?
        #TODO generalise!
        ngroups_per_assembly = self.disordered_structure.get_number_of_disorder_groups_per_assembly()
        nconfigs_per_assembly = [N**Z for N in ngroups_per_assembly]
        nconfigs = nconfigs_per_assembly[0] # defaults to first ngroups

        if correlated_assemblies:
            if nassemblies > 1:
                self.logger.debug('Treating the assemblies as correlated. \n'
                'This means that the same index from each assembly is chosen for each configuration.')
                if len(set(nconfigs_per_assembly)) != 1:
                    raise ValueError('Error: the number of configurations per assembly is not the same for all assemblies.\n'
                                        'Please set correlated_assemblies to False.')
        else:
            if nassemblies > 1:
                self.logger.debug('Treating the assemblies as uncorrelated '
                '(i.e. the assemblies are independent) and we will therefore generate lots of structures!')
                nconfigs = np.product(nconfigs_per_assembly)
        self.logger.debug(f'Generating {nconfigs} configs in primitive cell.')
        if nconfigs > 1e4:
            self.logger.warn(f'Warning: {nconfigs} is a large number of configs. '
                              'This may take a while!')

        
        all_configs = []
        for assembly in self.disordered_structure.disorder_assemblies:
            asslabel = assembly.label
            # the disorder groups for this assembly
            disorder_groups = assembly.disorder_groups
            # normally ngroups = ndisorder groups
            ngroups = len(disorder_groups)
            group_nops = [len(group.symmetry_operations) for group in disorder_groups]
            self.logger.debug(f'Group n ops {group_nops}')
            
            if len(set(group_nops)) != 1:
                raise ValueError('Error: the number of symmetry operations per group '
                                f' is not the same for all groups in assembly {asslabel}:\n'
                                f'{group_nops}'
                                    )


            # but when Zprime is less than 1, we need to be more careful!
            # if Zprime < 1: # means more symmetry operations than we have Z
            #     # in this case the actual number of groups (1/Zprime) will be larger than the number of disorder groups
            #     assert 1/Zprime >= ngroups

            #     ngroups = int(1 / Zprime)
            #     self.logger.warn(f'Special case for nops != Z. (i.e. Zprime < 1)'
            #         f'nops = {nops_all} '
            #         f'Z = {Z}, '
            #         f'Zprime = {Zprime}\n'
            #         'This is less well-tested so proceed with caution.'
            #         'If you encounter any issues, please contact the developers.')
            
            
            
            self.logger.debug(f'Assembly {asslabel} has {ngroups} groups')

            if ngroups == 1:
                # Special disorder site case with only one group
                # -> each symmetry operation generated a new config. 
                config_max = group_nops[0]
            else:
                config_max = ngroups

            # generate all possible configs for this assembly
            config_indices = itertools.product(list(range(config_max)), repeat=Z)
            assembly_configs = [self.get_config(config_idx, asslabel) for config_idx in config_indices]
            # update max_spacegroup_kinds tags
            for config in assembly_configs:
                labels = config.get_array('labels')
                u, spacegroup_kinds = np.unique(labels, return_inverse=True)
                config.set_array('spacegroup_kinds', spacegroup_kinds + max_spacegroup_kinds + 1)
            self.logger.debug(f'Assembly {asslabel} has {len(assembly_configs)} configs')
            if all_configs:
                if correlated_assemblies:
                    # then just add the config atoms to the existing ones!
                    assert len(assembly_configs) == len(all_configs)
                    for iconfig in range(len(assembly_configs)):
                        temp = all_configs[iconfig].copy()
                        temp.extend(assembly_configs[iconfig])
                        all_configs[iconfig] = temp
                else:
                    ref_configs = [config.copy() for config in all_configs]
                    all_configs = []
                    for iconfig in range(len(ref_configs)):
                        for jconfig in range(len(assembly_configs)):
                            temp = ref_configs[iconfig].copy()
                            temp.extend(assembly_configs[jconfig])
                            all_configs.append(temp)
            else:
                # must the be first assembly -- all_configs is empty so let's start from scratch
                for config in assembly_configs:
                    temp = atoms.copy()
                    temp.extend(config)
                    all_configs.append(temp)
            
        self.logger.debug(f"Generated a total of {len(all_configs)} primitive configs")
        return all_configs


    def get_supercell_configs(self,
                              supercell, 
                              maxiters = 5000, 
                              exclude_ordered = False, 
                              random_configs=False,
                              return_configs=False):
        '''
        loop over supercell cells,
        add in one of the ndisordergroups^Z configurations per cell
        if return_configs is True, return the configs as well as the supercells
        if molecular_crystal is True, then the structure is reloaded as a molecular crystal
             -- trying to keep the pieces in tact
        '''
        # pre-compute all primitive cells
        self.logger.debug('Pre-computing all primitive configurations...')
        # TODO: include a maxiters argument here?
        images = self.get_all_configs(exclude_ordered)
        if self.disordered_structure.molecular_crystal:
            self.logger.debug(f'Reloading {len(images)} images as molecular crystals')
            if len(images) > 256:
                # this seems to be get pretty slow for large numbers of images
                self.logger.warn('Warning: reloading molecular crystals is slow for this many images.\n'
                'Consider skipping this step with the --not_molecular_crystal flag!\n'
                'We will use a cheaper method to reload the images as molecular crystals instead\n'
                ' -- assuming the connectivity is the same for each image')
                cheap_method = True
            else:
                cheap_method = False
            images = reload_as_molecular_crystal(images, cheap=cheap_method)

        # some aliases
        Z = self.disordered_structure.Z
        cell = self.disordered_structure.cell
        na, nb, nc = supercell
        # we previously just took ndisordergroups = self.cif.ndisordergroups, but 
        # this gets tricky when we have multiple assemblies etc. Safer to just take this:
        ndisordergroups = int(np.round(np.exp(np.log(len(images))/Z)))
        # total number of configs given this supercell:
        ncombinations = ndisordergroups**(Z*na*nb*nc)
        self.logger.debug(f'Found {len(images)} nimages,  {ncombinations} configs for this supercell (Z={Z}, na={na}, nb={nb}, nc={nc}), ndisordergroups={ndisordergroups}')
        # how many configs to actually generate:
        if random_configs:
            # just take maxiters
            n_configs = maxiters
        else:
            # take whichever is smallest betwee
            # maxiters and ncombinations
            n_configs = min([maxiters, ncombinations])
        
        
        # this iterator generates all 
        # possible lists of 0s and 1s etc of length Z*na*nb*nc
        all_combinations = itertools.product(list(range(ndisordergroups)), repeat=Z*na*nb*nc)

        all_supercells = []
        all_configs = []
        if random_configs:
            self.logger.info(f'Generating {n_configs} random configurations in a {supercell} supercell:')
        else:
            self.logger.info(f'Generating {n_configs} out of the {ncombinations} possible configurations in the {supercell} supercell:')

        for i in tqdm(range(n_configs), disable=self.quiet):
            if random_configs:
                config = np.random.randint(ndisordergroups, size=Z*na*nb*nc)
            else:
                config = np.array(all_combinations.__next__())
            # make ncells copies of the relevant config
            self.logger.debug(f'         {config}')
            supercell_atoms = Atoms(cell = cell * supercell, pbc = True)
            for icell, c in enumerate(select_configs(config, supercell=supercell, maxgroups=ndisordergroups)):
                # make a copy of prim config c
                temp = images[c].copy()
                # work out what the translation vector should be:
                ia = icell % na
                ib = (icell // na) % nb
                ic = icell // (na*nb)
                R = cell.T.dot([ia,ib,ic])
                # translate primitive copy to correct place
                temp.translate(R)
                # add this block into the supercell
                supercell_atoms.extend(temp)
            # TODO
            # check for overlapping atoms
            # cutoff = 0.25
            # dists = pdist(supercell_atoms.get_positions(), 'sqeuclidean')
            # if (dists < cutoff**2).any():
            #     continue
            # else:
            all_supercells.append(supercell_atoms)
            all_configs.append(config)

            
            # if get_duplicate_atoms(supercell_atoms, cutoff=0.5, delete=False) is None:
            # else:
            #     continue
        if return_configs:
            return all_supercells, all_configs
        return all_supercells
    









