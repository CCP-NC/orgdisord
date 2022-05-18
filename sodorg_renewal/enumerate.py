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


def binary_to_idx(config):
    return sum([j*(2**i) for i,j in list(enumerate(reversed(config)))])

def select_configs(config, supercell):
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
    return [binary_to_idx(c) for c in  chunks(config, chunk)]

def reload_as_molecular_crystal(images):
    '''
    takes in a list of atoms,
    identifies molecules,
    unwraps molecules so that they're 'connected' 
    across periodic boundaries
    '''
    images_new = []
    for atoms in images:
        mols = Molecules.get(atoms)
        temp = mols[0].subset(atoms, use_cell_indices=True)
        for mol in mols[1:]:
            temp.extend(mol.subset(atoms, use_cell_indices=True))
        images_new.append(temp)
    return images_new


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in np.arange(0, len(lst), n):
        yield lst[i:i + n]

class OrderedfromDisordered:
    def __init__(self, cif, symprec=1e-4, verbose=False):
        '''
        cif must be an instance of CifParser from parse_cif_file.py

        '''
        self.cif = cif
        self.symprec = symprec
        self.verbose = verbose

    def _get_config_symbols_coords(self, config, onduplicates='error'):
        '''
        takes config as a binary list (i.e. list of ones and zeros)
        of length Z (number of disordered molecular units in cell)

        onduplicates must be one of "keep", "replace", "warn" or "error"
            - says what to do if duplicate atoms are created. defaults to error

        returns lists of [symbols, coordinates] for each generated site for the specified
        configuration
        '''
        # a few aliases
        Z = self.cif.Z
        nops_all = self.cif.nops
        ops  = self.cif.symops
        ndisordergroups = self.cif.ndisordergroups
        nassemblies     = self.cif.nassemblies
        asymmetric_scaled_coords = self.cif.asymmetric_scaled_coords
        asymmetric_symbols = self.cif.asymmetric_symbols
        groups = self.cif.disorder_groups
        
        if nops_all != Z:
            Zprime = Z / nops_all
            chunksize = int(1 / Zprime)
            assert chunksize == 2 # TODO: does this generalise beyond this? 
            assert nassemblies == 1 # TODO: is this always the case? 
            ops = [chunk[c] for chunk, c in zip(chunks(ops, chunksize), config)]

        # now loop over symmetry operations   
        sites = []
        symbols = []
        tags = []
        for iops, op in enumerate(ops):
            
            # loop over groups
            for igroup in range(ndisordergroups):
                
                if nops_all != Z:
                    iassembly = 0
                else:
                    iassembly = config[iops]
                
                group = np.array(groups[iassembly][igroup])
                if len(group) > 0:
                    unique_positions = asymmetric_scaled_coords[group]
                    unique_symbols   =       asymmetric_symbols[group]
                    kinds = range(len(unique_symbols))
                    # loop over site in each disorder group:
                    for kind, pos in enumerate(unique_positions):
                        # apply symmetry operation
                        rot, trans = op
                        site = np.mod(np.dot(rot, pos) + trans, 1.)                  

                        if not sites:
                            sites.append(site)
                            symbols.append(unique_symbols[kind])
                            tags.append(igroup + 1)
                            continue
                        t = site - sites
                        mask = np.all(
                            (abs(t) < self.symprec) | (abs(abs(t) - 1.0) < self.symprec), axis=1)
                        if np.any(mask):

                            # TODO: simplify this ? We don't need all these options i think
                            inds = np.argwhere(mask).flatten()
                            for ind in inds:
                                # then we would just add the same thing again -> skip
                                if kinds[ind] == kind:
                                    pass
                                elif onduplicates == 'keep':
                                    pass
                                elif onduplicates == 'replace':
                                    kinds[ind] = kind
                                elif onduplicates == 'warn':
                                    warnings.warn('scaled_positions %d and %d '
                                                    'are equivalent' %
                                                    (kinds[ind], kind))
                                elif onduplicates == 'error':
                                    raise SpacegroupValueError(
                                        'scaled_positions %d and %d are equivalent' %
                                        (kinds[ind], kind))
                                else:
                                    raise SpacegroupValueError(
                                        'Argument "onduplicates" must be one of: '
                                        '"keep", "replace", "warn" or "error".')
                        else:
                            sites.append(site)
                            symbols.append(unique_symbols[kind])
                            tags.append(igroup + 1)
        # print(config, groups, Z, len(symbols), nops_all)
        return [symbols, sites, tags]
    
    def get_config(self, config, exclude_ordered=False):
        '''
        Return a fully ordered atoms object for the specified config

        exclude_ordered (boolean) sets whether or not to include the non-disordered part of the crystal

        '''
        # set up base atoms object to which sites will be 
        # added. 
        if exclude_ordered:
            # just an empty box
            atoms = Atoms(cell=self.cif.cell, pbc=True)
        else: 
            atoms = self.cif.ordered_atoms.copy()
            # give ordered sites a tag of 0
            atoms.set_tags(0)


        # add in disordered sites, ordered according to config, tag = 1
        symbols, sites, tags = self._get_config_symbols_coords(config)
        for site, symbol, tag in zip(sites, symbols, tags):
            atom = Atom(symbol=symbol, position = self.cif.cell.T.dot(site), tag=tag)
            atoms.append(atom)

        return atoms
    def get_all_configs(self, exclude_ordered = False):
        if self.cif.ndisordergroups != 2:
            raise ValueError('Error: we cannot yet handle cases of ngroups != 2')
        # generate all possible ndisordergroups^Z combinations
        Z = self.cif.Z
        all_combinations = np.array(list(itertools.product(list(range(self.cif.ndisordergroups)), repeat=Z)))
        return [self.get_config(config, exclude_ordered) for config in all_combinations]


    def get_supercell_configs(self,
                              supercell, 
                              maxiters = 5000, 
                              exclude_ordered = False, 
                              random_configs=False):
        '''
        loop over supercell cells,
        add in one of the ndisordergroups^Z configurations per cell
        '''
        # some aliases
        Z = self.cif.Z
        cell = self.cif.cell
        na, nb, nc = supercell
        # total number of configs given this supercell:
        ncombinations = 2**(Z*na*nb*nc)
        # how many configs to actually generate:
        if random_configs:
            # just take maxiters
            n_configs = maxiters
        else:
            # take whichever is smallest betwee
            # maxiters and ncombinations
            n_configs = min([maxiters, ncombinations])
        
        # pre-compute all primitive cells
        images = self.get_all_configs(exclude_ordered)
        images_mol = reload_as_molecular_crystal(images)
        
        # this iterator generates all 
        # possible lists of 0s and 1s of length Z*na*nb*nc
        all_combinations = itertools.product(list(range(2)), repeat=Z*na*nb*nc)
        
        all_supercells = []
        print(f'Generating {n_configs} in {supercell} supercell')
        for i in tqdm(range(n_configs), disable=not self.verbose):
            if random_configs:
                config = np.random.randint(2, size=Z*na*nb*nc)
            else:
                config = all_combinations.__next__()
            # config = gen_random_config(Z*na*nb*nc)
            # make ncells copies of the relevant config
            supercell_atoms = Atoms(cell = cell * supercell, pbc = True)
            for icell, c in enumerate(select_configs(config, supercell=supercell)):
                # make a copy of prim config c
                temp = images_mol[c].copy()
                # work out what the translation vector should be:
                ia = icell % na
                ib = (icell // na) % nb
                ic = icell // (na*nb)
                R = cell.T.dot([ia,ib,ic])
                # translate primitive copy to correct place
                temp.translate(R)
                # add this block into the supercell
                supercell_atoms.extend(temp)
            # check for overlapping atoms
            # cutoff = 0.25
            # dists = pdist(supercell_atoms.get_positions(), 'sqeuclidean')
            # if (dists < cutoff**2).any():
            #     continue
            # else:
            all_supercells.append(supercell_atoms)

            
            # if get_duplicate_atoms(supercell_atoms, cutoff=0.5, delete=False) is None:
            # else:
            #     continue
        
        return all_supercells
    









