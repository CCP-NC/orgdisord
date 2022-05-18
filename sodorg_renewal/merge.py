'''This module implements several structure merging strategies 
to reduce a list of ASE atoms objects to the symmetry inequivalent ones.'''

# import modules
from tabnanny import verbose
import numpy as np
import warnings

# try to import tqdm if not define dummy version
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **kwargs):
            return x


def merge_structures(images, algo='symm', symops=None, use_disordered_only=True, symprec=1e-3, verbose=False, oxidation_states=None):
    '''Merge a list of ASE atoms objects to the symmetry inequivalent ones.
    
    Args:
        images (list): A list of ASE atoms objects.
        algo (str): The merging algorithm. Default is 'symm'.
        symops (list): A list of tuples: (rotation, translation).
                          Note: we skip the first one assuming it's the identity.
        use_disordered_only (bool): If True, use the disordered atoms (tag != 0)
                                    instead of the full structure when comparing.
        symprec (float): The tolerance for symmetry equivalence.
        verbose (bool): If True, print out the number of atoms in each structure.
        oxidation_states (dict): A dictionary of oxidation states for each element.
                                (Only needed for the 'ewald' algorithm.)
        
    Returns:
        list: A list of ASE atoms objects.
    '''
    
    # get the number of atoms in each structure
    natoms = [len(atoms) for atoms in images]
    
    # TODO does this help? 
    # sort the structures by the number of atoms
    # images = [images[i] 
    #           for i in sorted(range(len(natoms)),
    #                           key=natoms.__getitem__,
    #                           reverse=True)]

    # call the merging algorithm
    if algo == 'symm':
        # if symops is not provided, raise an error
        if symops is None:
            raise ValueError('symops must be provided for the symmetric merging algorithm')
        images = merge_symm(images, symops, symprec=symprec)
    elif algo == 'rematch': 
        images = merge_rematch(images, eps=symprec)
    elif algo == 'ewald':
        images = merge_ewald(images, oxidation_states, eps=symprec, parallel=True, verbose=verbose)

    else:
        raise ValueError('Unknown merging algorithm: {}'.format(algo))
    
    return images





def coords_match_symmops(coords_ref, coords, symops, symprec=1e-4):
        '''
        Check if two sets of coordinates match under any of the
        symmetry operations of the crystal
        and under permutations of the order of the atoms. 

        Args:
            coords_ref (array): A reference numpy array of fractional coordinates.
            coords (array): A numpy array of fractional coordinates.
            symops (list): A list of tuples: (rotation, translation).
                           Note: we skip the first one assuming it's the identity.
            symprec (float): The tolerance for the similarity.
        Returns:
            bool: True if the coordinates match, else False


        '''
        thesame = False
        # loop over all symmetry operations
        for symop in symops:
            rot, trans = symop
            a1 = np.mod(coords_ref, 1.0)
            a2 = np.mod((np.dot(rot, coords.T).T + trans), 1.0)
            if coords_match(a1, a2, symprec = symprec):
                thesame = True
                # break out of trying all the symmetry 
                # operations, we're done...
                break

        return thesame

def compare_ref_unmatched(ref, unmatched, symops, symprec=1e-4, parallel = True):
        '''
        compare lists of fractional coordinates with reference
        return the indices of the unmatched atoms
        '''
        if parallel:
            from multiprocessing import Pool, cpu_count
            # we can only run efficiently on about 100 structures per core
            # fewer than that and we'll just use serial
            ncores = min([cpu_count(), len(unmatched)//100])
        if parallel and ncores > 1:
            pool = Pool(processes=ncores)
            argument_list = [(ref, unmatched_i[1], symops, symprec) for unmatched_i in unmatched]
            jobs = [pool.apply_async(func=coords_match_symmops, args=(*argument,)) for argument in argument_list]
            pool.close()
            mask = []
            for job in jobs:
                mask.append(job.get())
        else:
            mask = [coords_match_symmops(ref, unmatched_i[1], symops, symprec) for unmatched_i in unmatched]
        matching_inds = np.where(mask)[0]

        return matching_inds

def merge_symm(supercell_images, symops, symprec = 1e-4, use_disordered_only =True):

        if use_disordered_only:
            frac_positions = []
            for atoms in supercell_images:
                disordered_idx = np.where(atoms.get_tags() != 0)[0]
                pos_f = atoms.get_scaled_positions()[disordered_idx]
                frac_positions.append(pos_f)
        else:
            frac_positions = [atoms.get_scaled_positions() for atoms in supercell_images]
        
        
        all_groups = []
        tagged_images = [[i,np.array(sites)] for i, sites in enumerate(frac_positions)]
        unmatched = tagged_images.copy()
        # set up progress bar
        pbar = tqdm(desc='Checking symmetry-equivalence', total=len(tagged_images))
        # loop over unmatched images
        while len(unmatched) > 0:
            i, ref_coords = unmatched.pop(0)
            matches = [i]
            inds = compare_ref_unmatched(ref_coords, unmatched, symops, symprec=symprec)
            matches.extend([unmatched[i][0] for i in inds])
            unmatched = [unmatched[i] for i in range(len(unmatched)) if i not in inds]
            all_groups.append(matches)
            # update the progress bar
            for match in matches:
                pbar.update(1)
        pbar.close()
        
        # print(all_groups)
        return [(supercell_images[g[0]], len(g)) for g in all_groups]



def rematcher(re, features1, features2, eps):
    '''Compare two sets of SOAP features and return True if the
       features match.
       
       Args:
            re (RE): The RE object.
            features1 (list): A list of SOAP features.
            features2 (list): A list of SOAP features.
            eps (float): The tolerance for the similarity.
        Returns:
            bool: True if the features match, else False

       '''
    re_kernel = re.create([features1, features2])
    val = re_kernel[0][1]

    return abs(val - 1) < eps 


def merge_rematch(supercell_images, eps=1e-2):
        '''
        This merging algorithm creates a fingerprint for each structure 
        bases on the SOAP representation and then used the 
        ReMatch kernel  to try to find equivalences. 

        It contains a number of parameters that would need to be fine tuned
        In my testing, it works pretty well for small systems, but 
        it overzealous at merging larger systems...

        More work is necessary! 
        '''
        from dscribe.descriptors import SOAP
        from dscribe.kernels import REMatchKernel

        nimages = len(supercell_images)
        # assumes all supercells have the same atoms set... is this always the case? 
        disordered_species_set = set([atom.symbol for atom in supercell_images[0] if atom.tag != 0 ] )


        desc = SOAP(species=disordered_species_set,
                    rcut=7.5,
                    nmax=2, 
                    lmax=1, 
                    sigma=0.5, 
                    periodic=True, 
                    crossover=True, 
                    sparse=False)

        print(f'generating descriptors for {nimages} atoms objects')
        features = [desc.create(sub, n_jobs=4) for sub in supercell_images]

        print('performing structure matching')
        # Any metric supported by scikit-learn will work: e.g. a Gaussian.
        re = REMatchKernel(metric="rbf", gamma=1, alpha=1, threshold=1e-6)
        
        unmatched = list(range(nimages))
        all_groups = []
        while len(unmatched) > 0:
            i = unmatched.pop(0)
            matches = [i]
            inds = [j for j in unmatched if rematcher(re, features[i], features[j], eps)]

            matches += inds
            unmatched = [u for u in unmatched if u not in matches]
            all_groups.append(matches)
        # return tuple of one structure per group and the group multiplicity
        return [(supercell_images[g[0]], len(g)) for g in all_groups]

def coords_match(a1, a2, symprec = 1e-4):
    """
    Returns True if a1 and a2 
    are equal apart from a reordering of the 
    rows. 
    
    """
    cmp = np.abs(a1[:, None] - a2) < symprec
    eq = np.all(cmp, axis=-1)
    number_of_matches = np.argwhere(eq).shape[0]
    if number_of_matches == a1.shape[0]:
        return True
    else:
        return False



def merge_ewald(supercell_images, oxidation_states, eps=1e-2, parallel = True, verbose=False):
    '''
    This merging algorithm calculates the ewald sum of each structure 
    and then uses the ewald sum to try to find equivalences.

    Args:
        supercell_images (list): A list of Atoms objects.
        oxidation_states (dictionary): A dictionary of oxidation states for each element.
             e.g. oxidation_states = {'N': -3.0, 'H': 1.0, 'C': 1.0, 'O': -2.0}
        eps (float): The tolerance eV for structure similarity
        parallel (bool): If True, use multiprocessing to speed up the calculation.

    Returns:
        list: A list of tuples containing the merged structure and the multiplicity.
    '''
    from pymatgen.analysis.ewald import EwaldSummation
    from pymatgen.io.ase import AseAtomsAdaptor
    
    

    # check if oxidation states are provided
    if oxidation_states is None:
        raise ValueError('Oxidation states must be provided')

    print('Adding oxidation states')
    # convert to pymatgen structure
    ewalds = []
    for i, atoms in enumerate(supercell_images):
        structure = AseAtomsAdaptor.get_structure(atoms)
        structure.add_oxidation_state_by_element(oxidation_states)
        # structures.append(structure)
        acc_factor=4.0
        ewalds.append(EwaldSummation(structure, acc_factor=acc_factor))
    print('Calculating Ewald energies')
    if parallel:
        from multiprocessing import Pool, cpu_count
        ncores = cpu_count()
        print(f'using {ncores} cores')

        pool = Pool(processes=ncores)

        energies = []
        for energy in tqdm(pool.imap(func=ewald_energy,
                           iterable=ewalds),
                           total=len(ewalds)):
            energies.append(energy)
    else:
        print('Csing serial')
        print(f'Calculating Ewald energies for {len(ewalds)} structures')
        energies = [ewald_energy(ewald) for ewald in tqdm(ewalds)]
    
    # compare the ewald sums
    unmatched = list(range(len(supercell_images)))
    all_groups = []
    while len(unmatched) > 0:
        i = unmatched.pop(0)
        matches = [i]
        inds = [j for j in unmatched if abs(energies[i] - energies[j]) < eps]
        matches += inds
        unmatched = [u for u in unmatched if u not in matches]
        all_groups.append(matches)
    # return tuple of one structure per group and the group multiplicity
    if verbose:
        print(all_groups)
    return [(supercell_images[g[0]], len(g)) for g in all_groups]

def ewald_energy(ewald):
        '''
        acc_factor (float): No. of significant figures each sum is
                converged to.
                (defaults to 12.0 in pymatgen. Here we just need a quick and dirty answer)
        '''
        acc_factor=3.0
        # print(structure)
        # ewald = EwaldSummation(structure, acc_factor=acc_factor)
        return ewald.total_energy