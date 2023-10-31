"""This module implements several structure merging strategies 
to reduce a list of ASE atoms objects to the symmetry inequivalent ones."""


# import modules
import numpy as np
from typing import List, Tuple

import logging
from multiprocessing import Pool, cpu_count

# try to import tqdm if not define dummy version
try:
    from tqdm import tqdm
except ImportError:

    def tqdm(x, **kwargs):
        return x


logger = logging.getLogger("orgdisord.merge")


def merge_structures(
    images,
    algo="symm",
    symops=None,
    use_disordered_only=True,
    symprec=1e-3,
    oxidation_states=None,
    return_group_indices=False,
    quiet=False,
    check_species=True,
):
    """Merge a list of ASE atoms objects to the symmetry inequivalent ones.

    Args:
        images (list): A list of ASE atoms objects.
        algo (str): The merging algorithm. Default is 'symm'.
        symops (list): A list of tuples: (rotation, translation).
                          Note: we skip the first one assuming it's the identity.
        use_disordered_only (bool): If True, use the disordered atoms (tag != 0)
                                    instead of the full structure when comparing.
        symprec (float): The tolerance for symmetry equivalence.
        oxidation_states (dict): A dictionary of oxidation states for each element.
                                (Only needed for the 'ewald' algorithm.)
        return_group_indices (bool): If True, return the indices of the merged structures
                                    in addition to the merged structures.
        quiet (bool): If True, suppress progress bar.
        check_species (bool): If True, check that the species match before comparing
                              the coordinates.

    Returns:
        Tuple of one structure per group and the group multiplicity. If return_group_indices is True,
        the indices of the merged structures in the original list are also returned.

    """

    logger.debug("\n\n")
    logger.debug("-----------------------------")
    logger.debug("---  MERGING STRUCTURES  ---")
    logger.debug("-----------------------------")
    logger.debug("Algorithm: {}".format(algo))

    if not check_species:
        logger.warning("WARNING: not checking species before comparing coordinates")

    # get the number of atoms in each structure
    natoms = [len(atoms) for atoms in images]

    # TODO does this help?
    # sort the structures by the number of atoms
    # images = [images[i]
    #           for i in sorted(range(len(natoms)),
    #                           key=natoms.__getitem__,
    #                           reverse=True)]

    # call the merging algorithm
    if algo == "symm":
        # if symops is not provided, raise an error
        if symops is None:
            raise ValueError(
                "symops must be provided for the symmetric merging algorithm"
            )
        groups = merge_symm(
            images,
            symops,
            use_disordered_only=use_disordered_only,
            symprec=symprec,
            quiet=quiet,
            check_species=check_species,
        )
    elif algo == "rematch":
        groups = merge_rematch(images, eps=symprec, quiet=quiet)
    elif algo == "ewald":
        groups = merge_ewald(
            images, oxidation_states, eps=symprec, parallel=True, quiet=quiet
        )

    else:
        raise ValueError("Unknown merging algorithm: {}".format(algo))

    logger.debug(f"Merged into {len(groups)} groups:")
    for group in groups:
        logger.debug(f"    {group}")

    merged_images = [(images[g[0]], len(g)) for g in groups]
    if return_group_indices:
        return merged_images, groups
    else:
        return merged_images


def coords_match_symmops(
    ref: Tuple[int, np.ndarray, List],
    to_match: Tuple[int, np.ndarray, List],
    symops,
    symprec=1e-4,
    check_species=True,
):
    """
    Check if two sets of coordinates match under any of the
    symmetry operations of the crystal
    and under permutations of the order of the atoms.

    Args:
        ref ((int, np.array, list)): A reference tuple of (index, coords, symbols) .
        to_match ((int, np.array, list)): A tuple of (index, coords, symbols) to check against the ref.
        symops (list): A list of tuples: (rotation, translation).
                       Note: we skip the first one assuming it's the identity.
        symprec (float): The tolerance for the similarity.

    Returns:
        bool: True if the coordinates match, else False


    """
    thesame = False
    _, coords_ref, symbols_ref = ref
    _, coords_to_match, symbols_to_match = to_match

    # check the species match
    if check_species:
        # Set of species must be the same (ignoring number and order)
        if set(symbols_ref) != set(symbols_to_match):
            return False
        else:
            # check that the number of each species matches
            # Note that this still allows for different orderings
            # of the species. Later we will check that the orderings
            # match, but we can't do that until we've checked the
            # symmetry operations to see which positions map to which
            # positions
            for symbol in set(symbols_ref):
                if symbols_ref.count(symbol) != symbols_to_match.count(symbol):
                    return False

    # loop over all symmetry operations
    a1 = np.mod(coords_ref, 1.0)
    for symop in symops:
        rot, trans = symop
        a2 = np.mod((np.dot(rot, coords_to_match.T).T + trans), 1.0)
        matching_coords, mapping = coords_match(
            a1, a2, symprec=symprec, return_indices=True
        )
        if matching_coords:
            if check_species:
                # make sure the species match when reordered
                # according to the positions
                if not all(
                    [
                        symbols_ref[mapping[i]] == symbols_to_match[i]
                        for i in range(len(symbols_ref))
                    ]
                ):
                    logger.debug("species don't match despite coords matching")
                    logger.debug("mapping: {}".format(mapping))
                    continue
            thesame = True

            # break out of trying all the symmetry
            # operations, we're done...
            break

    return thesame


def compare_ref_unmatched(
    ref, unmatched, symops, symprec=1e-4, parallel=True, check_species=True
):
    """
    compare lists of fractional coordinates with reference
    return the indices of the unmatched atoms

    """

    if parallel:
        # we can only run efficiently on about 100 structures per core
        # fewer than that and we'll just use serial
        ncores = min([cpu_count(), len(unmatched) // 100])
    if parallel and ncores > 1:
        pool = Pool(processes=ncores)
        argument_list = [
            (ref, unmatched_i, symops, symprec, check_species)
            for unmatched_i in unmatched
        ]
        jobs = [
            pool.apply_async(func=coords_match_symmops, args=(*argument,))
            for argument in argument_list
        ]
        pool.close()
        mask = []
        for job in jobs:
            mask.append(job.get())
    else:
        mask = [
            coords_match_symmops(ref, unmatched_i, symops, symprec, check_species)
            for unmatched_i in unmatched
        ]
    matching_inds = np.where(mask)[0]

    return matching_inds


def merge_symm(
    supercell_images,
    symops,
    symprec=1e-4,
    use_disordered_only=True,
    quiet=False,
    check_species=True,
):
    frac_positions = []
    symbols = []
    for atoms in supercell_images:
        current_symbols = atoms.get_chemical_symbols()
        current_positions = atoms.get_scaled_positions()
        if use_disordered_only:
            disordered_idx = np.where(atoms.get_tags() != 0)[0]
            frac_positions.append(current_positions[disordered_idx])
            symbols.append([current_symbols[i] for i in disordered_idx])
        else:
            frac_positions.append(current_positions)
            symbols.append(current_symbols)

    all_groups = []
    # create a list of (index, frac_positions, symbols) tuples
    unmatched = [(i, frac_positions[i], symbols[i]) for i in range(len(symbols))]

    # set up progress bar
    pbar = tqdm(
        desc="Checking symmetry-equivalence", total=len(unmatched), disable=quiet
    )
    # loop over unmatched images
    while len(unmatched) > 0:
        ref = unmatched.pop(0)
        matches = [ref[0]]
        inds = compare_ref_unmatched(
            ref, unmatched, symops, symprec=symprec, check_species=check_species
        )
        matches.extend([unmatched[i][0] for i in inds])
        unmatched = [unmatched[i] for i in range(len(unmatched)) if i not in inds]
        all_groups.append(matches)
        # update the progress bar
        for match in matches:
            pbar.update(1)
    pbar.close()

    # return groups
    return all_groups


def rematcher(re, features1, features2, eps):
    """Compare two sets of SOAP features and return True if the
    features match.

    Args:
         re (RE): The RE object.
         features1 (list): A list of SOAP features.
         features2 (list): A list of SOAP features.
         eps (float): The tolerance for the similarity.

     Returns:
         bool: True if the features match, else False

    """

    re_kernel = re.create([features1, features2])
    val = re_kernel[0][1]

    return abs(val - 1) < eps


def merge_rematch(supercell_images, eps=1e-3, quiet=False, parallel=True):
    """
    This merging algorithm creates a fingerprint for each structure
    bases on the SOAP representation and then used the
    ReMatch kernel  to try to find equivalences.

    It contains a number of parameters that would need to be fine tuned
    In my testing, it works pretty well for small systems, but
    it overzealous at merging larger systems...

    More work is necessary!

    """
    # try to import dscribe
    # if it's not available, raise an error, and suggest installing orgdisord with the
    # ML option
    try:
        from dscribe.descriptors import SOAP
        from dscribe.kernels import REMatchKernel
    except ImportError:
        raise ImportError(
            "The merge_rematch algorithm requires the dscribe package. "
            "Please install orgdisord with the ML option."
            "e.g. \npip install orgdisord[ML]"
        )

    ncores = 1
    if parallel:
        # we can only run efficiently on about 100 structures per core
        # fewer than that and we'll just use serial
        ncores = min([cpu_count(), len(supercell_images) // 100 + 1])
    logger.info(f"Running ReMatch algo using {ncores} cores")
    nimages = len(supercell_images)
    # assumes supercells [0] and [-1] will cover all unique species
    disordered_species_set = set(
        [
            [atom.symbol for atom in atoms if atom.tag != 0]
            for atoms in [supercell_images[0], supercell_images[-1]]
        ]
    )

    desc = SOAP(
        species=disordered_species_set,
        rcut=7.5,
        nmax=2,
        lmax=1,
        sigma=0.5,
        periodic=True,
        crossover=True,
        sparse=False,
    )
    logger.info(f"generating descriptors for {nimages} atoms objects")
    features = desc.create(supercell_images, n_jobs=ncores, verbose=not quiet)

    logger.info("performing structure matching")
    # Any metric supported by scikit-learn will work: e.g. sigmoid.
    re = REMatchKernel(
        metric="rbf", gamma=1, alpha=1, threshold=1e-3, kernel_params={"n_jobs": ncores}
    )

    unmatched = list(range(nimages))
    all_groups = []
    while len(unmatched) > 0:
        i = unmatched.pop(0)
        matches = [i]
        inds = [j for j in unmatched if rematcher(re, features[i], features[j], eps)]

        matches += inds
        unmatched = [u for u in unmatched if u not in matches]
        all_groups.append(matches)
    # return groups
    return all_groups


def coords_match(a1, a2, symprec=1e-4, return_indices=False):
    """
    Returns True if a1 and a2
    are equal apart from a reordering of the
    rows.

    Args:
        a1 (np.array): An array of coordinates.
        a2 (np.array): An array of coordinates.
        symprec (float): The tolerance for the similarity.
        return_indices (bool): If True, return the mapping indices of the matching atoms.
                               i.e. an array of indices such that a1[indices] == a2

    """

    # Normalise a1 and a2 to the unit cell
    a1 = np.mod(a1, 1.0)
    a2 = np.mod(a2, 1.0)

    cmp = np.abs(a1[:, None] - a2) < symprec
    eq = np.all(cmp, axis=-1)
    number_of_matches = np.argwhere(eq).shape[0]
    is_match = number_of_matches == a1.shape[0]
    if return_indices:
        indices = np.argwhere(eq)[:, 1]
        return is_match, indices
    return is_match


def merge_ewald(
    supercell_images, oxidation_states, eps=1e-2, parallel=True, quiet=False
):
    """
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
    """
    from pymatgen.analysis.ewald import EwaldSummation
    from pymatgen.io.ase import AseAtomsAdaptor

    # check if oxidation states are provided
    if oxidation_states is None:
        raise ValueError("Oxidation states must be provided")

    logger.info("Adding oxidation states")
    # convert to pymatgen structure
    ewalds = []
    for i, atoms in enumerate(supercell_images):
        structure = AseAtomsAdaptor.get_structure(atoms)
        structure.add_oxidation_state_by_element(oxidation_states)
        # structures.append(structure)
        acc_factor = 4.0
        ewalds.append(EwaldSummation(structure, acc_factor=acc_factor))
    logger.info("Calculating Ewald energies")
    if parallel:
        ncores = cpu_count()
        logger.info(f"using {ncores} cores")

        pool = Pool(processes=ncores)

        energies = []
        for energy in tqdm(
            pool.imap(func=ewald_energy, iterable=ewalds),
            total=len(ewalds),
            disable=quiet,
        ):
            energies.append(energy)
    else:
        logger.info("Serial calculation")
        logger.info(f"Calculating Ewald energies for {len(ewalds)} structures")
        energies = [ewald_energy(ewald) for ewald in tqdm(ewalds, disable=quiet)]

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

    # return groups
    return all_groups


def ewald_energy(ewald):
    """
    acc_factor (float): No. of significant figures each sum is
            converged to.
            (defaults to 12.0 in pymatgen. Here we just need a quick and dirty answer)
    """

    acc_factor = 3.0
    # print(structure)
    # ewald = EwaldSummation(structure, acc_factor=acc_factor)
    return ewald.total_energy
