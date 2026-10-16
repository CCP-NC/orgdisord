"""This module implements several structure merging strategies
to reduce a list of ASE atoms objects to the symmetry inequivalent ones."""

import logging
from multiprocessing import Pool, cpu_count

import numpy as np
from ase.build.supercells import lattice_points_in_supercell
from scipy.spatial import cKDTree

# try to import tqdm if not define dummy version
try:
    from tqdm import tqdm
except ImportError:

    def tqdm(x, **kwargs):
        return x


logger = logging.getLogger("orgdisord.merge")


class NonDiscreteSitesError(ValueError):
    """The atoms do not sit on a discrete set of sites closed under the symmetry
    operations, so canonical labelling cannot be used."""


def _wrap(x):
    """Wrap fractional coordinates into [0, 1).

    A plain ``np.mod(x, 1.0)`` returns exactly 1.0 for tiny negative inputs
    (e.g. -1e-17), which a periodic ``cKDTree`` rejects.
    """
    x = np.mod(x, 1.0)
    return np.where(x >= 1.0, 0.0, x)


def supercell_symops(symops, supercell=(1, 1, 1)):
    """Convert unit-cell symmetry operations to the supercell frame.

    The supercell lattice is P @ (unit-cell lattice). A unit-cell operation
    ``u' = u R^T + t`` becomes ``x' = x (P R^T P^-1) + t P^-1`` in supercell
    fractional coordinates x = u P^-1. Only operations that map the supercell
    lattice onto itself (``P R^T P^-1`` integer) are symmetries of the supercell
    as a periodic crystal, so the others are dropped. The pure unit-cell
    translations that are distinct modulo the supercell are then composed with
    the remaining operations, which makes the set a closed group.

    Args:
        symops (list): List of (rotation_matrix, translation_vector) tuples in
            unit-cell fractional coordinates.
        supercell (None, length-3 sequence or 3x3 array): the diagonal repeat
            (na, nb, nc) or the integer transformation matrix P.

    Returns:
        list: List of (rotation_sc, translation_sc) tuples in supercell
        fractional coordinates.
    """
    if supercell is None:
        supercell = (1, 1, 1)
    arr = np.asarray(supercell)
    if arr.shape not in [(3,), (3, 3)]:
        raise ValueError(f"supercell must be length 3 or 3x3, got shape {arr.shape}")
    if not np.allclose(arr, np.round(arr)):
        raise ValueError(f"Supercell must contain integers only, got {supercell}")
    P = np.round(arr).astype(int)
    if P.ndim == 1:
        P = np.diag(P)
    det = int(round(np.linalg.det(P)))
    if det == 0:
        raise ValueError(f"Supercell matrix has zero determinant: {supercell}")
    if det < 0:
        P = -P
    if np.array_equal(P, np.eye(3, dtype=int)):
        return list(symops)

    Pinv = np.linalg.inv(P)
    shifts = lattice_points_in_supercell(P)

    out = []
    dropped = 0
    for rot, trans in symops:
        rot_sc = (P @ np.asarray(rot, dtype=float).T @ Pinv).T
        if not np.allclose(rot_sc, np.round(rot_sc), atol=1e-8):
            dropped += 1
            continue
        rot_sc = np.round(rot_sc)
        trans_sc = np.asarray(trans, dtype=float) @ Pinv
        for s in shifts:
            out.append((rot_sc, _wrap(trans_sc + s)))
    if dropped:
        logger.warning(
            f"{dropped} of {len(symops)} symmetry operations do not map the "
            "supercell lattice onto itself and were ignored: configurations "
            "related only by these are treated as distinct."
        )
    return out


def _extract_positions_and_symbols(images, use_disordered_only=True):
    """Return wrapped fractional positions and symbols for each image.

    If ``use_disordered_only`` is True but an image has no tagged (tag != 0)
    atoms, all of its atoms are used.
    """
    frac_positions = []
    symbols = []
    for atoms in images:
        current_symbols = np.array(atoms.get_chemical_symbols())
        current_positions = atoms.get_scaled_positions()
        if use_disordered_only:
            disordered_idx = np.where(atoms.get_tags() != 0)[0]
            if len(disordered_idx) == 0:
                disordered_idx = np.arange(len(atoms))
            frac_positions.append(_wrap(current_positions[disordered_idx]))
            symbols.append(current_symbols[disordered_idx])
        else:
            frac_positions.append(_wrap(current_positions))
            symbols.append(current_symbols)
    return frac_positions, symbols


def merge_structures(
    images,
    algo="symm",
    symops=None,
    supercell=None,
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
        algo (str): The merging algorithm. Default is 'symm': canonical labelling,
                    falling back (with a warning) to the much slower pairwise
                    comparison if the atoms are not on a discrete symmetric set of
                    sites. 'canonical' is the same without the fallback (raises
                    NonDiscreteSitesError); 'pairwise' always compares pairs.
                    Other options are 'rematch' and 'ewald'.
        symops (list): A list of tuples: (rotation, translation) in unit-cell coordinates.
        supercell (None, length-3 sequence or 3x3 integer matrix): The supercell the
                    images were built in: the diagonal repeat (na, nb, nc) or the
                    transformation matrix P (supercell lattice = P @ unit-cell
                    lattice). Used to convert symmetry operations to the supercell
                    frame; must match how the images were generated. Translations by
                    unit-cell lattice vectors are treated as equivalences.
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
    logger.debug(f"Algorithm: {algo}")

    if not check_species:
        logger.warning("WARNING: not checking species before comparing coordinates")

    if len(images) == 0:
        return ([], []) if return_group_indices else []

    if algo in ("symm", "canonical", "pairwise"):
        if symops is None:
            raise ValueError(
                f"symops must be provided for the '{algo}' merge algorithm"
            )
        sc_symops = supercell_symops(symops, supercell)
        groups = None
        if algo in ("symm", "canonical"):
            try:
                groups = merge_canonical(
                    images,
                    sc_symops,
                    symprec=symprec,
                    use_disordered_only=use_disordered_only,
                    check_species=check_species,
                )
            except NonDiscreteSitesError as e:
                if algo == "canonical":
                    raise
                logger.warning(
                    f"Canonical merge not applicable ({e}); falling back to the "
                    "much slower pairwise symmetry merge."
                )
        if groups is None:
            groups = merge_symm(
                images,
                sc_symops,
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
        raise ValueError(f"Unknown merging algorithm: {algo}")

    logger.debug(f"Merged into {len(groups)} groups:")
    for group in groups:
        logger.debug(f"    {group}")

    merged_images = [(images[g[0]], len(g)) for g in groups]
    if return_group_indices:
        return merged_images, groups
    else:
        return merged_images


def merge_canonical(
    supercell_images,
    symops,
    symprec=1e-4,
    use_disordered_only=True,
    check_species=True,
):
    """Group structures into symmetry-inequivalent orbits using canonical labelling.

    Builds a shared discrete site template across all images, closes it under the
    symmetry operations, and maps each image to an occupancy vector over the template.
    The canonical key for each image is the lexicographical minimum of the occupancy
    vector across all symmetry operations in the group.

    Args:
        supercell_images (list): List of ASE Atoms objects.
        symops (list): List of (rotation, translation) symmetry operations in
            supercell fractional coordinates (see ``supercell_symops``).
        symprec (float): Tolerance for site mapping, in fractional coordinates
            of the supercell.
        use_disordered_only (bool): If True, compare only disordered atoms (tag != 0).
        check_species (bool): If True, differentiate atom species.

    Returns:
        list of list of int: Grouped indices of equivalent images.

    Raises:
        NonDiscreteSitesError: if the atoms are not on a discrete set of sites
            closed under ``symops``.
    """
    if not supercell_images:
        return []

    pos, sym = _extract_positions_and_symbols(
        supercell_images, use_disordered_only=use_disordered_only
    )

    # 1. Build unique site template from all images
    sites_list = []
    tree = None
    for p in pos:
        for pt in p:
            if tree is None:
                sites_list.append(pt)
                tree = cKDTree(np.array(sites_list), boxsize=1.0)
            else:
                d, _ = tree.query(pt, p=np.inf, distance_upper_bound=symprec)
                if not np.isfinite(d):
                    sites_list.append(pt)
                    tree = cKDTree(np.array(sites_list), boxsize=1.0)

    sites = np.array(sites_list)

    # 2. Close template under all symmetry operations
    max_sites = 25 * max(len(p) for p in pos)
    iterations = 0
    while True:
        iterations += 1
        if iterations > 30 or len(sites) > max_sites:
            raise NonDiscreteSitesError(
                f"Template closure did not converge (sites={len(sites)}, iter={iterations}). "
                "Symmetry operations may not form a closed group on this lattice."
            )
        added = 0
        for r, t in symops:
            cand = _wrap(sites @ r.T + t)
            for pt in cand:
                d, _ = tree.query(pt, p=np.inf, distance_upper_bound=symprec)
                if not np.isfinite(d):
                    sites = np.vstack([sites, pt])
                    tree = cKDTree(sites, boxsize=1.0)
                    added += 1
        if added == 0:
            break

    tree = cKDTree(sites, boxsize=1.0)
    if tree.query_pairs(2 * symprec, p=np.inf):
        raise NonDiscreteSitesError(
            "Template sites are closer than 2 * symprec; site assignment is ambiguous."
        )
    nsites = len(sites)

    # 3. Permutation of template sites under each symmetry operation
    perms = []
    for r, t in symops:
        trans_sites = _wrap(sites @ r.T + t)
        d, idx = tree.query(trans_sites, p=np.inf, distance_upper_bound=symprec)
        if not np.all(np.isfinite(d)) or len(np.unique(idx)) != len(idx):
            raise NonDiscreteSitesError(
                "Template sites do not permute bijectively under symops."
            )
        perms.append(idx)
    perms = np.array(perms)  # shape: (n_symops, nsites)

    # 4. Map chemical symbols to integer codes
    all_species = sorted({s for s_list in sym for s in s_list})
    codes = {s: k + 1 for k, s in enumerate(all_species)}

    # 5. Compute canonical key for each image
    keys = []
    for p, s in zip(pos, sym):
        d, idx = tree.query(p, p=np.inf, distance_upper_bound=symprec)
        if not np.all(np.isfinite(d)) or len(np.unique(idx)) != len(idx):
            raise NonDiscreteSitesError(
                "Image atoms do not map uniquely onto template sites."
            )

        labels = np.zeros(nsites, dtype=np.int16)
        if check_species:
            labels[idx] = [codes[x] for x in s]
        else:
            labels[idx] = 1

        permuted = labels[perms]
        min_key = min(row.tobytes() for row in permuted)
        keys.append(min_key)

    # 6. Group indices by canonical key
    groups_dict = {}
    for i, key in enumerate(keys):
        groups_dict.setdefault(key, []).append(i)

    return list(groups_dict.values())


def merge_symm(
    supercell_images,
    symops,
    symprec=1e-4,
    use_disordered_only=True,
    quiet=False,
    check_species=True,
):
    """Pairwise symmetry merge using periodic KD-tree matching.

    Args:
        supercell_images (list): List of ASE Atoms objects.
        symops (list): List of (rotation, translation) symmetry operations.
        symprec (float): Tolerance in fractional coordinates.
        use_disordered_only (bool): If True, use only disordered atoms.
        quiet (bool): If True, suppress progress bar.
        check_species (bool): If True, check species.

    Returns:
        list of list of int: Grouped indices of equivalent images.
    """
    pos, sym = _extract_positions_and_symbols(
        supercell_images, use_disordered_only=use_disordered_only
    )

    unmatched = list(range(len(pos)))
    all_groups = []

    pbar = tqdm(
        desc="Checking symmetry-equivalence", total=len(unmatched), disable=quiet
    )

    while len(unmatched) > 0:
        ref_idx = unmatched.pop(0)
        matches = [ref_idx]
        tree_ref = cKDTree(pos[ref_idx], boxsize=1.0)
        ref_sym = sym[ref_idx]
        n_ref = len(ref_sym)

        remaining_unmatched = []
        for other_idx in unmatched:
            oth_p = pos[other_idx]
            oth_s = sym[other_idx]
            if len(oth_s) != n_ref:
                remaining_unmatched.append(other_idx)
                continue

            if check_species and sorted(ref_sym) != sorted(oth_s):
                remaining_unmatched.append(other_idx)
                continue

            is_match = False
            for rot, trans in symops:
                a2 = _wrap(oth_p @ rot.T + trans)
                d, mapping = tree_ref.query(a2, p=np.inf, distance_upper_bound=symprec)
                if not np.all(np.isfinite(d)):
                    continue
                if len(np.unique(mapping)) != len(mapping):
                    continue
                if check_species and not np.array_equal(oth_s, ref_sym[mapping]):
                    continue
                is_match = True
                break

            if is_match:
                matches.append(other_idx)
            else:
                remaining_unmatched.append(other_idx)

        unmatched = remaining_unmatched
        all_groups.append(matches)
        for _ in matches:
            pbar.update(1)

    pbar.close()
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
    except ImportError as err:
        raise ImportError(
            "The merge_rematch algorithm requires the dscribe package. "
            "Please install orgdisord with the ML option."
            "e.g. \npip install orgdisord[ML]"
        ) from err

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
    for atoms in supercell_images:
        structure = AseAtomsAdaptor.get_structure(atoms)
        structure.add_oxidation_state_by_element(oxidation_states)
        # structures.append(structure)
        acc_factor = 4.0
        ewalds.append(EwaldSummation(structure, acc_factor=acc_factor))
    logger.info("Calculating Ewald energies")
    if parallel:
        ncores = cpu_count()
        logger.info(f"using {ncores} cores")

        with Pool(processes=ncores) as pool:
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

    # print(structure)
    # ewald = EwaldSummation(structure, acc_factor=acc_factor)
    return ewald.total_energy
