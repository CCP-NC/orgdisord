"""Module to enumerate ordered structures given a cif object containing disorder information."""

import itertools
import logging

import numpy as np
from ase import Atom, Atoms
from ase.build.supercells import find_optimal_cell_shape, lattice_points_in_supercell
from ase.data import covalent_radii
from ase.neighborlist import NeighborList
from ase.spacegroup.spacegroup import SpacegroupValueError
from scipy.optimize import linear_sum_assignment
from tqdm import tqdm

from orgdisord.utils import random_product, reload_as_molecular_crystal

logger = logging.getLogger("orgdisord.enumerate")


def _build_neighborlist_for_atoms(atoms, cutoff_scale=1.25):
    """
    Build an ASE NeighborList using covalent radii * cutoff_scale.
    This returns the neighborlist object; note that we DO NOT automatically
    interpret every neighbor as a bond — caller can filter pairs (e.g. skip H-H).
    """
    numbers = atoms.get_atomic_numbers().astype(int)
    # safe fallback handled by covalent_radii from ASE
    radii = np.array(
        [float(covalent_radii[z]) if z < len(covalent_radii) else 0.7 for z in numbers]
    )
    per_atom_cutoffs = radii * cutoff_scale
    nl = NeighborList(
        per_atom_cutoffs.tolist(), self_interaction=False, bothways=True, skin=0.0
    )
    nl.update(atoms)
    return nl


def connected_components_indices(atoms, cutoff_scale=1.25, exclude_hh=True):
    """
    Return list-of-lists of atom indices per connected component, using covalent cutoffs.
    If exclude_hh is True, H-H neighbor pairs are ignored when forming connectivity.
    """
    nl = _build_neighborlist_for_atoms(atoms, cutoff_scale=cutoff_scale)
    N = len(atoms)
    visited = np.zeros(N, dtype=bool)
    comps = []
    numbers = atoms.get_atomic_numbers().astype(int)

    for i in range(N):
        if visited[i]:
            continue
        queue = [i]
        comp = []
        visited[i] = True
        while queue:
            idx = queue.pop(0)
            comp.append(idx)
            neighs, offsets = nl.get_neighbors(idx)
            for nb in neighs:
                # Skip H-H edges if requested
                if exclude_hh and numbers[idx] == 1 and numbers[nb] == 1:
                    continue
                if not visited[nb]:
                    visited[nb] = True
                    queue.append(nb)
        comps.append(comp)

    return comps


def component_centroids(atoms, comps, mass_weighted=False):
    """
    Return centroids (Nx3) for components `comps`.
    If mass_weighted True, use atomic masses.
    """
    pos = atoms.get_positions()
    centroids = []
    if mass_weighted:
        masses = atoms.get_masses()
    for comp in comps:
        coords = pos[comp]
        if mass_weighted:
            w = masses[comp][:, None]
            c = (coords * w).sum(axis=0) / w.sum()
        else:
            c = coords.mean(axis=0)
        centroids.append(c)
    return np.array(centroids)


def _nearest_image_vectors(centroids, ref_centroids, cell, search=1):
    """
    For every (centroid, reference) pair return the integer lattice shift n that
    minimises |ref - (centroid + n @ cell)|, and that minimal distance.

    The search is centred on the rounded fractional difference, so it is exact
    however many cells apart the two points are (`search` only covers the
    residual error for strongly skewed cells).

    Returns:
        shifts (C, R, 3) int array, dist (C, R) float array
    """
    cell = np.asarray(cell)
    inv_cell = np.linalg.inv(cell)
    diff = ref_centroids[None, :, :] - centroids[:, None, :]  # (C, R, 3)
    base = np.round(diff @ inv_cell).astype(int)
    rng = range(-search, search + 1)
    window = np.array(list(itertools.product(rng, rng, rng)), dtype=int)  # (W, 3)
    cand = base[:, :, None, :] + window[None, None, :, :]  # (C, R, W, 3)
    resid = diff[:, :, None, :] - cand @ cell
    dist = np.linalg.norm(resid, axis=-1)  # (C, R, W)
    best = dist.argmin(axis=-1)  # (C, R)
    ic, ir = np.indices(best.shape)
    return cand[ic, ir, best], dist[ic, ir, best]


def align_components_to_reference(
    atoms, reference_centroids, cell, cutoff_scale=1.25, search=1, mass_weighted=False
):
    """
    Align connected components (molecules) from `atoms` to reference_centroids by
    integer-cell translations. Returns a new ASE Atoms copy.

    Each component is matched to a distinct reference centroid by minimising the
    total nearest-image distance (optimal assignment), then shifted by the whole
    lattice vector that puts its centroid in the same lattice image as its
    reference. Components only ever move by whole cell vectors, so the periodic
    structure is unchanged.

    The result deliberately is NOT re-wrapped afterwards: the same physical
    molecular site must sit in the same lattice image in every config, which is
    what makes configs safe to mix tile by tile in a supercell.

    Components left unmatched (more components than references) are not moved.
    """
    A = atoms.copy()
    comps = connected_components_indices(A, cutoff_scale=cutoff_scale)
    if len(comps) == 0:
        return A
    centroids = component_centroids(A, comps, mass_weighted=mass_weighted)
    ref = np.asarray(reference_centroids)
    shifts, dist = _nearest_image_vectors(centroids, ref, cell, search=search)

    rows, cols = linear_sum_assignment(dist)
    new_positions = A.get_positions().copy()
    for ic, ir in zip(rows, cols):
        new_positions[comps[ic]] += shifts[ic, ir] @ np.asarray(cell)
    A.set_positions(new_positions)
    return A


def make_whole_molecule(atoms, cutoff_scale=1.25, preserve_pbc=False, exclude_hh=True):
    """
    Unwrap connected molecules so each connected component is contiguous in Cartesian space.

    - atoms: ASE Atoms (can be periodic). Function works from current positions and cell.
    - cutoff_scale: multiplier for (r_i + r_j) to determine bond cutoff between atoms i and j.
    - preserve_pbc: if True, preserve atoms.pbc; otherwise set pbc=False on returned Atoms.
    - exclude_hh: if True, ignore H-H pairs when building connectivity graph (prevents false H-H bonds).

    Returns a copy of `atoms` with positions moved by integer multiples of cell vectors
    so connected components are placed compactly (no pieces split across periodic images).
    """
    if len(atoms) == 0:
        return atoms.copy()

    A = atoms.copy()
    cell = np.asarray(A.get_cell())  # 3x3 array
    orig = A.get_positions().copy()
    numbers = A.get_atomic_numbers().astype(int)
    N = len(A)

    # Build neighborlist (per-atom cutoffs)
    nl = _build_neighborlist_for_atoms(A, cutoff_scale=cutoff_scale)

    # shifts[i] is the accumulated integer lattice shift applied to atom i.
    # The neighbour list returns, for each bond i-j, the exact image offset o
    # such that orig[j] + o @ cell is bonded to orig[i], so a bonded neighbour
    # simply inherits shifts[i] + o. This is exact however far a molecule
    # extends from its seed atom (even over several cells).
    shifts = np.zeros((N, 3), dtype=int)
    visited = np.zeros(N, dtype=bool)
    for start in range(N):
        if visited[start]:
            continue

        # BFS through this connected component
        queue = [start]
        visited[start] = True
        while queue:
            i = queue.pop(0)
            neighs, offsets = nl.get_neighbors(i)
            for j, off in zip(neighs, offsets):
                # Skip H-H pairs if configured
                if exclude_hh and numbers[i] == 1 and numbers[j] == 1:
                    continue
                if visited[j]:
                    continue
                shifts[j] = shifts[i] + off
                visited[j] = True
                queue.append(j)

    # write back positions and clear PBC (unless preserve_pbc is True)
    A.set_positions(orig + shifts @ cell)
    if not preserve_pbc:
        A.set_pbc([False, False, False])

    return A


def wrap_components_into_cell(atoms):
    """
    Shift each connected component (molecule) by whole cell vectors so that its
    centroid lies inside the cell. Molecules are kept whole; atoms of molecules
    larger than (or straddling) the cell may still lie outside [0, 1).

    Only use this on a single reference structure: wrapping several configs
    independently is discontinuous for molecules whose centroid is near a cell face.

    Returns a copy of `atoms`.
    """
    A = atoms.copy()
    cell = np.asarray(A.get_cell())
    inv_cell = np.linalg.inv(cell)
    pos = A.get_positions()
    for comp in connected_components_indices(A):
        frac = pos[comp].mean(axis=0) @ inv_cell
        pos[comp] += -np.floor(frac) @ cell
    A.set_positions(pos)
    return A


def binary_to_idx(config, maxgroups):
    return sum([j * (maxgroups**i) for i, j in list(enumerate(reversed(config)))])


def select_configs(config, nsuper, maxgroups):
    """
    from the supercell config, return indices of each
    component primitive cell
    e.g.
    a 2x2x1 supercell of a Z=4 system (nsuper=4)
    would have configs of length 16.
    This function chunks them into 4
    sets of 4 and returns the integer representation
    (0-15) of each 'binary' number.
    So for example:
    >>> select_configs([1,1,1,1, 0,0,0,0, 1,1,0,0, 0,0,0,1], nsuper=4, maxgroups=2)
    [15, 0, 12, 1]

    """
    chunk = len(config) // nsuper
    # should equal Z
    return [binary_to_idx(c, maxgroups) for c in chunks(config, chunk)]


def chunks(lst, n, offset=0):
    """Yield successive n-sized chunks from lst."""
    if offset != 0:
        # rotate copy of list by offset
        temp = lst[offset:] + lst[:offset]
    else:
        temp = lst
    for i in np.arange(0, len(lst), n):
        yield temp[i : i + n]


class OrderedfromDisordered:
    def __init__(self, disordered_structure, symprec=1e-4, quiet=False):
        """
        Args:
            disordered_structure (DisorderedStructure)
            symprec (float): tolerance for symmetry finding
            quiet (bool): if True, suppresses warnings

        """
        self.quiet = quiet
        self.disordered_structure = disordered_structure

        self.symprec = symprec
        self.logger = logging.getLogger("orgdisord.enumerate")
        self.logger.debug("\n\n")
        self.logger.debug("------------------------------")
        self.logger.debug("--- ENUMERATING STRUCTURES ---")
        self.logger.debug("------------------------------")

    def _get_config_symbols_coords(self, config, asslabel):
        """
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

        """
        assembly = self.disordered_structure.get_assembly(asslabel)
        disorder_groups = assembly.disorder_groups
        # disorder_group_labels = sorted(group.label for group in disorder_groups)

        # a few aliases
        # how many disorder groups in this assembly?
        ndisordergroups = len(disorder_groups)

        sites = []
        symbols = []
        tags = []
        labels = []
        # loop over groups
        for igroup in range(ndisordergroups):
            group = disorder_groups[igroup]
            assert len(group.atoms) > 0
            group_symmops = group.symmetry_operations

            if len(group_symmops) == len(config):
                # choose the subset of symmops for this group given config
                group_config_symmops = [
                    group_symmops[i] for i, c in enumerate(config) if c == igroup
                ]
            else:
                # Zprime < 1 so config instead needs to index into subset of group symops
                # e.g. config [0,0,1,0] for would select those subsets of symops for this group
                # probably need to select the inverse for the next group...

                group_config_symmops = [
                    group_symmops[c][i] for i, c in enumerate(config)
                ]

            # tag = igroup + 1 if nassemblies == 1 else 1e2*(iassembly+1) + igroup + 1
            tag = group.tag
            [symbols_g, sites_g, tags_g, labels_g] = self._get_group_config_sites(
                tag, group, group_config_symmops
            )
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
        for op in group_config_symmops:
            # g = config[iops]
            # # Which group should we pick?
            # # if Zprime <1 or disorder_group_labels[igroup] < 0:
            # #     grouplabel = disorder_group_labels[igroup]
            # # else:
            # grouplabel = disorder_group_labels[igroup]
            # group = assembly.get_disorder_group(grouplabel)
            unique_positions = group.atoms.get_scaled_positions()
            unique_symbols = group.atoms.get_chemical_symbols()
            unique_labels = group.atoms.get_array("labels")
            # loop over site in each disorder group:
            for kind, pos in enumerate(unique_positions):
                # apply symmetry operation
                rot, trans = op
                site = np.mod(np.dot(rot, pos) + trans, 1.0)

                if not sites:
                    sites.append(site)
                    symbols.append(unique_symbols[kind])
                    tags.append(tag)
                    labels.append(unique_labels[kind])
                    continue
                t = site - sites
                mask = np.all(
                    (abs(t) < self.symprec) | (abs(abs(t) - 1.0) < self.symprec), axis=1
                )
                if np.any(mask):
                    inds = np.argwhere(mask).flatten()
                    if len(inds) > 1:
                        raise SpacegroupValueError(
                            "Found multiple equivalent sites!".format()
                        )
                    else:
                        pass
                else:
                    sites.append(site)
                    symbols.append(unique_symbols[kind])
                    tags.append(tag)
                    labels.append(unique_labels[kind])
        return [symbols, sites, tags, labels]

    def get_config(self, config, asslabel):
        """
        Return a fully ordered atoms object for the specified config.
        Sites are placed by the selected symmetry operations and wrapped into the
        primitive cell, so molecules may be split across periodic boundaries here.
        get_all_configs() makes the final structures whole again.
        """
        # set up base atoms object to which sites will be
        # added (just an empty box).
        atoms = Atoms(cell=self.disordered_structure.cell, pbc=True)

        # add in disordered sites, ordered according to config, (tag = 1e2*(iassembly+1) + igroup + 1)
        symbols, sites, tags, labels = self._get_config_symbols_coords(
            config, asslabel=asslabel
        )
        for site, symbol, tag in zip(sites, symbols, tags):
            atom = Atom(symbol=symbol, position=atoms.cell.T.dot(site), tag=tag)
            atoms.append(atom)
        atoms.set_array("labels", np.array(labels))

        return atoms

    def get_all_configs(self, exclude_ordered=False):
        """Generate all possible configs for the primitive cell and return them as
        ASE Atoms objects (pbc=False) with whole molecules, each in the same lattice
        image as the corresponding molecule of the first config (whose molecules are
        wrapped into the cell by centroid). Atoms may lie outside the cell.
        """
        cell = self.disordered_structure.cell
        Z = self.disordered_structure.Z
        # Zprime = self.cif.Zprime
        nassemblies = self.disordered_structure.get_number_of_assemblies()
        correlated_assemblies = self.disordered_structure.correlated_assemblies

        max_spacegroup_kinds = 0
        # set up base atoms object to which sites will be
        # added.
        if exclude_ordered:
            # just an empty box
            atoms = Atoms(cell=cell, pbc=True)
        else:
            atoms = self.disordered_structure.ordered_atoms.copy()
            if len(atoms) > 0:
                sg = self.disordered_structure.spacegroup
                max_spacegroup_kinds = sg.tag_sites(atoms.get_scaled_positions()).max()
                # give ordered sites a tag of 0
                atoms.set_tags(0)

        # what's the final number of configs we will generate?
        # TODO generalise!
        ngroups_per_assembly = (
            self.disordered_structure.get_number_of_disorder_groups_per_assembly()
        )
        nconfigs_per_assembly = [N**Z for N in ngroups_per_assembly]
        nconfigs = nconfigs_per_assembly[0]  # defaults to first ngroups

        if correlated_assemblies:
            if nassemblies > 1:
                self.logger.debug(
                    "Treating the assemblies as correlated. \n"
                    "This means that the same index from each assembly is chosen for each configuration."
                )
                if len(set(nconfigs_per_assembly)) != 1:
                    raise ValueError(
                        "Error: the number of configurations per assembly is not the same for all assemblies.\n"
                        "Please set correlated_assemblies to False."
                    )
        else:
            if nassemblies > 1:
                self.logger.debug(
                    "Treating the assemblies as uncorrelated "
                    "(i.e. the assemblies are independent) and we will therefore generate lots of structures!"
                )
                nconfigs = np.prod(nconfigs_per_assembly)
        self.logger.debug(f"Generating {nconfigs} configs in primitive cell.")
        if nconfigs > 1e4:
            self.logger.warning(
                f"Warning: {nconfigs} is a large number of configs. "
                "This may take a while!"
            )

        all_configs = []
        for assembly in self.disordered_structure.disorder_assemblies:
            asslabel = assembly.label
            # the disorder groups for this assembly
            disorder_groups = assembly.disorder_groups
            # normally ngroups = ndisorder groups
            ngroups = len(disorder_groups)
            group_nops = [len(group.symmetry_operations) for group in disorder_groups]
            self.logger.debug(f"Group n ops {group_nops}")

            if len(set(group_nops)) != 1:
                raise ValueError(
                    "Error: the number of symmetry operations per group "
                    f" is not the same for all groups in assembly {asslabel}:\n"
                    f"{group_nops}"
                )

            # but when Zprime is less than 1, we need to be more careful!
            # if Zprime < 1: # means more symmetry operations than we have Z
            #     # in this case the actual number of groups (1/Zprime) will be larger than the number of disorder groups
            #     assert 1/Zprime >= ngroups

            #     ngroups = int(1 / Zprime)
            #     self.logger.warning(f'Special case for nops != Z. (i.e. Zprime < 1)'
            #         f'nops = {nops_all} '
            #         f'Z = {Z}, '
            #         f'Zprime = {Zprime}\n'
            #         'This is less well-tested so proceed with caution.'
            #         'If you encounter any issues, please contact the developers.')

            self.logger.debug(f"Assembly {asslabel} has {ngroups} groups")

            if ngroups == 1:
                # Special disorder site case with only one group
                # -> each symmetry operation generated a new config.
                config_max = group_nops[0]
            else:
                config_max = ngroups

            # generate all possible configs for this assembly
            config_indices = itertools.product(list(range(config_max)), repeat=Z)
            assembly_configs = [
                self.get_config(config_idx, asslabel) for config_idx in config_indices
            ]
            # update max_spacegroup_kinds tags
            for config in assembly_configs:
                labels = config.get_array("labels")
                u, spacegroup_kinds = np.unique(labels, return_inverse=True)
                config.set_array(
                    "spacegroup_kinds", spacegroup_kinds + max_spacegroup_kinds + 1
                )
            self.logger.debug(
                f"Assembly {asslabel} has {len(assembly_configs)} configs"
            )
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

        # Make molecules whole and put every config's molecules in the same lattice
        # image as the first config's. Only the reference is wrapped (by centroid);
        # the others are aligned to it by nearest image and NOT re-wrapped. Wrapping
        # each config by its own centroid would put a disordered molecule whose
        # centroid sits near a cell face one lattice vector away in different
        # configs, so mixed supercells would overlap that molecule with itself.
        # Done once per config here so get_supercell_config can simply copy,
        # translate and tile these images. Note the images are non-periodic
        # (pbc=False): reload_as_molecular_crystal relies on that so its
        # atoms.wrap() does not tear the molecules apart again.
        if not all_configs:
            return all_configs
        reference = wrap_components_into_cell(make_whole_molecule(all_configs[0]))
        cell = np.asarray(reference.get_cell())
        ref_comps = connected_components_indices(reference)
        ref_centroids = component_centroids(reference, ref_comps, mass_weighted=False)

        aligned_configs = []
        for cfg in all_configs:
            cfg = make_whole_molecule(cfg)
            cfg = align_components_to_reference(cfg, ref_centroids, cell)
            aligned_configs.append(cfg)
        all_configs = aligned_configs

        self.logger.debug(f"Generated {len(all_configs)} primitive configs")
        return all_configs

    def get_supercell_configs(
        self,
        supercell=None,
        target_size=None,
        target_shape="sc",
        score_key="length",
        maxiters=5000,
        exclude_ordered=False,
        random_configs=False,
        return_configs=False,
        fix_ratio=False,
        ratio=None,
        ratio_tolerance=0.01,
    ):
        """
        loop over supercell cells,
        add in one of the ndisordergroups^Z configurations per cell
        if return_configs is True, return the configs as well as the supercells
        if molecular_crystal is True, then the structure is reloaded as a molecular crystal
             -- trying to keep the pieces in tact

        Args:
            supercell (int, list, tuple, or 3x3 array): supercell dimensions e.g. [2,2,2],
                            or 3x3 integer transformation matrix P (supercell cell = P @ primitive cell),
                            or an integer number of unit cells (same as target_size).
                            Only one of supercell and target_size may be given.
            target_size (int): number of primitive cells. The supercell whose shape is closest to
                            target_shape is used.
            target_shape (str): 'sc' (simple cubic) or 'fcc'. Shape used when searching for a target_size supercell.
            score_key (str): 'length' or 'metric'. Score used when searching for a target_size supercell.
            maxiters (int): maximum number of structures to generate. If random_configs is True,
                            then this is the number of random configs to generate.
                            If random_configs is False and the total number of configs is less than maxiters,
                            then all configs will be generated.
                            Default is 5000.
            exclude_ordered (bool): if True, exclude the ordered part of supercell in the returned structures.
            random_configs (bool): if True, choose the configs randomly from all allowed configs.
            return_configs (bool): if True, return the configs (e.g. (0,1,0,1)) as well as the supercells
            fix_ratio (bool): if True, then the ratio(s) of the disorder components is/are fixed. Will use the group
                              occupancies to determine the ratio(s) unless ratio is specified.
            ratio (float): Ratio(s) of the disorder components. If None, then the ratio(s) will be determined from the
                            group occupancies. If ratio is specified, then the group occupancies will be ignored.
                            If ratio is a float, then the ratio of the first group will be fixed to this value.
                            If ratio is a list, then the ratio of each group will be fixed to the corresponding value.
                            If ratio is a list of lists, then the ratio of each group will be fixed to the corresponding
                            value in the inner list. The outer list corresponds to the assembly.
                            Default is None.
            ratio_tolerance (float): Tolerance for the ratio(s) of the disorder components. Default is 0.01. See np.isclose atol parameter.

        Returns:
            supercells (list): list of supercells
            configs (list): list of configs (optional)
        """
        # pre-compute all primitive cells
        self.logger.debug("Pre-computing all primitive configurations...")
        self.ratio_tolerance = ratio_tolerance
        if ratio and not fix_ratio:
            self.logger.warning(
                "Warning! ratio is specified but fix_ratio is False. Ignoring ratio."
            )

        ########################################
        # Get all primitive configs
        self.images = self.get_all_configs(exclude_ordered)
        nimages_prim = len(self.images)
        if self.disordered_structure.molecular_crystal:
            self.logger.debug(f"Reloading {nimages_prim} images as molecular crystals")
            cheap_method = False
            if nimages_prim > 256:
                # this seems to be get pretty slow for large numbers of images
                self.logger.warning(
                    "Warning: reloading molecular crystals is slow for this many images.\n"
                    "Consider skipping this step with the --not_molecular_crystal flag!\n"
                    "We will use a cheaper method to reload the images as molecular crystals instead\n"
                    " -- assuming the connectivity is the same for each image"
                )
                cheap_method = True

            self.images = reload_as_molecular_crystal(self.images, cheap=cheap_method)
        ########################################

        # some aliases
        Z = self.disordered_structure.Z
        self.cell = self.disordered_structure.cell

        if target_size is not None and supercell is not None:
            raise ValueError("Only one of supercell and target_size may be given.")
        if target_size is not None:
            spec = target_size
        elif supercell is not None:
            spec = supercell
        else:
            spec = [1, 1, 1]

        # Determine transformation_matrix P
        if isinstance(spec, (int, np.integer)) and not isinstance(spec, bool):
            if spec < 1:
                raise ValueError(f"target_size must be >= 1. Got: {spec}")
            msg = (
                f"Could not find a supercell of {spec} primitive cells with "
                f"target_shape={target_shape!r}; try a different target size or shape."
            )
            try:
                P = find_optimal_cell_shape(
                    self.cell,
                    target_size=int(spec),
                    target_shape=target_shape,
                    score_key=score_key,
                )
            except ValueError as err:
                raise ValueError(msg) from err
            if P is None:
                raise ValueError(msg)
            P = np.asarray(P, dtype=int)
        else:
            arr = np.asarray(spec)
            if arr.shape not in [(3,), (3, 3)]:
                raise ValueError(
                    f"supercell must be an integer target size, a length-3 sequence, or a 3x3 matrix. Got: {spec}"
                )
            if not np.allclose(arr, np.round(arr)):
                raise ValueError(f"supercell must contain integers only. Got: {spec}")
            arr = np.round(arr).astype(int)
            P = np.diag(arr) if arr.ndim == 1 else arr

        det = int(round(np.linalg.det(P)))
        if det == 0:
            raise ValueError("Supercell transformation matrix has zero determinant.")
        if det < 0:
            P = -P

        self.transformation_matrix = P
        self.nsuper = abs(det)
        self.supercell_matrix = P @ np.asarray(self.cell)
        self.lattice_points_frac = lattice_points_in_supercell(P)
        self.lattice_points_cart = self.lattice_points_frac @ self.supercell_matrix

        # we previously just took ndisordergroups = self.cif.ndisordergroups, but
        # this gets tricky when we have multiple assemblies etc. Safer to just take this:
        self.ndisordergroups = int(np.round(np.exp(np.log(nimages_prim) / Z)))
        # total number of configs given this supercell:
        ncombinations = self.ndisordergroups ** (Z * self.nsuper)
        self.logger.debug(
            f"""
        Found {nimages_prim} primitive configurations,
              {ncombinations} configs for this supercell (Z={Z}, nsuper={self.nsuper}, matrix=\n{P}),
              {self.ndisordergroups} ndisordergroups"""
        )

        options = list(range(self.ndisordergroups))
        size = Z * self.nsuper
        if random_configs:
            all_combinations = random_product(options, repeat=size)
        else:
            # this iterator generates all
            # possible lists of e.g. 0s and 1s etc of length Z*nsuper
            all_combinations = itertools.product(options, repeat=size)

        if fix_ratio:
            # we will generate configs with a fixed ratio of ordered to disordered
            # we will do this by generating a random list of 0s and 1s
            # and then sorting them by the number of 1s they have
            # this should give us a list of configs with a roughly fixed ratio
            # of ordered to disordered
            # we will then take the first n_configs of these
            self.ratio = self._process_ratio(ratio)
            self.logger.info(f"Fixing ratio to {self.ratio}")

            # filter all_combinations by the ratio
            # all_combinations = self._ratio_filter(all_combinations, ratio)

            # apply filter to all_combinations
            all_combinations = filter(self._ratio_filter, all_combinations)

        # how many configs to actually generate:
        if random_configs:
            # just take maxiters
            n_configs = maxiters
        else:
            # take whichever is smallest betwee
            # maxiters and ncombinations
            n_configs = min([maxiters, ncombinations])
            if n_configs < ncombinations:
                self.logger.warning(
                    f"""
                Warning: maxiters ({maxiters}) is smaller than the total number of configs ({ncombinations})
                We will only generate {n_configs} configs.
                You might want to consider increasing maxiters or using the --random flag"""
                )

        all_supercells = []
        all_configs = []
        if random_configs:
            self.logger.info(
                f"Generating {n_configs} random configurations in a {spec} supercell:"
            )
        else:
            self.logger.info(
                f"Generating {n_configs} out of the {ncombinations} possible configurations in the {spec} supercell:"
            )

        for _ in tqdm(range(n_configs), disable=self.quiet):
            config = next(all_combinations, False)
            if not config:
                self.logger.debug("No more configs left!")
                break
            # make ncells copies of the relevant config
            self.logger.debug(f"         {config}")
            all_supercells.append(self.get_supercell_config(config))
            # TODO
            # check for overlapping atoms
            # cutoff = 0.25
            # dists = pdist(supercell_atoms.get_positions(), 'sqeuclidean')
            # if (dists < cutoff**2).any():
            #     continue
            # else:
            all_configs.append(config)

        if fix_ratio:
            self.logger.info(
                f"Only {len(all_configs)}/{n_configs} configs had the specified ratio of {self.ratio}"
            )
        if return_configs:
            return all_supercells, all_configs
        return all_supercells

    def _process_ratio(self, ratio):
        """
        Process the ratio argument
        """
        if ratio is None:
            # use the group occupancies
            ratio = [
                [g.occupancy for g in ass.disorder_groups]
                for ass in self.disordered_structure.disorder_assemblies
            ]
        elif isinstance(ratio, float):
            # use the specified ratio
            # if only one float, then assume one assembly, two groups
            ratio = [[ratio, 1 - ratio]]
        elif isinstance(ratio, (list, tuple)):
            if isinstance(ratio[0], float):
                # use the specified ratios. Assume just one assembly
                ratio = [ratio]
            elif isinstance(ratio[0], int):
                # assume this is a list of integers, with the occupancies being
                # computed from the ratio.
                # use the specified ratios. Assume just one assembly
                ratio = [ratio / ratio.sum()]

            elif isinstance(ratio[0], list):
                # use the specified ratios as provided
                ratio = ratio
            else:
                raise ValueError("ratio must be a float, int or list")
        else:
            raise ValueError("ratio must be a float, int or list")
        # convert to np arrays
        ratio = np.array(ratio)
        return ratio

    def _ratio_filter(self, config):
        """
        Returns True if the config matches the ratio

        TODO: test/make this work for multiple assemblies!
        """

        computed_ratio = np.bincount(config, minlength=self.ndisordergroups)
        computed_ratio = computed_ratio / computed_ratio.sum()

        # flatten ratio array -- this this might not work for multiple assemblies
        ratio = self.ratio
        if len(ratio) > 1:
            raise NotImplementedError(
                "ratio filtering not implemented for multiple assemblies"
            )

        ratio = ratio[0]

        # check if the computed ratio is within self.ratio_tolerance of the desired ratio
        if np.allclose(ratio, computed_ratio, atol=self.ratio_tolerance):
            return True

        return False

        # return filtered_combination

    def get_supercell_config(self, config) -> Atoms:
        """
        get a supercell config

        Args:
            config (np.array): array of integers of length Z*nsuper
                                (where Z is the number of sites in the primitive cell
                                and nsuper is the number of tiled primitive cells)

        """
        supercell_atoms = Atoms(cell=self.supercell_matrix, pbc=True)
        for icell, c in enumerate(
            select_configs(config, nsuper=self.nsuper, maxgroups=self.ndisordergroups)
        ):
            # make a copy of prim config c
            # (primitive images are already whole and wrapped by get_all_configs)
            temp = self.images[c].copy()
            # translate primitive copy to its lattice point in the supercell
            temp.translate(self.lattice_points_cart[icell])
            # add this block into the supercell
            supercell_atoms.extend(temp)
        return supercell_atoms
