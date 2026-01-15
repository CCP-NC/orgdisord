"""Module for representing disordered crystal structures.

This module provides classes to represent and manipulate disordered crystal structures,
with support for multiple disorder assemblies and groups. It handles both molecular
and non-molecular crystals, tracking atomic positions, symmetry operations, and
occupancies.

Classes:
    DisorderGroup: Represents a group of disordered atoms sharing occupancy
    DisorderAssembly: Collection of related disorder groups
    DisorderedStructure: Complete representation of a disordered crystal structure
"""
from logging import warning
from orgdisord.utils import get_unique_atoms, standardise_cell
from dataclasses import dataclass
from typing import List, Tuple, Optional
from ase.spacegroup import Spacegroup, get_spacegroup
from ase.geometry import find_mic
from ase import Atoms
import numpy as np
import warnings


@dataclass
class DisorderGroup:
    """A group of atoms that represent one possible arrangement in a disorder.
    
    This class represents a set of atoms that share the same occupancy in a
    disordered structure. The atoms in a group are related by symmetry operations
    and together represent one possible configuration of the disorder.
    
    Attributes:
        label (str): Unique identifier for this disorder group
        atoms (Atoms): ASE Atoms object containing the atomic positions and properties
        symmetry_operations (List[Tuple[np.ndarray, np.ndarray]]): Symmetry operations 
            defining equivalent positions
        occupancy (Optional[float]): Site occupancy factor for this group
        special_disorder_symmetry (bool): Whether this group has special symmetry handling
        tag (Optional[int]): Numeric identifier used in visualization/analysis
    
    Requirements:
        The atoms object must have:
            - A defined unit cell
            - A "labels" array containing atom labels
            - Optional "occupancies" array if occupancy not specified
    """

    label: str
    # ASE atoms -- must have:
    #       cell
    #       positions
    #       symbols
    #       labels array
    atoms: Atoms
    symmetry_operations: List[Tuple[np.ndarray, np.ndarray]]
    occupancy: Optional[float] = None
    special_disorder_symmetry: bool = False
    tag: Optional[int] = None

    def __post_init__(self) -> None:
        """Validate and initialize the disorder group."""
        if not hasattr(self.atoms, "cell"):
            raise ValueError("Atoms object must have a defined unit cell")
        if not self.atoms.has("labels"):
            raise ValueError("Atoms object must have a 'labels' array")

        self._process_occupancies()
        
        # Set special symmetry flag based on label format
        self.special_disorder_symmetry = "-" in self.label

    def _process_occupancies(self) -> None:
        """
        Process the occupancies of the atoms in the group.
        """
        if self.occupancy:
            if not self.atoms.has("occupancies"):
                self.atoms.new_array(
                    "occupancies", np.ones(len(self.atoms)) * self.occupancy
                )
            else:
                # make sure the occupancies are the same as the group occupancy
                if not np.allclose(self.atoms.get_array("occupancies"), self.occupancy):
                    warnings.warn(
                        "Occupancies of atoms in group do not match group occupancy."
                    )
        else:
            assert self.atoms.has("occupancies")
            self.occupancy = self.atoms.get_array("occupancies").mean()

    def get_group_symmetry_operations(
        self,
        sg: Spacegroup,
        Z: int,
    ) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Get the symmetry operations of the group.

        Args:
            sg: ASE spacegroup object for the crystal in question.
            Z: Number of molecular units in the primitive cell.


        Returns:
            List of symmetry operations, each of which is a tuple of
            (rotation matrix, translation vector).
            For Zprime < 1, we chunk the returned list into 1/Zprime groups of
            symmetry operations.
        """

        # list of all symmetry operations
        all_symmetry_operations = sg.get_symop()

        # HACK for now:
        if self.special_disorder_symmetry:
            # assert Z != len(all_symmetry_operations)
            nsubgroups = len(all_symmetry_operations) // Z
            print(f"Found {nsubgroups} subgroups")
            ## chunk into subgroups
            group_symmetry_operations = [
                all_symmetry_operations[i::nsubgroups] for i in range(nsubgroups)
            ]
            print(f"Group symmetry operations: {group_symmetry_operations}")
        else:
            ## otherwise use all symmetry operations
            group_symmetry_operations = all_symmetry_operations
        return group_symmetry_operations

    def get_labels(self) -> List[str]:
        """
        Get the labels of the atoms in the group.
        """
        return self.atoms.get_array("labels")

    def __repr__(self) -> str:
        repr = f"Disorder group: {self.label} "
        if self.special_disorder_symmetry:
            repr += "has special disorder symmetry and contains "
            repr += f"{[len(symops) for symops in self.symmetry_operations]} symm. ops. and "
        else:
            repr += f"contains {len(self.symmetry_operations)} symm. ops., "
        repr += f"a group occupancy of {self.occupancy} and "
        repr += f"{len(self.atoms)} sites:\n"
        site_labels = self.get_labels()
        occupancies = self.atoms.get_array("occupancies")
        for i, atom in enumerate(self.atoms):
            repr += f"\tlabel: {site_labels[i]: >8}, "
            repr += f"species: {atom.symbol: >4}, "
            repr += f"occupancy: {occupancies[i]: >5.2f}\n"
        return repr


@dataclass
class DisorderAssembly:
    """A collection of disorder groups that together represent a specific 
    arrangement of a disordered structure.

    Attributes:
        label (str): Unique identifier for this disorder assembly
        disorder_groups (List[DisorderGroup]): List of disorder groups in this assembly
        ngroups (int): Number of disorder groups in the assembly
        tag (Optional[int]): Numeric identifier used in visualization/analysis
    """

    label: str
    disorder_groups: List[DisorderGroup]
    ngroups: int = None
    tag: Optional[int] = None

    def __post_init__(self):
        """
        Post-initialization.
        """
        self.ngroups = len(self.disorder_groups)
        assert self.ngroups > 0

    def __repr__(self) -> str:
        repr = f"Disorder assembly: {self.label}\n"
        repr += "Contains the following groups:\n"
        for group in self.disorder_groups:
            repr += f"{group}\n"
        return repr

    def get_disorder_group(self, label) -> DisorderGroup:
        """
        Get the disorder group of the assembly with
        the given label.
        """
        for group in self.disorder_groups:
            if group.label == label:
                return group
        raise ValueError(f"No disorder group with label {label}")


@dataclass
class DisorderedStructure:
    """A disordered structure consisting of multiple disorder assemblies.

    This class represents a crystal structure that exhibits disorder, which can
    be due to the presence of different atomic or molecular species at certain
    positions.

    Attributes:
        ordered_atoms (Atoms): The atoms in the structure, in a standard orientation
        Z (int): Number of molecular units in the primitive cell
        spacegroup (Spacegroup): The spacegroup of the structure
        disorder_assemblies (List[DisorderAssembly]): List of disorder assemblies in 
            the structure
        molecular_crystal (bool): Should the system be treated as a molecular crystal?
        correlated_assemblies (bool): Are the assemblies in the structure correlated 
            with each other?
    """

    # The list of atoms in the structure.
    # these have a symbol and a position
    ordered_atoms: Atoms
    Z: int
    # The spacegroup of the structure.
    spacegroup: Spacegroup
    disorder_assemblies: List[DisorderAssembly]
    # Should we treat the system as a molecular crystal?
    molecular_crystal: bool = True
    # are the assemblies in the structure correlated with each other?
    correlated_assemblies: bool = False

    def __post_init__(self):
        """
        Post-initialization.
        """
        # get cell
        self.cell = self.ordered_atoms.cell

        # print the assemblies
        print("Disordered structure:")
        for assembly in self.disorder_assemblies:
            print(assembly)

        # warning if any assemblies has no groups
        for assembly in self.disorder_assemblies:
            if assembly.ngroups == 0:
                warning.warn(f"Warning: assembly {assembly.label} has no groups!")

        # drop empty assemblies
        self.disorder_assemblies = [
            assembly for assembly in self.disorder_assemblies if assembly.ngroups > 0
        ]

        # if the assemblies are correlated, then we need to make sure that
        # the number groups in each assembly is the same
        if self.correlated_assemblies:
            ngroups = self.disorder_assemblies[0].ngroups
            for assembly in self.disorder_assemblies:
                assert assembly.ngroups == ngroups

    def get_number_of_assemblies(self) -> int:
        """
        Get the number of disorder assemblies in the structure.
        """
        return len(self.disorder_assemblies)

    def get_number_of_disorder_groups_per_assembly(self) -> List[int]:
        """
        Get the number of disorder groups in each assembly in the structure.
        """
        return [assembly.ngroups for assembly in self.disorder_assemblies]

    def __repr__(self) -> str:
        rep = "DisorderedStructure:\n"
        rep += "\tordered_atoms: {}\n".format(self.ordered_atoms)
        rep += "\tcell: {}\n".format(self.cell)
        rep += "\tZ: {}\n".format(self.Z)
        rep += "\tspacegroup: {}\n".format(self.spacegroup)
        rep += "\tdisorder_assemblies: \n\t{}\n".format(self.disorder_assemblies)
        return rep

    def get_assembly(self, label) -> DisorderAssembly:
        """
        Get the disorder assembly with
        the given label.
        """
        for assembly in self.disorder_assemblies:
            if assembly.label == label:
                return assembly
        raise ValueError(f"No disorder assembly with label {label}")


def from_disorder_components(
    atoms_maj,
    atoms_min,
    tolerance=1e-4,
    symprec=1e-4,
    ratio=None,
    group_occupancies=None,
):
    """Create a disordered structure from two components.

    Args:
        atoms_maj (Atoms): ASE Atoms object representing the major component
        atoms_min (Atoms): ASE Atoms object representing the minor component
        tolerance (float, optional): Tolerance for comparing coordinates. Defaults to 1e-4.
        symprec (float, optional): Tolerance for symmetry operations. Defaults to 1e-4.
        ratio (float, optional): Ratio of major to minor components. Defaults to None.
        group_occupancies (list, optional): List of occupancies for each disorder group.
            Defaults to None.

    Returns:
        DisorderedStructure: A disordered structure object
    """
    # # Add labels array if not present
    # if not atoms_maj.has('labels'):
    #     atoms_maj.set_array('labels', [f"{sym}{i}" for i, sym in enumerate(atoms_maj.get_chemical_symbols())])
    # if not atoms_min.has('labels'):
    #     atoms_min.set_array('labels', [f"{sym}{i}" for i, sym in enumerate(atoms_min.get_chemical_symbols())])

    # # Get spacegroup and standardise cell
    # sg = get_spacegroup(atoms_maj, symprec=symprec)
    atoms_maj = standardise_cell(atoms_maj, symprec=symprec)
    atoms_min = standardise_cell(atoms_min, symprec=symprec)

    # set the group occupancies
    group_occupancies = group_occupancies or [0.5, 0.5]
    # if ratio is set, use that instead
    if ratio:
        group_occupancies = [ratio, 1 - ratio]

    atoms_ref = atoms_maj.copy()
    cell = atoms_ref.cell

    ## Since the atoms objects are already in the correct order,
    ## we can just use the diffence in positions to determine the non-changing atoms
    ## If this were not the case, we'd need to reorder the atoms first!
    maj_pos = atoms_maj.get_positions()
    min_pos = atoms_min.get_positions()

    # get minimum image convention distances between maj and min atoms
    # maj_min_dist = np.linalg.norm(maj_pos - min_pos, axis = 1)
    _, maj_min_dist = find_mic(maj_pos - min_pos, cell)

    # get the indices of the atoms that are different
    diff_idx = np.where(maj_min_dist > tolerance)[0]

    # Or check where the symbol list is different
    diff_idx_symbols = np.where(
        [
            maj != min
            for maj, min in zip(
                atoms_maj.get_chemical_symbols(), atoms_min.get_chemical_symbols()
            )
        ]
    )[0]

    # combine the two
    diff_idx = np.unique(np.concatenate([diff_idx, diff_idx_symbols]))

    maj_atoms = atoms_maj[diff_idx].copy()
    min_atoms = atoms_min[diff_idx].copy()

    # ordered atoms are those that are the same in both structures
    ordered_atoms = atoms_maj[
        [i for i in range(len(atoms_maj)) if i not in diff_idx]
    ].copy()

    print(
        f"""
    The major and minor components have {len(atoms_maj)} and {len(atoms_min)} atoms respectively.
    The number of atoms that are different is {len(diff_idx)}.
    The number of atoms that are the same is {len(ordered_atoms)}.
    The number of atoms in the ordered structure is {len(ordered_atoms)}.
    The number of atoms in the major component is {len(maj_atoms)}.
    The number of atoms in the minor component is {len(min_atoms)}.
    """    )    

    sg = get_spacegroup(ordered_atoms, symprec=symprec)

    # get just the asymmetric cell atoms
    maj_atoms, maj_Z = get_unique_atoms(maj_atoms, get_spacegroup(atoms_maj, symprec=symprec), symprec=symprec)
    min_atoms, min_Z = get_unique_atoms(min_atoms, get_spacegroup(atoms_min, symprec=symprec), symprec=symprec)

    if maj_Z != min_Z:
        raise ValueError(
            f"""
        The major and minor components seem to have different number of Z units!
        maj: {maj_Z} vs min: {min_Z}"""
        )
    Z = maj_Z

    # set occupancies
    ordered_atoms.set_array("occupancies", np.ones(len(ordered_atoms)))
    maj_atoms.set_array("occupancies", group_occupancies[0] * np.ones(len(maj_atoms)))
    min_atoms.set_array("occupancies", group_occupancies[1] * np.ones(len(min_atoms)))

    # symmetry operations of reference structure
    sg = get_spacegroup(ordered_atoms, symprec=symprec)
    symops = sg.get_symop()
    nsymops = len(symops)
    ## even indexed symops#
    # TODO THIS NEEDS TO BE GENERALISED!
    if nsymops != Z:
        warnings.warn(
            f"""
        The number of symmetry operations ({nsymops}) 
        is not equal to the number of molecular units (Z = {Z}).
        This is not well covered by the code at the moment.
        Proceed with caution!"""
        )
        symops = symops[0 :: nsymops // Z]

    # make list of DisorderGroup objects to pass into DisorderAssembly
    disorder_groups = []
    for i, atoms in enumerate([maj_atoms, min_atoms]):
        disorder_group = DisorderGroup(
            label=str(i + 1),
            atoms=atoms,
            symmetry_operations=symops,
            occupancy=group_occupancies[i],
            tag=i + 1,
        )
        disorder_groups.append(disorder_group)

    # create disorder assembly
    da = DisorderAssembly(
        label="A",
        disorder_groups=disorder_groups,
        tag=0,
    )

    # create disordered structure
    ds = DisorderedStructure(
        ordered_atoms=ordered_atoms,
        Z=Z,
        spacegroup=sg,
        disorder_assemblies=[da],
        molecular_crystal=True,
    )

    return ds
