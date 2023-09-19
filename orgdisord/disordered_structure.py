"""
Contains the class DisorderedStructure, which is used to represent a
disordered structure.
"""
from logging import warning
from orgdisord.utils import molecule_collide, get_unique_atoms, standardise_cell
from dataclasses import dataclass
from typing import List, Tuple, Union
from ase import Atoms, Atom
from ase.cell import Cell
from ase.spacegroup import Spacegroup, get_spacegroup
from ase.geometry import find_mic
import numpy as np
import warnings


@dataclass
class DisorderGroup:
    """
    A class for representing a disorder group.
    This contains the fractional coordinates of the atoms in the group,
    their symbols and their symmetry operations.
    """

    label: str
    # ASE atoms -- must have:
    #       cell
    #       positions
    #       symbols
    #       occupancies array
    #       labels array
    atoms: Atoms
    # list of symmetry operations
    symmetry_operations: List[Tuple[np.ndarray, np.ndarray]]
    occupancy: float = None
    special_disorder_symmetry: bool = False
    tag: int = None

    def __post_init__(self):
        # check that atoms object has some custom attributes
        assert hasattr(self.atoms, "cell")
        assert self.atoms.has("labels")

        self._process_occupancies()
        # if self.occupancy:
        #     if not self.atoms.has('occupancies'):
        #         self.atoms.new_array('occupancies', [self.occupancy] * len(self.atoms))
        # else:
        #     assert self.atoms.has('occupancies')
        #     self.occupancy = self.atoms.get_array('occupancies').mean()
        # Get subset of group symmetry operations
        # self.symmetry_operations = self.get_group_symmetry_operations()
        # if label has a - in it, then it's a special disorder group
        if "-" in self.label:
            self.special_disorder_symmetry = True
        else:
            self.special_disorder_symmetry = False

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
            ## chunk into subgroups
            group_symmetry_operations = [
                all_symmetry_operations[i::nsubgroups] for i in range(nsubgroups)
            ]
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
    """
    A class for representing a disorder assembly.
    An assembly is a list of DisorderGroups
    """

    label: str
    disorder_groups: List[DisorderGroup]
    ngroups: int = None
    tag: int = None

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
    """
    A class for representing a disordered structure.
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
    tolerance=1e-2,
    symprec=1e-3,
    ratio=None,
    group_occupancies=None,
):
    """
    Given the major and minor components of a disordered structure,
    (i.e. two fully ordered structures representing the two disordered structures)
    return a DisorderedStructure object.

    Limitations: At the moment it only works for 1 assembly and two groups.

    We also need to generalise the choice of symmetry operation subset!

    Args:
        atoms_maj (Atoms): The P1, ordered major component of the disordered structure.
        atoms_min (Atoms): The P1, ordered minor component of the disordered structure.
        tolerance (float): The tolerance for determining if two atoms in the same position. (Å)
        symprec   (float): Tolerance used for finding the spacegroup using spglib).
        ratio     (float): The ratio of the occupancies of the two groups.
                           For example 0.75 means that the major component
                           has an occupancy of 0.75.
                           If None, then the ratio is determined from the
                           occupancies of the atoms in the two structures
                           or the ``group_occupancies``.
        group_occupancies (list): The occupancies of the two groups.
                                  If None, then the occupancies are
                                  determined from the occupancies of
                                  the atoms in the two structures.

    Returns:
        DisorderedStructure: The disordered structure.
    """
    # convert the structures to the standardised cell for this spacegroup
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

    maj_atoms = atoms_maj[diff_idx].copy()
    min_atoms = atoms_min[diff_idx].copy()

    # ordered atoms are those that are the same in both structures
    ordered_atoms = atoms_maj[
        [i for i in range(len(atoms_maj)) if i not in diff_idx]
    ].copy()

    sg = get_spacegroup(ordered_atoms, symprec=symprec)

    # get just the asymmetric cell atoms
    maj_atoms, maj_Z = get_unique_atoms(maj_atoms, sg, symprec=symprec)
    min_atoms, min_Z = get_unique_atoms(min_atoms, sg, symprec=symprec)

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

    # return disordered structure
    return ds
