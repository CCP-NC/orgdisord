'''
Contains the class DisorderedStructure, which is used to represent a
disordered structure.
'''
from logging import warning
from sodorg_renewal.utils import molecule_collide
from dataclasses import dataclass
from typing import List, Tuple, Union
from ase import Atoms, Atom
from ase.cell import Cell
from ase.spacegroup import Spacegroup
from ase.geometry import find_mic
import numpy as np
import warnings

@dataclass
class DisorderGroup:
    '''
    A class for representing a disorder group.
    This contains the fractional coordinates of the atoms in the group,
    their symbols and their symmetry operations.
    '''
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
    special_disorder_symmetry: bool = False
    tag: int = None

    def __post_init__(self):
        # check that atoms object has some custom attributes
        assert hasattr(self.atoms, 'cell')
        assert self.atoms.has('labels') and self.atoms.has('occupancies')
        # Get subset of group symmetry operations
        # self.symmetry_operations = self.get_group_symmetry_operations()
        # if label has a - in it, then it's a special disorder group
        if '-' in self.label:
            self.special_disorder_symmetry = True
        else:
            self.special_disorder_symmetry = False
    
    def get_group_symmetry_operations(self,
                                      sg: Spacegroup,
                                      Z: int,
                                      ) -> List[Tuple[np.ndarray, np.ndarray]]:
        '''
        Get the symmetry operations of the group.

        Args:
            sg: ASE spacegroup object for the crystal in question.
            Z: Number of molecular units in the primitive cell.
            

        Returns:
            List of symmetry operations, each of which is a tuple of
            (rotation matrix, translation vector).
            For Zprime < 1, we chunk the returned list into 1/Zprime groups of
            symmetry operations.
        '''

        # list of all symmetry operations
        all_symmetry_operations = sg.get_symop()


        # HACK for now:
        if self.special_disorder_symmetry:
            # assert Z != len(all_symmetry_operations)
            nsubgroups = len(all_symmetry_operations) // Z
            ## chunk into subgroups
            group_symmetry_operations = [all_symmetry_operations[i::nsubgroups] 
                                            for i in range(nsubgroups)]
        else:
            ## otherwise use all symmetry operations
            group_symmetry_operations = all_symmetry_operations
        return group_symmetry_operations



    def get_labels(self) -> List[str]:
        '''
        Get the labels of the atoms in the group.
        '''
        return self.atoms.get_array('labels')
    def __repr__(self) -> str:
        repr = f'Disorder group: {self.label} '
        if self.special_disorder_symmetry:
            repr += 'has special disorder symmetry and contains '
            repr += f'{[len(symops) for symops in self.symmetry_operations]} symm. ops. and '
        else:
            repr += f'contains {len(self.symmetry_operations)} symm. ops. and '
        repr += f'{len(self.atoms)} sites:\n'
        site_labels = self.atoms.get_array('labels')
        occupancies = self.atoms.get_array('occupancies')
        for i, atom in enumerate(self.atoms):
            repr += f"\tlabel: {site_labels[i]: >8}, "
            repr += f"species: {atom.symbol: >4}, "
            repr += f"occupancy: {occupancies[i]: >5.2f}\n"
        return repr
    
    


@dataclass
class DisorderAssembly:
    '''
    A class for representing a disorder assembly.
    An assembly is a list of DisorderGroups
    '''
    label: str
    disorder_groups: List[DisorderGroup]
    ngroups: int = None
    tag: int = None

    def __post_init__(self):
        '''
        Post-initialization.
        '''
        self.ngroups = len(self.disorder_groups)
        assert self.ngroups > 0


    def __repr__(self) -> str:
        repr = f'Disorder assembly: {self.label}\n'
        repr += 'Contains the following groups:\n'
        for group in self.disorder_groups:
            repr += f'{group}\n'
        return repr

    def get_disorder_group(self, label) -> DisorderGroup:
        '''
        Get the disorder group of the assembly with 
        the given label.
        '''
        for group in self.disorder_groups:
            if group.label == label:
                return group
        raise ValueError(f'No disorder group with label {label}')


@dataclass
class DisorderedStructure:
    '''
    A class for representing a disordered structure.
    '''
    # The list of atoms in the structure.
    # these have a symbol and a position
    ordered_atoms: Atoms
    cell: Cell
    Z: int
    # The spacegroup of the structure.
    spacegroup: Spacegroup
    disorder_assemblies: List[DisorderAssembly]
    # Should we treat the system as a molecular crystal?
    molecular_crystal: bool = True
    # are the assemblies in the structure correlated with each other?
    correlated_assemblies: bool = False
    
    def __post_init__(self):
        '''
        Post-initialization.
        '''
        # print the assemblies
        print('Disordered structure:')
        for assembly in self.disorder_assemblies:
            print(assembly)

        # warning if any assemblies has no groups
        for assembly in self.disorder_assemblies:
            if assembly.ngroups == 0:
                warning.warn(f'Warning: assembly {assembly.label} has no groups!')

        # drop empty assemblies
        self.disorder_assemblies = [assembly for assembly in self.disorder_assemblies
                                    if assembly.ngroups > 0]

        # if the assemblies are correlated, then we need to make sure that
        # the number groups in each assembly is the same
        if self.correlated_assemblies:
            ngroups = self.disorder_assemblies[0].ngroups
            for assembly in self.disorder_assemblies:
                assert assembly.ngroups == ngroups


    def get_number_of_assemblies(self) -> int:
        '''
        Get the number of disorder assemblies in the structure.
        '''
        return len(self.disorder_assemblies)

    def get_number_of_disorder_groups_per_assembly(self) -> List[int]:
        '''
        Get the number of disorder groups in each assembly in the structure.
        '''
        return [assembly.ngroups for assembly in self.disorder_assemblies]
    def __repr__(self) -> str:
        rep = 'DisorderedStructure:\n'
        rep += '\tordered_atoms: {}\n'.format(self.ordered_atoms)
        rep += '\tcell: {}\n'.format(self.cell)
        rep += '\tZ: {}\n'.format(self.Z)
        rep += '\tspacegroup: {}\n'.format(self.spacegroup)
        rep += '\tdisorder_assemblies: \n\t{}\n'.format(self.disorder_assemblies)
        return rep
    
    def get_assembly(self, label) -> DisorderAssembly:
        '''
        Get the disorder assembly with 
        the given label.
        '''
        for assembly in self.disorder_assemblies:
            if assembly.label == label:
                return assembly
        raise ValueError(f'No disorder assembly with label {label}')
