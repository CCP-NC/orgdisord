import numpy as np
from ase import Atoms, Atom
from soprano.properties.linkage import Molecules
from typing import List


def get_molecules(atoms: Atoms, use_cell_indices: bool=True):
    """
    Get the molecules from an atoms object.

    Args:
        atoms (Atoms): The atoms object.
        use_cell_indices (bool): Whether to connect molecules across periodic boundaries.

    Returns:
        list: A list of atoms objects, each corresponding to a molecule within Atoms.

    """
    molecules = []
    for mol in Molecules.get(atoms):
        molecules.append(mol.subset(atoms, use_cell_indices=use_cell_indices))
    return molecules


def reload_as_molecular_crystal(
                    images: List[Atoms],
                    parallel: bool=True,
                    cheap: bool=False,
                    wrap:bool=False):
    '''
    takes in a list of atoms,
    identifies molecules,
    unwraps molecules so that they're 'connected' 
    across periodic boundaries
    TODO: Test wrap in parallel mode and in cheap mode.

    Args:
        images (list): list of atoms objects
        parallel (bool): whether to use parallel processing. 
                        Only gives a speedup if there are 
                        ~100 structures per core. Fewer than that and we
                        just run in serial
                        Default: True
        cheap (bool): whether to use cheap mode. This is faster, but
                        doesn't work for all lists of structures. 
                        Avoid if possible.
                        It assumes that the molecules in each structure have the same indices.
                        (i.e. it only looks at the connectivity of the first structure and
                        assumes that the connectivity is the same for all structures).
                        Default: False.
        wrap (bool): whether to wrap molecules
                     Default: False

    Returns:
        images_new (list): list of atoms objects

    '''

    nimages = len(images)
    if cheap:
        # we assume the connectivity is the same for all images
        # definitely not a safe assumption, but it does speed things up a lot!
        mols = Molecules.get(images[0])
        images_new = []
        for atoms in images:
            atoms.wrap()
            temp = mols[0].subset(atoms, use_cell_indices=True)
            for mol in mols[1:]:
                temp.extend(mol.subset(atoms, use_cell_indices=True))
            if wrap:
                temp = wrap_molecule(temp)
            images_new.append(temp)
    else:
        if parallel:
            from multiprocessing import Pool, cpu_count
            # we can only run efficiently on about 100 structures per core
            # fewer than that and we'll just use serial
            ncores = min([cpu_count(), nimages//20])
        if parallel and ncores > 1:
            # pass wrap argument to unwrap_molecules
            # in a parallel pool
            with Pool(ncores) as p:
                jobs = [p.apply_async(unwrap_molecules, args=(atoms, wrap)) for atoms in images]
                images_new = [job.get() for job in jobs]
        
        else:
            images_new = [unwrap_molecules(image, wrap_each_molecule=wrap) for image in images]
    
    return images_new

def unwrap_molecules(atoms: Atoms, wrap_each_molecule: bool=False):
    atoms.wrap()
    mols = Molecules.get(atoms)
    temp = mols[0].subset(atoms, use_cell_indices=True)
    for mol in mols[1:]:
        mol_instance = mol.subset(atoms, use_cell_indices=True)
        if wrap_each_molecule:
            mol_instance = wrap_molecule(mol_instance)
        temp.extend(mol_instance)
    return temp

def wrap_molecule(atoms: Atoms,
                  max_dist: float=0.5):
    '''
    Translates molecule such that no bit of the 
    molecule sticks out passed max_dist in any direction.
    
    max_dist is in fractional units (default is 0.5).

    returns the atoms object with the molecule wrapped
    '''
    max_scaled_pos = np.max(atoms.get_scaled_positions(wrap=False), axis=0)
    mask = max_scaled_pos < max_dist
    T = atoms.cell.T.dot([1,1,1])
    T[mask] = 0
    atoms.translate(-T)
    return atoms


def molecule_collide(
                    atom: Atom,
                    atoms: Atoms,
                    vdw_scale:float = 1.0,
                    tolerance:float = 0.0):
    '''
    Checks if an atom is colliding with a group of atoms defined by the 
    initial_fractional_positions.


    We could naiively see if fractional position is within
    the bounding box of initial_fractional_positions, but that would
    not work in cases where the 'molecule' corresponding to initial_fractional_positions
    spans the cell boundary.
    '''
    # fractional_position = np.array(fractional_position)
    # cell = atoms.cell
    mols = Molecules.get(atoms, vdw_scale=vdw_scale)
    while len(mols) < 1:
        # we need to increase vdw_scale until we get a single molecule
        vdw_scale *= 1.1
        mols = Molecules.get(atoms, vdw_scale=vdw_scale)
    assert len(mols) == 1

    # now check that if we run the molecule finder with the same parameters, we get the same number of molecules
    atoms_plus = atoms.copy()
    atoms_plus.append(atom)
    mols_new = Molecules.get(atoms_plus, vdw_scale=vdw_scale)
    
    return len(mols_new) == len(mols)



    # # bounding box:
    # atoms_new = mols[0].subset(atoms, use_cell_indices=False)
    # box_positions = atoms_new.positions
    # min_pos = np.min(box_positions, axis=0)
    # max_pos = np.max(box_positions, axis=0)

    # # check the normalised fractional position
    # # against the bounding box
    # pos = cell.T.dot(fractional_position % 1)
    # if np.all(pos >= (min_pos-tolerance)) and np.all(pos <= (max_pos + tolerance)):
    #         return True
   
    # # if that wasn't a match, we still need to check its the periodic images
    # for i in range(-1,2):
    #     for j in range(-1,2):
    #         for k in range(-1,2):
    #             pos = cell.T.dot(fractional_position + np.array([i, j, k]))
    #             if np.all(pos >= (min_pos-tolerance)) and np.all(pos <= (max_pos + tolerance)):
    #                 return True
    # return False
    
