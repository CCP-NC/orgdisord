'''
Module containing class to help parsing cif file containing marked up disorder.
'''
from ase.io import read
from ase import Atoms
import warnings
from ase.spacegroup import get_spacegroup, crystal
import numpy as np
import logging


class CifParser:
    '''
    Class to parse cif file containing marked up disorder.
    '''
    def __init__(self, cif_file):
        self.cif_file = cif_file
        # use ASE to read in cif and info
        self.atoms = read(self.cif_file, store_tags=True)
        self.cell = self.atoms.cell
        # save the info dictionary
        self.info = self.atoms.info
        self.logger = logging.getLogger("sodorg.parse_cif_file")
        self.logger.debug("\n\n")
        self.logger.debug("------------------------------")
        self.logger.debug("---    PARSING CIF FILE    ---")
        self.logger.debug("------------------------------")
        self.logger.debug(f"{self.cif_file}")
        

        # Disorder assemblies and groups
        # Make sure the info dictionary contains the disorder tags
        if not '_atom_site_disorder_group' in self.info:
            raise ValueError(f'Cif file {self.cif_file} does not contain disorder group tags. Please edit the cif file appropriately.')
            # TODO more helpful error message -- guidance on using these tags
        if not '_atom_site_disorder_assembly' in self.info:
            warnings.warn(f"Cif file {self.cif_file} does not contain disorder assembly tags. Best practice is to use these to disambiguate disorder assemblies/groups.")
            assemblies_present = False
        else:
            assemblies_present = True

        # split into assembly and disorder groups
        self.disorder_groups = self.split_disorder_groups(assemblies_present=assemblies_present)
        # structure is {assembly label: {disorder group label: [indices]}}
        self.nassemblies = len(self.disorder_groups.keys())

        if self.nassemblies == 0:
            # this should never happen now that we check above... 
            raise ValueError(f'Cif file {self.cif_file} does not contain any marked up disorder assemblies/groups')

        # for now raise error if there are different numer of disorder groups in different assemblies
        ngroups = [len(group) for asslabel, group in self.disorder_groups.items()]
        if len(set(ngroups)) > 1:
            self.logger.warn(f'Cif file {self.cif_file} contains different numbers'
            ' of disorder groups in different assemblies. We therefore assume that they are not correlated.\n'
            'This is a bit experimental, so double check your results!')        


        # spacegroup info:
        sg = self.atoms.info['spacegroup']
        if sg.centrosymmetric:
            warnings.warn('Warning: crystal read as centrosymmetric, but we will proceed as if it is not!'
                        'This is usually fine in my experience... ')
            sg._centrosymmetric = 0
        self.sg = sg
        self.nops = sg.nsymop
        # symmetry operations — tuples of (rotation, translation)
        self.symops = sg.get_symop()

        self.logger.debug(f'Found {self.nops} symmetry operations. Spacegroup {sg}')


        # TODO: better way to get Z?
        try:
            Z = self.info['_cell_formula_units_z']
        except:
            Z = 4
            self.logger.warn(f"WARNING: Couldn't parse Z -- taking Z = {Z} for now and proceeding")
        self.Z = Z
        self.Zprime = Z / self.nops

        # scaled coordinates of the asymmetric sites:
        asymmetric_x = self.info['_atom_site_fract_x']
        asymmetric_y = self.info['_atom_site_fract_y']
        asymmetric_z = self.info['_atom_site_fract_z']
        self.asymmetric_scaled_coords = np.array([asymmetric_x, asymmetric_y, asymmetric_z]).T
        # element symbols of the asymmetric sites:
        self.asymmetric_symbols = np.array(self.info['_atom_site_type_symbol'])
        # labels of the asymmetric sites:
        self.asymmetric_labels = np.array(self.info['_atom_site_label'])
        # occupancies of the asymmetric sites:
        occupancies = np.array(self.info['_atom_site_occupancy'])
        self.occupancies = occupancies
        # let's define the ordered sites as those with 100% occupancy
        self.ordered_mask = occupancies == 1

        # report the disorder groups to user
        for asslabel, groups in self.disorder_groups.items():
            self.logger.debug(f'Assembly {asslabel} contains {len(groups)} disorder groups:')
            for grouplabel, group in groups.items():
                self.logger.debug(f'  Group {grouplabel} contains {len(group)} sites:')
                for ig in group:
                    self.logger.debug(f"     label: {self.info['_atom_site_label'][ig]: >8}, index: {ig: 5d},  species: {self.asymmetric_symbols[ig]: >4}, occupancy: {occupancies[ig]: >5.2f}")
                tol = 1e-4
                if any(abs(1 - occupancies[group]) < tol):
                    self.logger.warning(f"WARNING: Disorder group {grouplabel} from assembly {asslabel} contains site(s) with full occupancy. Please check the log file and your cif file carefully!")
        self.logger.debug("Check to make sure these are the groups you were expecting!")
        self.logger.debug("Otherwise, try editing the cif file and rerunning.")

        # build out the full crystal with the ordered sites
        self.gen_ordered_atoms() # (saves to self.ordered_atoms)


    def gen_ordered_atoms(self, symprec=1e-4):
        '''
        if there are any fully ordered sites, 
        build out the full crystal with the ordered sites
        otherwise return empty atoms object
        '''
        mask = self.ordered_mask
        symbols = self.asymmetric_symbols[mask]
        labels = self.asymmetric_labels[mask]
        scaled_coords = self.asymmetric_scaled_coords[mask]
        if len(symbols) > 0:
            ordered_atoms = crystal(symbols=symbols,
                                    basis=scaled_coords,
                                    occupancies=None,
                                    spacegroup=self.atoms.info['spacegroup'],
                                    cell=self.atoms.cell,
                                    onduplicates='warn',
                                    symprec=symprec,
                                    pbc=True,
                                    primitive_cell=False)
            # now map the labels onto the corresponding atoms
            kinds = ordered_atoms.get_array('spacegroup_kinds')
            all_labels = np.array([labels[kind] for kind in kinds])
            ordered_atoms.set_array('labels', all_labels)
            
            self.ordered_atoms = ordered_atoms
        
        # temp_atoms = Atoms(symbols=symbols, scaled_positions = scaled_coords, cell = self.cell)
        # if len(temp_atoms) > 0:
        #     ordered_atoms = crystal(temp_atoms, spacegroup=sg.no,
        #                             setting=sg.setting, primitive_cell=False, occupancies=None )
        #     self.ordered_atoms = ordered_atoms
            
        else:
            # empty box of the right shape
            self.ordered_atoms = Atoms(cell=self.atoms.cell, pbc=True)

    def split_assembly_groups(self, atoms):
        '''
        
        Note: atoms.info must contain _atom_site_disorder_assembly
        Returns indices grouped by the sorted assembly label


        From the IUCR definition of _atoms_site_disorder_assembly:    
        > A code which identifies a cluster of atoms that show long-range
            positional disorder but are locally ordered. Within each such
            cluster of atoms, _atom_site_disorder_group is used to identify
            the sites that are simultaneously occupied. This field is only
            needed if there is more than one cluster of disordered atoms
            showing independent local order.
                    
        '''
        
        disorder_assemblies = atoms.info['_atom_site_disorder_assembly']
        labels = list(set(disorder_assemblies))
        if '.' in labels:
            labels.remove('.')
        labels.sort()
        self.logger.debug(f'Found the following disorder assembly labels: {labels}')
        
        indices = range(len(disorder_assemblies))
        # dictionary of indices grouped by the sorted assembly label
        indices_by_assembly = {label: [] for label in labels}
        for i, label in enumerate(disorder_assemblies):
            if label in labels:
                indices_by_assembly[label].append(indices[i])
        return indices_by_assembly
        
        # assembly_groups = [[idx for idx in indices if disorder_assemblies[idx] == label] for label in labels]

        # return assembly_groups


    def split_disorder_groups(self, assemblies_present=True):
        '''
        Note: atoms.info must contain _atom_site_disorder_assembly and _atom_site_disorder_group
        Returns indices of asymmetric atoms grouped by sorted disorder group within each assembly group

        From the IUCR definition of _atoms_site_disorder_group:
        >  A code which identifies a group of positionally disordered atom
            sites that are locally simultaneously occupied. Atoms that are
            positionally disordered over two or more sites (e.g. the hydrogen
            atoms of a methyl group that exists in two orientations) can
            be assigned to two or more groups. Sites belonging to the same
            group are simultaneously occupied, but those belonging to
            different groups are not. A minus prefix (e.g. '-1') is used to
            indicate sites disordered about a special position.

            *** This data item would not in general be used in a
            macromolecular data block. ***
        '''
        site_disorder_groups = self.info['_atom_site_disorder_group']
        labels = list(set(site_disorder_groups))
        if '.' in labels:
            labels.remove('.')
        labels.sort()

        self.logger.debug(f'Found the following disorder group labels: {labels}')
        
        
        if assemblies_present:
            assembly_groups = self.split_assembly_groups(self)
            self.logger.debug(f'Found {len(assembly_groups.keys())} assembly group(s):')
            for label, indices in assembly_groups.items():
                self.logger.debug(f'  Assembly {label}: {len(indices)} sites: {indices}')
        else:
            # no assemblies present, so just group by disorder group and give label 'A'
            assembly_groups = {'A': [idx for idx in range(len(site_disorder_groups))]}
        
        disorder_groups = {asslabel: {label: [] for label in labels} for asslabel in assembly_groups.keys()}
        
        for asslabel, indices in assembly_groups.items():
            for i, label in enumerate(site_disorder_groups):
                if label in labels:
                    if len(indices) > 0:
                        if i in indices:
                            disorder_groups[asslabel][label].append(i)
                    else:
                        del disorder_groups[asslabel][label]
        

            self.logger.debug(f'Found {len(disorder_groups[asslabel])} disorder groups in assembly group {asslabel}:')
            self.logger.debug(f'    Python indices: {disorder_groups[asslabel]}')
        
        # warning if any are in an assebmly but not in a disorder group
        for asslabel, indices in assembly_groups.items():
            for i in indices:
                if site_disorder_groups[i] not in labels:
                    self.logger.warning(f'Warning: Site {i} in assembly {asslabel} has no disorder group. Double check the CIF file.')

        return disorder_groups

