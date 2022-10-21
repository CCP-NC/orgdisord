'''
Module containing class to help parsing cif file containing marked up disorder.

'''
from turtle import position
from ase.io import read
from ase import Atoms, Atom
import warnings
from ase.spacegroup import get_spacegroup, crystal
from sodorg_renewal.utils import molecule_collide
from sodorg_renewal.disordered_structure import DisorderedStructure, DisorderGroup, DisorderAssembly

import numpy as np
import logging


class CifParser:
    '''
    Class to parse cif file containing marked up disorder.
    '''
    def __init__(self, cif_file, correlated_assemblies=False, molecular_crystal=True, symprec = 1e-4):
        self.cif_file = cif_file
        self.symprec = symprec
        # use ASE to read in cif and info
        self.atoms = read(self.cif_file, store_tags=True)
        self.cell = self.atoms.cell
        self.correlated_assemblies = correlated_assemblies
        self.molecular_crystal = molecular_crystal
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
        if '_atom_site_occupancy' in self.info:
            occupancies = np.array(self.info['_atom_site_occupancy'])
        else:
            occupancies = np.ones(len(self.asymmetric_symbols))
            self.logger.warn("Warning: No occupancies found in CIF file. Assuming all occupancies are 1.")

        self.occupancies = occupancies
        # let's define the ordered sites as with disorder group '.'
        # make sure we're comparing strings
        self.ordered_mask = np.array(self.info['_atom_site_disorder_group']).astype(str) == '.'
        # self.ordered_mask = occupancies == 1

        # split into assembly and disorder groups
        self.disorder_groups = self.split_disorder_groups(assemblies_present=assemblies_present)
        # structure is {assembly label: {disorder group label: [indices]}}
        self.nassemblies = len(self.disorder_groups.keys())

        all_kinds = self.atoms.get_array('spacegroup_kinds')
        all_disorder_tags = np.zeros(len(all_kinds))
        iass = 0
        for asslabel, groups in self.disorder_groups.items():
            iass += 1e3
            for grouplabel, group in groups.items():
                for idx in group:
                    all_disorder_tags[np.where(all_kinds == idx)[0]] = iass + grouplabel
        
        self.atoms.set_tags(all_disorder_tags)




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
            self.logger.debug('The crystal read as centrosymmetric, but we will proceed as if it is not.'
                        'This is usually fine in my experience... ')
            sg._centrosymmetric = 0
        self.sg = sg
        self.nops = sg.nsymop
        # symmetry operations — tuples of (rotation, translation)
        self.symops = sg.get_symop()

        self.logger.debug(f'Found {self.nops} symmetry operations. Spacegroup {sg}')


        # TODO: is there a better way to get Z?
        ## warn if _cell_formula_units_Z is not present in info dict
        if not '_cell_formula_units_z' in self.info:
            self.logger.warn(
                f'Cif file {self.cif_file} does not contain _cell_formula_units_z tag. '
                'Assuming Z = number of symmetry operations. '
                'This is not correct in many cases - please edit the cif file appropriately.'
                )
            Z = self.nops
        else:
            Z = self.info['_cell_formula_units_z']
        
        self.Z = Z
        self.Zprime = Z / self.nops

        self.logger.debug(f'Z = {Z}, Zprime = {self.Zprime}')



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
        self.ordered_atoms = self.gen_ordered_atoms()

        # ds = self.get_disordered_structure()
        # self.logger.debug(f'{ds}')





    def gen_ordered_atoms(self):
        '''
        if there are any fully ordered sites, 
        build out the full crystal with the ordered sites
        otherwise return empty atoms object
        '''
        mask = self.ordered_mask

        if any(mask):
            # -> we have some ordered sites
            self.logger.debug("Found ordered sites. Building full crystal with ordered sites.")
            # any ordered sites with occupancy different to 1? 
            partially_occupied = [idx for idx, occ in enumerate(self.occupancies) if occ != 1 and mask[idx]]
            if len(partially_occupied) > 0:
                self.logger.warning(f'Warning: Sites {self.asymmetric_labels[partially_occupied]} '
            'are in the ordered sites group but have occupancy != 1.\n'
            'They will be included in the generated structures but will probably be wrong...\n'
            'Please check your CIF file carefully!')

            symbols = self.asymmetric_symbols[mask]
            labels = self.asymmetric_labels[mask]
            scaled_coords = self.asymmetric_scaled_coords[mask]
            ordered_atoms = crystal(symbols=symbols,
                                    basis=scaled_coords,
                                    occupancies=None,
                                    spacegroup=self.atoms.info['spacegroup'],
                                    cell=self.atoms.cell,
                                    onduplicates='warn',
                                    symprec=self.symprec,
                                    pbc=True,
                                    primitive_cell=False)
            # now map the labels onto the corresponding atoms
            kinds = ordered_atoms.get_array('spacegroup_kinds')
            all_labels = np.array([labels[kind] for kind in kinds])
            ordered_atoms.set_array('labels', all_labels)
        else:
            # empty box of the right shape
            ordered_atoms = Atoms(cell=self.atoms.cell, pbc=True)
        return ordered_atoms

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
        
        # warn if any sites are neither in any group nor in the ordered sites
        unassigned_sites =  [idx for idx in range(len(site_disorder_groups)) 
                                 if site_disorder_groups[idx] == '.' and not self.ordered_mask[idx]]
        if len(unassigned_sites) > 0:
            self.logger.warning(f'Warning: Sites {self.asymmetric_labels[unassigned_sites]} '
            'are neither in any disorder group nor in the ordered sites! '
            'These will *not* be included in the generated structures. '
            'Please check your CIF file carefully!')
        
        if '.' in np.array(site_disorder_groups, dtype = 'U')[self.occupancies != 1]:
            self.logger.warn(f'Warning: Some sites are not in a disorder group but have occupancy != 1.' )
        
        
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
                    if i in indices:
                        disorder_groups[asslabel][label].append(i)
        
            self.logger.debug(f'Found {len(disorder_groups[asslabel])} disorder groups in assembly group {asslabel}:')
            self.logger.debug(f'    Python indices: {disorder_groups[asslabel]}')
        
        # delete any empty disorder assemblies or disorder groups
        disorder_groups = {
            ass_k: {
                group_k: group_v
                for group_k, group_v in ass_v.items() if group_v
            } 
            for ass_k, ass_v in disorder_groups.items() if ass_v
        }
        
        # warning if any are in an assembly but not in a disorder group
        for asslabel, indices in assembly_groups.items():
            for i in indices:
                if site_disorder_groups[i] not in labels:
                    self.logger.warning(f'Warning: Site {i} ({self.asymmetric_labels[i]}) in assembly {asslabel} has no disorder group. '
                    'This site will be considered *ordered*. \n'
                    'Double check the CIF file.')

        return disorder_groups


    # def is_assembly_asymmetric(self, assembly_label):
    #     '''
    #     Assembly is considered asymmetric if any of the disorder groups 
    #     have a negative label.
    #     This isn't used anywhere -- I'm not sure if it's the correct approach...
        
    #     '''
    #     # make sure we've already got the disorder groups
    #     if not hasattr(self, 'disorder_groups'):
    #         self.disorder_groups = self.split_disorder_groups(assemblies_present=True)
        
    #     disorder_group_labels = sorted(self.disorder_groups[assembly_label].keys())
    #     asymmetric_disorder_groups = [l for l in disorder_group_labels if l < 0]
        
    #     # check any of the disorder group labels is negative
    #     asymmetric_assembly = False
    #     if len(disorder_group_labels) == 1:
    #         if disorder_group_labels[0] < 0:
    #             # we have a single asymmetric disorder group
    #             asymmetric_assembly = True
    #         else:
    #             self.logger.warn(f'Warning: assembly {assembly_label} has just a single non-asymmetric disorder group.')
    #     else:
    #         if len(asymmetric_disorder_groups) > 0:
    #             # check if all the disorder group labels are negative
    #             if len(asymmetric_disorder_groups) == len(disorder_group_labels):
    #                 asymmetric_assembly = True
    #             else:
    #                 self.logger.warn(f'Warning: assembly {assembly_label} has a mix of positive and negative disorder group labels -- not sure what to do here!' )
    #     return asymmetric_assembly



    def get_disordered_structure(self):
        '''
        Retruns a DisorderedStructure object.
        '''
        # get ordered part
        ordered_atoms = self.ordered_atoms

        # get disordered part
        disorder_assemblies = []

        
        # loop over assemblies
        for iass, (asslabel, groups) in enumerate(self.disorder_groups.items()):
            asstag = (iass+1)*1e3 if self.nassemblies > 1 else 0

            # loop over groups
            disorder_groups = []
            for igroup, (grouplabel, group) in enumerate(groups.items()):
                grouptag = igroup + 1 + asstag
                # create ASE Atoms object for group
                symbols = self.asymmetric_symbols[group]
                labels = self.asymmetric_labels[group]
                scaled_coords = self.asymmetric_scaled_coords[group]
                occupancies = self.occupancies[group]
                atoms = Atoms(symbols = symbols,
                              scaled_positions = scaled_coords,
                              cell=self.cell, pbc=True)
                atoms.set_array('labels', labels)
                atoms.set_array('occupancies', occupancies)

                group_occupancy = np.mean(occupancies)
                
                # create DisorderedGroup object
                disorder_group = DisorderGroup(
                            label=str(grouplabel),
                            atoms =  atoms,
                            symmetry_operations=self.symops,
                            occupancy=group_occupancy,
                            tag=grouptag)

                group_symops = disorder_group.get_group_symmetry_operations(self.sg, self.Z)
                # set the symmetry operations for the group
                disorder_group.symmetry_operations = group_symops
                # add to list of groups
                disorder_groups.append(disorder_group)
            # create DisorderedAssembly object out of list of DisorderedGroup objects
            disorder_assembly = DisorderAssembly(label = asslabel,
                                                disorder_groups = disorder_groups,
                                                tag = asstag)
            # add to list of assemblies
            disorder_assemblies.append(disorder_assembly)

        # create DisorderedStructure object to be returned
        disordered_structure = DisorderedStructure(
                            ordered_atoms = ordered_atoms,
                            Z = self.Z,
                            spacegroup = self.sg,
                            disorder_assemblies = disorder_assemblies,
                            molecular_crystal=self.molecular_crystal,
                            correlated_assemblies = self.correlated_assemblies,
        )
        return disordered_structure