"""Console script for sodorg_renewal."""
import sys
import click
from parse_cif_file import CifParser
from enumerate import OrderedfromDisordered
from merge import merge_structures
import spglib
import time
from ase.io import write


@click.command()
@click.argument('cif_file', type=click.Path(exists=True))
# add option to specify supercell
@click.option('--supercell', '-s', nargs=3, type=click.INT, default=[1,1,1], help='Supercell size e.g. 2 1 1')
# add option to specify max number of iterations
@click.option('--maxiters', '-N', type=click.INT, default=100, help='Maximum number of configurations generated')
# add option to include ordered configurations
@click.option('--exclude_ordered', is_flag=True, default=False, help='Exclude ordered part of structure in output structures?')
# add option to generate random configurations
@click.option('--random', '-r', is_flag=True, default=False, help='Generate random configurations?')
# add option to merge structures
@click.option('--merge', '-m', is_flag=True, default=False, help='Merge equivalent structures?')
#add option to specify algorithm
@click.option('--algo', '-a', type=click.Choice(['symm', 'rematch', 'ewald']), default='symm', help='Algorithm used for merging equivalent structures.')
# add option to specify use of disordered only in merging
@click.option('--use_disordered_only', '-d', is_flag=True, default=True, help='Use only the disordered part of structure in merging?')
# add option to specify symmetry precision
@click.option('--symprec', '-p', type=click.FLOAT, default=1e-4, help='Symmetry precision used for merging equivalent structures.')
# add option to output to specified file
@click.option('--output', '-o', type=click.Path(exists=False), default=None, help='Output file name. Must be a format able to hold multiple structures e.g. .cif, .traj, .xyz and written by ASE')
# add option to specify verbosity
@click.option('--verbose', '-v', is_flag=True, default=False, help='Print verbose output?')
# add option to view structures using ASE gui
@click.option('--view', is_flag=True, default=False, help='View structures using ASE gui?')
# Add option to specify oxidation states
@click.option('--ox', type=click.Tuple([str, int]), multiple=True, default=None)


# ox = {'C':1, 'O':-2}
# # combine all ox options into a dictionary
# def get_oxidations(ox):
#     ox_dict = {}
#     if ox:
#         for ox_name, ox_val in ox:
#             ox_dict[ox_name] = ox_val
#     return ox_dict


def main(cif_file, 
         supercell, 
         maxiters, 
         exclude_ordered, 
         random,
         merge,
         symprec, 
         verbose,
         algo,
         use_disordered_only,
         output,
         ox,
         view
         ):

    """Command line interface for sodorg_renewal.
    """

    print('Parsing disorder in cif file...')
    # parse cif file containing disordered structure
    cif = CifParser(cif_file)
    symops = cif.symops

    # enumerate ordered configurations
    od = OrderedfromDisordered(cif, verbose=verbose)
    images = od.get_supercell_configs(
                    supercell = supercell,
                    maxiters = maxiters,
                    exclude_ordered = exclude_ordered,
                    random_configs=random)

    
    
    # merge symmetrically equivalent configurations?
    if merge:
        oxidation_states = {}
        if algo == 'ewald':
            if ox:
                for ox_name, ox_val in ox:
                    oxidation_states[ox_name] = ox_val
        print(oxidation_states)

        print('Merging structures...')
        start = time.time()
        groups = merge_structures(
                            images,
                            algo=algo,
                            symops=symops,
                            use_disordered_only = use_disordered_only,
                            symprec=symprec,
                            verbose=verbose,
                            oxidation_states=oxidation_states)
        
        stop = time.time()
        duration = stop-start
        print(f'Merging took {duration:8.2f} s and found {len(groups)} groups')
        if verbose:
            for g in groups:
                print(f"Spacegroup: {spglib.get_spacegroup(g[0])}\t multiplicity: {g[1]}")
        
        # overwrite images with merged images
        images = [g[0] for g in groups]

    if output:
        # write images to file

        # for i, atoms in enumerate(images):
        #     write(f'{output}/image-{i:04d}.xyz', atoms)
        write(output, images)

    if view:
        # view images using ASE gui
        from ase.visualize import view
        view(images)


    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
