"""Console script for sodorg_renewal."""
import sys
import os
import click
from sodorg_renewal.parse_cif_file import CifParser
from sodorg_renewal.enumerate import OrderedfromDisordered
from sodorg_renewal.merge import merge_structures
import spglib
import time
from ase.io import write
from ase.units import _amu
import pandas as pd
import numpy as np
import time
import logging
logging.captureWarnings(True)

# TODO: sort out DF order/ what info to include
# TODO sort out DF printing/saving to file

# new log header
LOGHEADER = '''


###############################################################################
#                          -----  SODORG-PY  -----                            #
#                                                          J. Kane Shenton    #
#                                                                   CCP-NC    #
###############################################################################
'''
@click.command()
@click.pass_context
@click.argument('cif_file', type=click.Path(exists=True))
# add option to specify supercell
@click.option('--supercell', '-s', nargs=3, type=click.INT, default=[1,1,1], help='Supercell size e.g. 2 1 1')
# add option to specify max number of iterations
@click.option('--maxiters', '-N', type=click.INT, default=512, help='Maximum number of configurations generated')
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
# option to suppress writing out structures
@click.option('--no_write', is_flag=True, default=False, help='Suppress writing out structures?')
# add option to output to specified file
@click.option('--prefix', type=click.STRING, default='sodorg', help='Prefix for output file names.')
# option format options
@click.option('--format', '-f', type=click.Choice(['cif', 'xyz', 'poscar', 'cell']), default='xyz', help='Format of output strucuture files.')

# add option to specify verbosity
@click.option('--quiet', '-q', is_flag=True, default=False, help='Suppress output?')
# add option to specify log file
@click.option('--log', '-l', type=click.Path(exists=False), default="sodorg.log", help='Log file name.')

# add option to view structures using ASE gui
@click.option('--view', is_flag=True, default=False, help='View structures using ASE gui?')
# Add option to specify oxidation states
@click.option('--ox', type=click.Tuple([str, int]), multiple=True, default=None, help='Specify oxidation states for each species e.g. --ox H 1 --ox O -2')


# ox = {'C':1, 'O':-2}
# # combine all ox options into a dictionary
# def get_oxidations(ox):
#     ox_dict = {}
#     if ox:
#         for ox_name, ox_val in ox:
#             ox_dict[ox_name] = ox_val
#     return ox_dict


def main(ctx,
         cif_file, 
         supercell, 
         maxiters, 
         exclude_ordered, 
         random,
         merge,
         symprec, 
         quiet,
         algo,
         use_disordered_only,
         no_write,
         prefix,
         format,
         ox,
         log,
         view
         ):

    """Command line interface for sodorg_renewal.
    """
    # set up logging
    logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s \t %(message)s')
    # add console handler
    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    consolefmt = logging.Formatter('%(message)s')
    console.setFormatter(consolefmt)
    logger = logging.getLogger("sodorg")
    # add the handler to the main logger
    logger.addHandler(console)
    logger.info(LOGHEADER)
    
    # log the Click command line arguments and options used
    logger.debug("Run with these command line arguments: ")
    for param, value in ctx.params.items():
        logger.debug(f"          {param}: {value}")

    df = pd.DataFrame()

    logger.info(f'Parsing disorder in cif file: {cif_file}')
    # parse cif file containing disordered structure
    cif = CifParser(cif_file)
    symops = cif.symops

    # enumerate ordered configurations
    od = OrderedfromDisordered(cif, quiet=quiet)
    images, configs = od.get_supercell_configs(
                        supercell = supercell,
                        maxiters = maxiters,
                        exclude_ordered = exclude_ordered,
                        random_configs=random,
                        return_configs=True,
                        )

    # -- DATA FRAME -- #
    # store information to dataframe
    df['Unique index'] = np.arange(len(images))
    df['Configuration'] = configs
    ratios = [1-sum(config) / len(config) for config in configs]
    # TODO should we define this as the reverse? 
    df['Ratio of maj:min'] = ratios
    
    df['Supercell'] = "x".join(map(str, supercell))
    # TODO: units
    df['Free Energy per atom / kJ/mol'] = np.NaN
    df['Enthalpy per atom / kJ/mol'] = np.NaN
    df['Entropy per atom / kJ/mol'] = np.NaN
    # volume of simulation cell
    df['Volume/A^3'] = [image.get_volume() for image in images]
    # number of atoms
    df['N atoms'] = [len(image) for image in images]
    # density
    masses = [sum(image.get_masses()) for image in images]
    density_units = (_amu / 1e-30) * 1e-3 # > g/cm^3
    density = [(mass / volume) * density_units for mass, volume in zip(masses, df['Volume/A^3'])]
    df['Density/g/mol^3'] = density
    #/-- DATA FRAME -- #
    
    
    # merge symmetrically equivalent configurations?
    if merge:
        oxidation_states = {}
        if algo == 'ewald':
            if ox:
                for ox_name, ox_val in ox:
                    oxidation_states[ox_name] = ox_val
            logger.info(f"Oxidation states used: {oxidation_states}")

        logger.info('Merging structures...')
        start = time.time()
        groups, group_indices = merge_structures(
                                images,
                                algo=algo,
                                symops=symops,
                                use_disordered_only = use_disordered_only,
                                symprec=symprec,
                                oxidation_states=oxidation_states,
                                return_group_indices=True,
                                quiet = quiet)
        
        stop = time.time()
        duration = stop-start
        logger.info(f'Merging took {duration:8.2f} s and found {len(groups)} groups')
        unique_spacegroups = [spglib.get_spacegroup(g[0]) for g in groups]
        for sg, g in zip(unique_spacegroups, groups):
            logger.info(f"Spacegroup: {sg:<20} multiplicity: {g[1]}")
        
        # overwrite images with merged images
        images = [g[0] for g in groups]

        # tag images with group indices in df
        for ig, group_inds in enumerate(group_indices):
            for g in group_inds:
                df.loc[g, 'Unique index'] = ig
        
        # Update dataframe with merged information
        # merge df by Unique index
        df = df.groupby('Unique index').aggregate(np.mean).reset_index()
        # add multiplicity to df
        df['multiplicity'] = [g[1] for g in groups]
        # add spacegroup to df
        df['spacegroup'] = unique_spacegroups
        # add groups to df
        df['Group indices'] = group_indices
    
    # log the df
    logging.info('\n\t'+ df.to_string().replace('\n', '\n\t')) 
    # save the df to csv
    df.to_csv(prefix+'.csv', index=False)
    if not no_write:
        # write structures to file(s)
        directory = f"{prefix}-results"
        logger.info(f'Writing structures to file(s) in directory: {directory}')
        # check if directory already exists -- if so rename it (backup)
        if os.path.isdir(directory):
            logger.warning(f'Directory {directory} already exists. Backing up existing directory.')
            # backup existing directory
            backup_dir = f"{directory}-backup-{time.strftime('%Y%m%d-%H%M%S')}"
            if os.path.isdir(backup_dir):
                raise ValueError(f'Backup directory {backup_dir} already exists. Please remove it.')
            os.rename(directory, backup_dir)

        os.mkdir(directory)
        nleading_zeros = len(str(len(images)))
        for i, image in enumerate(images):
            filename = f"{prefix}-{str(i).zfill(nleading_zeros)}.{format}"
            image.write(f"{directory}/{filename}")
            logger.info(f'Wrote structure to file: {directory}/{filename}')

    # write out pandas dataframe containing the groups and multiplicities


    if view:
        # view images using ASE gui
        from ase.visualize import view
        view(images)


    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
