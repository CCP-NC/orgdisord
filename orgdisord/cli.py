"""Console script for orgdisord."""
import sys
import os
import click
from orgdisord.disordered_structure import from_disorder_components
from orgdisord.parse_cif_file import CifParser
from orgdisord.enumerate import OrderedfromDisordered
from orgdisord.merge import merge_structures
from orgdisord.utils import get_new_labels
import spglib
import time
from ase.io import read
from ase.units import _amu
from ase.units import kB, kJ, mol  # kB is the Boltzmann constant in eV/K
from ase.utils import atoms_to_spglib_cell
import pandas as pd
import numpy as np
import logging

logging.captureWarnings(True)

# TODO: sort out DF order/ what info to include
# TODO sort out DF printing/saving to file

# new log header
LOGHEADER = """


###############################################################################
#                          -----  orgdisord  -----                            #
#                                                          J. Kane Shenton    #
#                                                                   CCP-NC    #
###############################################################################
"""


@click.group()
def cli():
    pass


@cli.command("enumerate", context_settings={"show_default": True})
@click.pass_context
@click.argument("files", nargs=-1, type=click.Path(exists=True))
# add option to specify supercell
@click.option(
    "--supercell",
    "-s",
    nargs=3,
    type=click.INT,
    default=[1, 1, 1],
    help="Supercell size e.g. 2 1 1",
)
# add option to specify max number of iterations
@click.option(
    "--maxiters",
    "-N",
    type=click.INT,
    default=512,
    help="Maximum number of configurations generated",
)
# add option to include ordered configurations
@click.option(
    "--exclude_ordered",
    is_flag=True,
    default=False,
    help="Exclude ordered part of structure in output structures?",
)
# add option to generate random configurations
@click.option(
    "--random",
    "-r",
    is_flag=True,
    default=False,
    help="Generate random configurations?",
)
# add option to merge structures
@click.option(
    "--merge", "-m", is_flag=True, default=False, help="Merge equivalent structures?"
)
# add option to specify algorithm
@click.option(
    "--algo",
    "-a",
    type=click.Choice(["symm", "rematch", "ewald"]),
    default="symm",
    help="Algorithm used for merging equivalent structures.",
)
# add option to specify use of disordered only in merging
@click.option(
    "--use_disordered_only",
    "-d",
    is_flag=True,
    default=True,
    help="Use only the disordered part of structure in merging?",
)
# add option to specify symmetry precision
@click.option(
    "--symprec",
    "-p",
    type=click.FLOAT,
    default=1e-4,
    help="Symmetry precision used for merging equivalent structures.",
)
# option to suppress writing out structures
@click.option(
    "--no_write", is_flag=True, default=False, help="Suppress writing out structures?"
)
# add option to output to specified file
@click.option(
    "--prefix",
    type=click.STRING,
    default="orgdisord",
    help="Prefix for output file names.",
)
# option format options
@click.option(
    "--format",
    "-f",
    type=click.Choice(["cif", "xyz", "poscar", "cell"]),
    default="xyz",
    help="Format of output strucuture files.",
)

# add option to specify verbosity
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress output?")
# add option to specify log file
@click.option(
    "--log",
    "-l",
    type=click.Path(exists=False),
    default="orgdisord.log",
    help="Log file name.",
)

# add option to view structures using ASE gui
@click.option(
    "--view", is_flag=True, default=False, help="View structures using ASE gui?"
)
# Add option to specify oxidation states
@click.option(
    "--ox",
    type=click.Tuple([str, int]),
    multiple=True,
    default=None,
    help="Specify oxidation states for each species e.g. --ox H 1 --ox O -2",
)
# option to not reload as molecular crystal
@click.option(
    "--not_molecular_crystal",
    is_flag=True,
    default=False,
    help="Do not reload as molecular crystal?",
)
# correlated assemblies?
@click.option(
    "--correlated_assemblies",
    "-c",
    is_flag=True,
    default=False,
    help="Are the disorder assemblies correlated?"
    "i.e. are the disorder groups within each assembly to be selected with the same index? "
    "(The number of disorder groups in each assembly must be the same in such cases.)",
)
# fix ratio flag
@click.option(
    "--fix-ratio",
    is_flag=True,
    default=False,
    help="Fix the ratio of the disorder groups.",
)
# manually specify ratios
@click.option(
    "--ratio",
    nargs=1,
    type=click.FLOAT,
    default=None,
    help="Manually specify the ratio of the disorder groups. "
    "Must be used with --fix-ratio."
    "e.g. --ratio 0.75 would constrain the ratios of disorder group 1 and 2 to be 3:1"
    "NOTE: This is only implemented for 1-assembly, 2-group structures.",
)
# ratio tolerance
@click.option(
    "--ratio-tol",
    type=click.FLOAT,
    default=0.01,
    help="Tolerance for how close the ratio of disorder groups needs to be. "
    "See atol in numpy.isclose."
    "Only used if --fix-ratio is specified.",
)
# check species
@click.option(
    "--ignore-species",
    is_flag=True,
    default=False,
    help="Ignore species when merging structures?",
)
def main(
    ctx,
    files,
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
    view,
    not_molecular_crystal,
    correlated_assemblies,
    fix_ratio,
    ratio,
    ratio_tol,
    ignore_species,
):
    """Command line interface for orgdisord."""
    # set up logging
    logging.basicConfig(
        filename=log, level=logging.DEBUG, format="%(asctime)s \t %(message)s"
    )
    # add console handler
    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    consolefmt = logging.Formatter("%(message)s")
    console.setFormatter(consolefmt)
    logger = logging.getLogger("orgdisord")
    # add the handler to the main logger
    logger.addHandler(console)
    logger.info(LOGHEADER)

    # log the Click command line arguments and options used
    logger.debug("Run with these command line arguments: ")
    for param, value in ctx.params.items():
        logger.debug(f"          {param}: {value}")

    df = pd.DataFrame()
    nfiles = len(files)
    if nfiles == 1:
        logger.info(f"Parsing disorder in cif file: {files[0]}")
        # parse cif file containing disordered structure
        cif = CifParser(
            files[0],
            correlated_assemblies=correlated_assemblies,
            molecular_crystal=not not_molecular_crystal,
        )

        # get the disordered structure
        disordered_structure = cif.get_disordered_structure()

    elif nfiles == 2:
        logger.info(f"Parsing disorder in files: {files[0]} and {files[1]}")
        logger.info(
            "Assuming the two files contain two fully ordered structures, one for each disorder group."
        )

        disorder_components = []
        for f in files:
            labels = None
            if f.endswith(".cif"):
                # if they are CIF files, parse the extra tags
                atoms = read(f, store_tags=True)
            else:
                atoms = read(f)  ## hope ASE can read it

            # -- labels -- #
            if not atoms.has("labels"):
                if "_atom_site_label" in atoms.info:
                    labels = atoms.info["_atom_site_label"]
                elif "labels" in atoms.info:
                    labels = atoms.info["labels"]
                else:
                    labels = get_new_labels(atoms)
                atoms.set_array("labels", np.array(labels))

            # -- occupanices -- #
            if not atoms.has("occupancies"):
                if "_atom_site_occupancy" in atoms.info:
                    occupancies = atoms.info["_atom_site_occupancy"]
                elif "occupancies" in atoms.info:
                    occupancies = atoms.info["occupancies"]
                else:
                    occupancies = None

                if occupancies:
                    occupancies = np.array(occupancies)
                    atoms.set_array("occupancies", occupancies)

            disorder_components.append(atoms)
        disordered_structure = from_disorder_components(
            disorder_components[0],
            disorder_components[1],
            tolerance=symprec,  # use the same tolerance as for symmetry finder
            symprec=symprec,
            ratio=ratio,  # TODO make this more general
            group_occupancies=None,
        )
        logger.info(disordered_structure)
    else:
        raise ValueError(f"Expected 1 or 2 files, got {nfiles}.")

    symops = disordered_structure.spacegroup.get_symop()
    logger.info("Enumerating ordered configurations.")
    # enumerate ordered configurations
    od = OrderedfromDisordered(disordered_structure, quiet=quiet)
    images, configs = od.get_supercell_configs(
        supercell=supercell,
        maxiters=maxiters,
        exclude_ordered=exclude_ordered,
        random_configs=random,
        return_configs=True,
        fix_ratio=fix_ratio,
        ratio=ratio,
        ratio_tolerance=ratio_tol,
    )
    if len(images) == 0:
        logger.warning("No ordered configurations found.")
        return

    # -- DATA FRAME -- #
    # store information to dataframe
    df["Unique index"] = np.arange(len(images))
    # start out with equal weights
    df["Multiplicity"] = np.ones(len(images), dtype=int)
    df["Configuration"] = configs

    # count the number of each group
    ratios = [np.bincount(config, minlength=np.max(configs) + 1) for config in configs]
    # ratios = [1-sum(config) / len(config) for config in configs]
    # TODO should we define this as the reverse?
    df["Ratio"] = ratios

    df["Supercell"] = "x".join(map(str, supercell))
    # TODO: units
    df["Free Energy / kJ/mol"] = np.NaN
    df["Enthalpy per atom / kJ/mol"] = np.NaN
    df["Entropy per atom / kJ/mol"] = np.NaN
    # volume of simulation cell
    df["Volume/A^3"] = [image.get_volume() for image in images]
    # number of atoms
    df["N atoms"] = [len(image) for image in images]
    # density
    masses = [sum(image.get_masses()) for image in images]
    density_units = (_amu / 1e-30) * 1e-3  # > g/cm^3
    density = [
        (mass / volume) * density_units
        for mass, volume in zip(masses, df["Volume/A^3"])
    ]
    df["Density/g/mol^3"] = density
    # /-- DATA FRAME -- #

    # merge symmetrically equivalent configurations?
    if merge:
        oxidation_states = {}
        if algo == "ewald":
            if ox:
                for ox_name, ox_val in ox:
                    oxidation_states[ox_name] = ox_val
            logger.info(f"Oxidation states used: {oxidation_states}")

        logger.info("Merging structures...")
        start = time.time()
        groups, group_indices = merge_structures(
            images,
            algo=algo,
            symops=symops,
            use_disordered_only=use_disordered_only,
            symprec=symprec,
            oxidation_states=oxidation_states,
            return_group_indices=True,
            quiet=quiet,
            check_species=not ignore_species,
        )

        stop = time.time()
        duration = stop - start
        logger.info(f"Merging took {duration:8.2f} s and found {len(groups)} groups")
        unique_spacegroups = [spglib.get_spacegroup(atoms_to_spglib_cell(g[0])) for g in groups]
        for sg, g in zip(unique_spacegroups, groups):
            logger.info(f"Spacegroup: {sg:<20} multiplicity: {g[1]}")

        # overwrite images with merged images
        images = [g[0] for g in groups]

        # tag images with group indices in df
        for ig, group_inds in enumerate(group_indices):
            for g in group_inds:
                df.loc[g, "Unique index"] = ig

        # Update dataframe with merged information

        # For the Supercell and Ratio columns, we will assume they are the same for each group
        # so we can just take the first value.
        # For the Configuration column, we will take the first value since that corresponds to the
        # configuration of the first image in the group.
        # For the other columns, we will take the mean
        agg_rules = {"Supercell": "first", "Ratio": "first", "Configuration": "first"}
        for col in df.columns:
            if col not in agg_rules:
                agg_rules[col] = np.mean
        # remove 'Unique index' from agg_rules
        agg_rules.pop("Unique index")

        # merge df by Unique index
        df = df.groupby("Unique index").agg(agg_rules).reset_index()
        # add multiplicity to df
        df["Multiplicity"] = [g[1] for g in groups]
        # add spacegroup to df
        df["Spacegroup"] = unique_spacegroups
        # add groups to df
        df["Group indices"] = group_indices

    # log the df
    logging.info("\n\t" + df.to_string().replace("\n", "\n\t"))
    # save the df to csv
    df.to_csv(prefix + ".csv", index=False)
    if not no_write:
        # write structures to file(s)
        directory = f"{prefix}-results"
        logger.info(f"Writing structures to file(s) in directory: {directory}")
        # check if directory already exists -- if so rename it (backup)
        if os.path.isdir(directory):
            logger.warning(
                f"Directory {directory} already exists. Backing up existing directory."
            )
            # backup existing directory
            backup_dir = f"{directory}-backup-{time.strftime('%Y%m%d-%H%M%S')}"
            if os.path.isdir(backup_dir):
                raise ValueError(
                    f"Backup directory {backup_dir} already exists. Please remove it."
                )
            os.rename(directory, backup_dir)

        os.mkdir(directory)
        nleading_zeros = len(str(len(images)))
        for i, image in enumerate(images):
            filename = f"{prefix}-{str(i).zfill(nleading_zeros)}.{format}"
            kwargs = {}
            # cif can save labels if we manually specify them
            if format == "cif":
                kwargs["labels"] = [image.get_array("labels")]

            elif format == "cell":
                # castep cell writer looks for the castep_labels array
                # so let's copy across the labels array
                image.set_array("castep_labels", image.get_array("labels"))
            elif format == "xyz":
                # I had some issues if there are occupancies in the xyz file
                # so let's remove them
                image.set_array("occupancies", None)

            image.write(f"{directory}/{filename}", **kwargs)

            logger.info(f"Wrote structure to file: {directory}/{filename}")

    # write out pandas dataframe containing the groups and multiplicities

    if view:
        # view images using ASE gui
        from ase.visualize import view

        view(images)

    return 0


@cli.command(context_settings={"show_default": True})
@click.pass_context
@click.argument("csv_file", type=click.Path(exists=True))
@click.option(
    "--prefix",
    type=click.STRING,
    default="orgdisord",
    help="Prefix for output file names.",
)
# add option to specify verbosity
@click.option("--quiet", "-q", is_flag=True, default=False, help="Suppress output?")
# add option to specify log file
@click.option(
    "--log",
    "-l",
    type=click.Path(exists=False),
    default="orgdisord.analyse.log",
    help="Log file name.",
)
# specify temperature range and number of steps
@click.option(
    "--temperatures",
    "-t",
    nargs=2,
    type=click.FLOAT,
    default=[10, 300],
    help="Temperature range in K (inclusive of end points).",
)
@click.option(
    "--steps", "-s", type=click.INT, default=50, help="Number of temperature steps."
)
# or set dT
@click.option(
    "--dt", type=click.FLOAT, default=None, help="Temperature step size in K."
)
def analyse(
    ctx,
    csv_file,
    quiet,
    prefix,
    temperatures,
    steps,
    dt,
    log,
):
    """Command line interface for orgdisord thermodynamics.

    TODO energy units!
    TODO: slight numerical inconsistencies wrt to Jonas' code
    """

    # set up logging
    logging.basicConfig(
        filename=log, level=logging.DEBUG, format="%(asctime)s \t %(message)s"
    )
    # add console handler
    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    consolefmt = logging.Formatter("%(message)s")
    console.setFormatter(consolefmt)
    logger = logging.getLogger("orgdisord.analyse")
    # add the handler to the main logger
    logger.addHandler(console)
    logger.info(LOGHEADER)

    # log the Click command line arguments and options used
    logger.debug("Run orgdisord analyse with these command line arguments: ")
    for param, value in ctx.params.items():
        logger.debug(f"          {param}: {value}")

    # read in csv file
    df = pd.read_csv(csv_file)

    # log the number of rows and columns
    logger.info(f"Read in {len(df)} rows and {len(df.columns)} columns from {csv_file}")

    # check that the df contains the required columns
    required_columns = [
        "Multiplicity",
        "Free Energy / kJ/mol",
        "formula units per cell",
    ]
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in {csv_file}")
    # check that if df doens't contain the ratios it at least contains the configuration
    if "Ratio of maj:min" not in df.columns:
        if "Configuration" not in df.columns:
            raise ValueError(f"Neither Ratio nor Configuration found in {csv_file}")
        else:
            logger.debug(
                "Computing the ratios from the configuration specification in the csv file."
            )
            # compute the ratio from the configuration
            # convert to list of ints
            configs = [c.strip("(").strip(")").split(",") for c in df["Configuration"]]
            configs = [[int(c) for c in c] for c in configs]
            ratios = [sum(config) / len(config) for config in configs]
            df["Ratio of maj:min"] = ratios

    # check that there are no missing values
    if df.isnull().values.any():
        raise ValueError(f"Missing values found in {csv_file}")

    # Make new df for the thermodynamic properties
    df_thermo = pd.DataFrame()
    ratios = df["Ratio of maj:min"]
    logger.debug(ratios)
    multiplicities = df["Multiplicity"]
    n_formula_units = df["formula units per cell"]
    energies = df["Free Energy / kJ/mol"] * kJ / mol  # convert to eV
    energies = energies / n_formula_units  # convert to eV/f.u.
    minE = min(energies)
    rel_energies = energies - minE

    # array of temperatures
    if dt:
        # (include end points)
        temperatures = np.arange(temperatures[0], temperatures[1] + dt, dt)
    else:
        temperatures = np.linspace(temperatures[0], temperatures[1], steps)
    # loop over temperatures
    Zs = []
    taus = []
    Ss = []
    Elatts = []
    Us = []
    delta_As = []
    for T in temperatures:
        # partition functions
        Z = get_partition_function(T, multiplicities, rel_energies)
        probabilities = np.array(
            [
                get_probability(multiplicity, energy, T, Z)
                for (multiplicity, energy) in zip(multiplicities, rel_energies)
            ]
        )
        tau = get_tau(ratios, probabilities)
        S = get_S(probabilities)
        Elatt = get_Elatt(probabilities, rel_energies, n_formula_units)
        U = Elatt  # TODO streamline this!
        deltaA = get_deltaA(U, T, S)

        Zs.append(Z)
        taus.append(tau)
        # convert back to kJ/mol
        Ss.append(S * 1e3 * (mol / kJ))
        Elatts.append(Elatt * (mol / kJ))
        Us.append(U * (mol / kJ))
        delta_As.append(deltaA * (mol / kJ))

    # add to df_thermo
    df_thermo["Temperature"] = temperatures
    df_thermo["Z"] = Zs
    df_thermo["tau"] = taus
    df_thermo["S_config (J/mol /K)"] = Ss
    # df_thermo['E_latt (kJ/mol per f.u.)'] = Elatts # TODO do we need both?
    df_thermo["Delta U (kJ/mol per f.u.)"] = Us
    df_thermo["Delta A (kJ/mol per f.u.)"] = delta_As

    # set df print options number of decimal places
    pd.set_option("display.precision", 5)
    logger.info(df_thermo)
    df_thermo.to_csv(prefix + "_thermo.csv", index=False)


def get_boltzmann_weight(multiplicity, energy, temperature):
    """Calculate the Boltzmann weight for a given temperature."""
    return multiplicity * np.exp(-energy / (kB * temperature))


def get_partition_function(temperature, multiplicities, energies):
    """Calculate the partition function for a given temperature.
    Args:
        temperature: temperature in K
        multiplicities: array/list of multiplicities (e.g. [2,2,8,2])
        energies: array/list of energies in kJ/mol

    """
    return np.sum(get_boltzmann_weight(multiplicities, energies, temperature))


def get_probability(multiplicity, energy, temperature, Z):
    """Calculate the probability of a configuration
    for a given temperature and partition function, Z."""
    return get_boltzmann_weight(multiplicity, energy, temperature) / Z


def get_tau(ratios, probabilities):
    """Calculate the tau factor for a given temperature."""
    return np.sum(np.array(ratios) * np.array(probabilities))


def get_S(probabilities):
    """Calculate the entropy S factor for
    a given temperature (implicit in probabilities)."""
    return -kB * np.sum(np.array(probabilities) * np.log(probabilities))


def get_Elatt(probabilities, energies, n_formula_units):
    """Calculate the lattice energy Elatt for a given temperature."""
    return np.sum(np.array(probabilities) * np.array(energies))


def get_deltaA(U, T, S):
    """Calculate the change in free energy deltaA for a given temperature.
    Args:
        U: lattice energy in kJ/mol
        T: temperature in K
        S: entropy in kJ/mol/K
    """
    return U - (T * S)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
