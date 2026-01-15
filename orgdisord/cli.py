"""Console script for orgdisord."""
import sys
import os 
import time
import logging
from typing import List, Optional, Tuple, Any
from ase import Atoms
import pandas as pd
import numpy as np
import typer
from pathlib import Path
from ase.visualize import view as ase_view
from orgdisord.disordered_structure import from_disorder_components
from orgdisord.parse_cif_file import CifParser
from orgdisord.enumerate import OrderedfromDisordered
from orgdisord.merge import merge_structures
from orgdisord.utils import (
    get_new_labels,
    get_partition_function,
    get_probability,
    get_tau,
    get_S,
    get_Elatt,
    get_deltaA,
)
import spglib
from ase.io import read
from ase.units import _amu
from ase.units import kB, kJ, mol

logging.captureWarnings(True)

# new log header
LOGHEADER = """
###############################################################################
#                          -----  orgdisord  -----                            #
#                                                          J. Kane Shenton    #
#                                                                   CCP-NC    #
###############################################################################
"""

app = typer.Typer(help="Tool for handling disordered molecular crystals")

def validate_format(value: str) -> str:
    """Validate and normalize the format string."""
    valid_formats = ["xyz", "cif", "poscar", "cell"]
    value = value.lower()
    if value not in valid_formats:
        raise typer.BadParameter(f"Format must be one of {', '.join(valid_formats)}")
    return value

@app.command(name='enumerate')
def enumerate_configs(
    structure_files: List[Path] = typer.Argument(..., help="Structure files to process", exists=True),
    supercell: Tuple[int, int, int] = typer.Option((1, 1, 1), "--supercell", "-s", help="Supercell size e.g. 2 1 1"),
    maxiters: int = typer.Option(512, "--maxiters", "-N", help="Maximum number of configurations generated"),
    exclude_ordered: bool = typer.Option(False, help="Exclude ordered part of structure in output structures?"),
    random: bool = typer.Option(False, "--random", "-r", help="Generate random configurations?"),
    merge: bool = typer.Option(False, "--merge", "-m", help="Merge equivalent structures?"),
    algo: str = typer.Option("symm", "--algo", "-a", help="Algorithm used for merging equivalent structures.",
                          callback=lambda x: x.lower()),
    use_disordered_only: bool = typer.Option(True, "--use-disordered-only", "-d", 
                                          help="Use only the disordered part of structure in merging?"),
    symprec: float = typer.Option(1e-4, "--symprec", "-p", help="Symmetry precision used for merging equivalent structures."),
    no_write: bool = typer.Option(False, help="Suppress writing out structures?"),
    prefix: str = typer.Option("orgdisord", help="Prefix for output file names."),
    format: str = typer.Option("xyz", "--format", "-f", help="Format of output structure files.",
                            callback=validate_format),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress output?"),
    log: Path = typer.Option("orgdisord.log", "--log", "-l", help="Log file name."),
    view_structures: bool = typer.Option(False, "--view", help="View structures using ASE gui?"),
    ox_species: Optional[List[str]] = typer.Option(None, "--ox-species", help="Species for oxidation states"),
    ox_values: Optional[List[int]] = typer.Option(None, "--ox-values", help="Oxidation state values"),
    molecular_crystal: bool = typer.Option(True, help="Treat the system as a molecular crystal."),
    correlated_assemblies: bool = typer.Option(False, "--correlated-assemblies", "-c", 
                                           help="Are the disorder assemblies correlated?"),
    fix_ratio: bool = typer.Option(False, "--fix-ratio", help="Fix the ratio of the disorder groups."),
    ratio: Optional[float] = typer.Option(None, help="Manually specify the ratio of the disorder groups."),
    ratio_tol: float = typer.Option(0.01, help="Tolerance for ratio comparison."),
    ignore_species: bool = typer.Option(False, "--ignore-species", help="Ignore species when merging structures?")
):
    """Enumerate ordered configurations from a disordered structure."""
    # Set up logging
    logging.basicConfig(
        filename=str(log),
        level=logging.DEBUG,
        format="%(asctime)s \t %(message)s",
    )
    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.DEBUG)

    consolefmt = logging.Formatter("%(message)s")
    console.setFormatter(consolefmt)
    logger = logging.getLogger("orgdisord")
    logger.handlers = []
    logger.addHandler(console)
    logger.info(LOGHEADER)

    # Validate input files count early - this must be outside any try block
    nfiles = len(structure_files)
    if nfiles not in (1, 2):
        logger.error(f"Expected 1 or 2 input structure files, got {nfiles}.")
        raise typer.Exit(code=1)

    try:
        # Load and parse input structures
        if nfiles == 1:
            disordered_structure, symops = _load_single_structure(
                structure_files[0],
                correlated_assemblies=correlated_assemblies,
                molecular_crystal=not molecular_crystal,
                symprec=symprec,
                logger=logger
            )
        else:  # nfiles == 2
            disordered_structure, symops = _load_two_structures(
                structure_files[0],
                structure_files[1],
                symprec=symprec,
                ratio=ratio,
                logger=logger
            )
            
    except Exception as e:
        logger.error(f"Failed to load structures: {str(e)}")
        raise typer.Exit(code=1)

    # Generate configurations
    try:
        logger.info("Enumerating ordered configurations...")
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
    except Exception as e:
        logger.error(f"Failed to enumerate configurations: {str(e)}")
        raise typer.Exit(code=1)
        
    if not images:
        logger.error("No valid configurations were generated.")
        raise typer.Exit(code=1)

    logger.info(f"Generated {len(images)} configurations.")

    # Process configurations
    try:
        df = _create_structure_dataframe(images, configs)
    except Exception as e:
        logger.error(f"Failed to create structure data: {str(e)}")
        raise typer.Exit(code=1)

    # Merge structures if requested
    merged_groups = None
    if merge:
        try:
            merged_images, merged_groups = _merge_structures(
                images,
                algo=algo,
                symops=symops,
                use_disordered_only=use_disordered_only,
                symprec=symprec,
                ox_species=ox_species,
                ox_values=ox_values,
                quiet=quiet,
                logger=logger,
                ignore_species=ignore_species
            )
            images = [g[0] for g in merged_images]
        except Exception as e:
            logger.error(f"Failed to merge structures: {str(e)}")
            raise typer.Exit(code=1)

    # Write output files
    if not no_write:
        try:
            output_path = _get_output_path(prefix, nfiles, merge)
            output_dir = _write_output_files(
                images=images,
                output_path=output_path,
                format=format,
                merged_groups=merged_groups,
                df=df,
                logger=logger
            )
            logger.info(f"Results written to {output_dir}")
        except Exception as e:
            logger.error(f"Failed to write output files: {str(e)}")
            raise typer.Exit(code=1)

    # View structures if requested
    if view_structures:
        try:
            view_my_structures(images)
        except Exception as e:
            logger.error(f"Failed to view structures: {str(e)}")
            # Don't exit on view failure, as it's not critical

    return 0

def _load_single_structure(file_path: Path, correlated_assemblies: bool, 
                         molecular_crystal: bool, symprec: float, logger) -> Tuple:
    """Load a single CIF file."""
    parser = CifParser(str(file_path),
                      correlated_assemblies=correlated_assemblies,
                      molecular_crystal=not molecular_crystal,
                      symprec=symprec)
    disordered_structure = parser.get_disordered_structure()
    symops = disordered_structure.spacegroup.get_symop()
    logger.info("CIF file loaded.")
    return disordered_structure, symops

def _load_two_structures(file_path1: Path, file_path2: Path, symprec: float, 
                        ratio: Optional[float], logger) -> Any:
    """Load two structure files representing major and minor components."""
    try:
        logger.debug(f"Loading major component from {file_path1}")
        atoms_maj = read(str(file_path1))
        logger.debug(f"Major component loaded with {len(atoms_maj)} atoms")
        
        # Ensure arrays are set
        if not atoms_maj.has('labels'):
            atoms_maj.set_array('labels', np.array([f"{sym}{i}" for i, sym in enumerate(atoms_maj.get_chemical_symbols())]))
        if not atoms_maj.has('tags'):
            atoms_maj.set_array('tags', np.zeros(len(atoms_maj), dtype=int))

        logger.debug(f"Loading minor component from {file_path2}")
        atoms_min = read(str(file_path2))
        logger.debug(f"Minor component loaded with {len(atoms_min)} atoms")
        
        # Ensure arrays are set
        if not atoms_min.has('labels'):
            atoms_min.set_array('labels', np.array([f"{sym}{i}" for i, sym in enumerate(atoms_min.get_chemical_symbols())]))
        if not atoms_min.has('tags'):
            atoms_min.set_array('tags', np.zeros(len(atoms_min), dtype=int))
        
        logger.debug("Creating disordered structure from components")
        disordered_structure = from_disorder_components(
            atoms_maj,
            atoms_min,
            tolerance=symprec,
            symprec=symprec,
            ratio=ratio,
            group_occupancies=None,
        )

        logger.debug("Disordered structure created successfully")
        logger.info(str(disordered_structure))
        symops = disordered_structure.spacegroup.get_symop()
        return disordered_structure, symops
    except Exception as e:
        import traceback
        logger.error(f"Failed to load or process structures: {str(e)}")
        logger.debug(f"Traceback: {traceback.format_exc()}")
        raise

def _create_structure_dataframe(images: List, configs: List) -> pd.DataFrame:
    """Create a DataFrame containing structure information."""
    df = pd.DataFrame()
    df["Config_idx"] = list(range(len(images)))
    df["Config"] = [str(c) for c in configs]
    df["Volume/A^3"] = [float(image.get_volume()) for image in images]

    # Calculate masses and densities
    masses = [sum(image.get_masses()) for image in images]
    density_units = (_amu / 1e-30) * 1e-3  # Convert to g/cm^3
    density = [
        float((mass / volume) * density_units)
        for mass, volume in zip(masses, df["Volume/A^3"])
    ]
    df["Density/g/cm^3"] = density
    return df

def _merge_structures(images: List, algo: str, symops: Any, use_disordered_only: bool,
                     symprec: float, ox_species: Optional[List[str]], 
                     ox_values: Optional[List[int]], quiet: bool, logger, ignore_species: bool = False) -> Tuple[Tuple[List[Atoms], List], List]:
    """Merge equivalent structures."""
    oxidation_states = {}
    if ox_species and ox_values:
        oxidation_states = dict(zip(ox_species, ox_values))
        logger.info(f"Oxidation states used: {oxidation_states}")

    logger.info("Merging structures...")
    start = time.time()
    try:
        merged_result = merge_structures(
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
        logger.info(f"Merged into {len(merged_result[1])} groups.")
        end = time.time()
        logger.info(f"Time taken: {end-start:.2f} s")
        return merged_result
    except Exception as e:
        logger.error(f"Failed to merge structures: {str(e)}")
        raise

def _get_output_path(prefix: str, nfiles: int, merge: bool) -> str:
    """Generate the output path based on input parameters."""
    output_path = prefix
    if nfiles == 2:
        output_path += "_from_maj_min"
    if merge:
        output_path += "_merged"
    return output_path

def _write_output_files(images: List[Atoms], output_path: str, format: str,
                       merged_groups: Optional[List], df: pd.DataFrame,
                       logger) -> str:
    """Write enumerated structures and associated data to files."""
    directory = f"{output_path}-results"
    logger.info(f"Writing structures to directory: {directory}")

    # Handle existing directory
    if os.path.isdir(directory):
        backup_suffix = time.strftime('%Y%m%d-%H%M%S')
        backup_dir = f"{directory}-backup-{backup_suffix}"
        
        # Keep trying backups with incrementing numbers if needed
        backup_num = 1
        while os.path.isdir(backup_dir):
            backup_dir = f"{directory}-backup-{backup_suffix}-{backup_num}"
            backup_num += 1
            
        os.rename(directory, backup_dir)
        logger.info(f"Backed up existing directory to {backup_dir}")

    os.makedirs(directory, exist_ok=True)

    # Write structure files
    n_images = len(images)
    nleading_zeros = len(str(n_images))
    format = format.lower()
    
    for i, image in enumerate(images):
        try:
            filename = f"{output_path}-{str(i).zfill(nleading_zeros)}.{format}"
            filepath = os.path.join(directory, filename)
            
            kwargs = {}
            if format == "cif" and hasattr(image, "arrays") and "labels" in image.arrays:
                kwargs["labels"] = image.get_array("labels")
            elif format == "cell" and hasattr(image, "arrays") and "labels" in image.arrays:
                image.arrays["castep_labels"] = image.get_array("labels")
            elif format == "xyz" and hasattr(image, "arrays") and "occupancies" in image.arrays:
                del image.arrays["occupancies"]
            
            image.write(filepath, **kwargs)
            logger.debug(f"Wrote structure {i+1}/{n_images} to {filepath}")
        except Exception as e:
            logger.error(f"Failed to write structure {i} to {filepath}: {str(e)}")
            continue

    # Write data files
    try:
        csv_path = os.path.join(directory, f"{output_path}.csv")
        df.to_csv(csv_path, index=False)
        logger.info(f"Wrote configuration data to {csv_path}")
        
        if merged_groups:
            # Create mapping from structure index to group
            group_map = {}
            for group_idx, group in enumerate(merged_groups):
                for struct_idx in group:
                    group_map[struct_idx] = group_idx
            
            # Add group information
            df["Group"] = [group_map.get(i, -1) for i in range(len(images))]
            merged_csv_path = os.path.join(directory, f"{output_path}_merged.csv")
            df.to_csv(merged_csv_path, index=False)
            logger.info(f"Wrote merged groups data to {merged_csv_path}")
            
    except Exception as e:
        logger.error(f"Failed to write data files: {str(e)}")
        # Continue even if data files fail, as we still have the structure files

    return directory

def view_my_structures(images, viewer=None):
    """View structures using ASE GUI or specified viewer.
    
    Handles potential XCB/X11 errors gracefully by checking for a valid display
    and setting appropriate environment variables before attempting visualization.
    """
    logger = logging.getLogger("orgdisord")
    
    # Check if we have a valid display
    if not os.environ.get('DISPLAY'):
        logger.error("No display available. Cannot open ASE GUI viewer.")
        logger.error("Try running with --no-view option if you are in a headless environment")
        return
        
    try:
        # Set environment variables to avoid XCB threading issues
        os.environ['ASE_USE_GLIB'] = 'false'  # Disable GLib integration
        os.environ['XCB_ALLOW_THREADS'] = '0'  # Disable XCB threading
        os.environ['QT_X11_NO_MITSHM'] = '1'  # Disable MIT-SHM extension
        
        if viewer:
            ase_view(images, viewer=viewer, block=True)
        else:
            ase_view(images, block=True)
            
    except (AssertionError, IOError) as e:
        if "xcb" in str(e).lower():
            logger.error("X11/XCB threading error occurred while trying to open the viewer.")
            logger.error("This is a known issue with X11 forwarding and threading.")
            logger.error("Try one of these solutions:")
            logger.error("1. Run with --no-view option")
            logger.error("2. Use a different viewer: --viewer 'povray' or 'matplotlib'")
            logger.error("3. Set XCB_ALLOW_THREADS=0 in your environment")
        else:
            logger.error(f"Failed to open viewer: {str(e)}")
    except Exception as e:
        logger.error(f"Unexpected error while opening viewer: {str(e)}")

@app.command()
def analyse(
    csv_file: Path = typer.Argument(..., help="CSV file to analyze", exists=True),
    prefix: str = typer.Option("orgdisord", help="Prefix for output file names."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Suppress output?"),
    log: Path = typer.Option("orgdisord.analyse.log", "--log", "-l", help="Log file name."),
    temperatures: Tuple[float, float] = typer.Option((10, 300), "--temperatures", "-t", 
                                                  help="Temperature range in K (inclusive of end points)."),
    steps: int = typer.Option(50, "--steps", "-s", help="Number of temperature steps."),
    dt: Optional[float] = typer.Option(None, help="Temperature step size in K.")
):
    """Analyse enumeration results.
    
    Calculate thermodynamic properties and analyze disorder patterns from
    enumerated configurations.
    """
    # set up logging
    logging.basicConfig(
        filename=str(log),
        level=logging.DEBUG,
        format="%(asctime)s \t %(message)s"
    )
    console = logging.StreamHandler()
    if quiet:
        console.setLevel(logging.WARNING)
    else:
        console.setLevel(logging.INFO)

    consolefmt = logging.Formatter("%(message)s")
    console.setFormatter(consolefmt)
    logger = logging.getLogger("orgdisord.analyse")
    logger.handlers = []  # Clear any existing handlers
    logger.addHandler(console)
    logger.info(LOGHEADER)

    logger.debug("Run orgdisord analyse with these command line arguments: ")
    kwargs = locals()
    for param, value in kwargs.items():
        logger.debug(f"          {param}: {value}")

    # read in csv file
    df = pd.read_csv(csv_file)
    logger.info(f"Read in {len(df)} rows and {len(df.columns)} columns from {csv_file}")

    # check required columns
    required_columns = [
        "Multiplicity",
        "Free Energy / kJ/mol",
        "formula units per cell",
    ]
    for col in required_columns:
        if col not in df.columns:
            logger.error(f"Required column {col} not found in {csv_file}")
            raise typer.Exit(code=1)

    # check ratios
    if "Ratio of maj:min" not in df.columns:
        if "Configuration" not in df.columns:
            logger.error(f"Neither Ratio nor Configuration found in {csv_file}")
            raise typer.Exit(code=1)
        else:
            logger.debug(
                "Computing the ratios from the configuration specification in the csv file."
            )
            configs = [c.strip("(").strip(")").split(",") for c in df["Configuration"]]
            configs = [[int(c) for c in c] for c in configs]
            ratios = [sum(config) / len(config) for config in configs]
            df["Ratio of maj:min"] = ratios

    # check for missing values
    if df.isnull().values.any():
        logger.error(f"Missing values found in {csv_file}")
        raise typer.Exit(code=1)

    # Make new df for thermodynamic properties
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
    temps_start, temps_end = temperatures
    if dt:
        # (include end points)
        temperatures = np.arange(temps_start, temps_end + dt, dt)
    else:
        temperatures = np.linspace(temps_start, temps_end, steps)

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
    df_thermo["Delta U (kJ/mol per f.u.)"] = Us
    df_thermo["Delta A (kJ/mol per f.u.)"] = delta_As

    # set df print options number of decimal places
    pd.set_option("display.precision", 5)
    logger.info(df_thermo)
    df_thermo.to_csv(prefix + "_thermo.csv", index=False)

def cli():
    """Entry point for the CLI."""
    app()
