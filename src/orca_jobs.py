# https://sites.google.com/site/orcainputlibrary/geometry-optimizations
# Use GGA DFT functionals if they are accurate enough (depends on your system), with the RI-J approximation (default)
# as that is often the fastest useful optimization one can do. Use of the RI-J approximation leads to minimal
# geometrical errors. Often the slightly higher accuracy from hybrid functionals is not worth the effort.
import datetime
import glob
import os
import subprocess
import time
import logging
logger = logging.getLogger(__name__)

import yaml


def mkdir(dir_name: str) -> None:
    """Creates a directory."""
    try:
        os.mkdir(dir_name)
        logger.info(f"Directory '{dir_name}' created successfully.")
    except FileExistsError:
        logger.warning(f"Warning: Directory '{dir_name}' already exists.")
    except PermissionError:
        logger.warning(f"Warning: Permission denied: Unable to create '{dir_name}'.")
    except Exception as e:
        logger.error(f"Warning: An error occurred: {e}")


def make_inp_from_xyz(xyz_filename: str, inp_destination_path: str, job_type: str, RI: str,
                      functional: str="BP86", basis_set: str="def2-SVP", newgto: str="",
                      dispersion_correction:str ="D3BJ", solvent:str ="", grid:str ="", charge: int=0,
                      freq: bool=False, NMR: bool=False, cores=6) -> None:
    """Produces an ORCA input file for geom. opt. from xyz file."""
    # The choices of keywords here are from: https://sites.google.com/site/orcainputlibrary/geometry-optimizations
    # RI: RI-J approximation for Coulomb integrals: speed calculations at cost of small error, used for GGA calcs.
    # RIJCOSX should be used with hybrid functionals.
    # D3BJ: Grimme's D3 dispersion correction, with Becke-Johnson damping.
    # TIGHTSCF: convergence tolerance level, this is recommended for geometry optimizations.
    # NormalSCF: used for single point calculations.
    # Here we're using 6 cores by default.
    if job_type == "Geometry Optimization":
        job_keywords = "TightSCF Opt"
    elif job_type == "Single Point Calculation":
        job_keywords = "NormalSCF"
    else:
        logger.error(f"Unknown job type: {job_type}")
        raise Exception(f"Unknown job type {job_type}")

    freq = "Freq" if freq else ""
    NMR = "NMR" if NMR else ""
    if grid != "":
        grid = f"! {grid}\n"
    if newgto != "":
        newgto = f"%basis\nnewgto {newgto} end\nend\n"

    header = (f"! {RI} {functional} {basis_set} {dispersion_correction} {job_keywords} {solvent} {freq} {NMR}\n"
              f"{newgto}"
              f"{grid}"
              f"%pal\nnprocs {cores}\nend\n\n* xyz {charge} 1\n")
    with open(xyz_filename, 'r') as xyz_file:
        # Remove the initial lines of the xyz file, leaving only the atoms and their coordinates:
        inp_contents = "\n".join(xyz_file.read().splitlines()[2:])

    with open(inp_destination_path, "w") as inp_file:
        inp_file.write(header + inp_contents + "\n*")


def orca_job(path_to_xyz_file: str, xyz_filename: str, destination_path: str, job_type: str,
             RI: str, functional: str="BP86", basis_set: str="def2-SVP", newgto: str="",
             dispersion_correction:str ="D3BJ", solvent:str ="", grid: str="", charge: int=0, freq: bool=False,
             NMR: bool=False) -> None:
    """Performs ORCA calculations based on the inputs and given xyz file."""
    # We are trying to stick to Linux-style path formatting, so replace the Windows \\ with /:
    path_to_xyz_file = path_to_xyz_file.replace("\\", "/")

    if job_type == "Geometry Optimization":
        full_filename = xyz_filename + "_geom_opt"
        # Note that xyz_filename should not contain the extension ".xyz"!! It should be the filename part only.
    elif job_type == "Single Point Calculation":
        full_filename = xyz_filename + "_single_pt"
    else:
        logger.error(f"Unknown job type: {job_type}")
        raise Exception(f"Unknown job type {job_type}")

    make_inp_from_xyz(
        functional=functional,
        basis_set=basis_set,
        newgto=newgto,
        dispersion_correction=dispersion_correction,
        RI=RI,
        job_type=job_type,
        grid=grid,
        solvent=solvent,
        charge=charge,
        freq=freq,
        NMR=NMR,
        xyz_filename=path_to_xyz_file,
        inp_destination_path=f"{destination_path}/{full_filename}.inp"
    )

    # If an existing .out file is found in the directory, check if it seems that it was from a successful calc.
    # If so, skip doing the calculation.
    if os.path.exists(f"{destination_path}/{full_filename}.out"):
        with open(f"{destination_path}/{full_filename}.out", 'r') as out_file:
            if out_file.read().splitlines()[-2].strip() == "****ORCA TERMINATED NORMALLY****":
                logger.info(f"Existing .out file from completed calculation found in {destination_path}, "
                            f"skipping calculation.")
                return
            else:
                logger.warning(f"Existing .out file from apparently unsuccessful calculation found "
                               f"in {destination_path}, redoing calculation.")

    # Load the absolute path to ORCA from config.yaml; this is necessary for calculations run in parallel.
    with open("config.yaml") as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        orca_path = cfg['orca_path']

    logger.info(f"Performing {job_type} on {xyz_filename}: ")
    start = time.time()
    orca_command = f"{orca_path} {destination_path}/{full_filename}.inp > {destination_path}/{full_filename}.out"
    subprocess.run(orca_command, shell=True)
    logger.info(f"complete. Total time: {datetime.timedelta(seconds=time.time() - start)}\n")


def orca_job_sequence(path_to_conf_search_xyz_files: str, destination_path: str,
                      geom_opt_arguments: dict, single_pt_arguments: dict, log_path="logs/log.log") -> None:
    """
    Perform sequential geometry optimization and single point calculations based on .xyz files in a given directory.
    """
    # Recursively search for the paths to all xyz files anywhere under the given path:
    xyz_file_paths = glob.glob(f"{path_to_conf_search_xyz_files}/**/*.xyz", recursive=True)

    logger.info(f"List of .xyz files found: {xyz_file_paths}.")
    logger.info(f"Geometry Optimization arguments: {geom_opt_arguments}")
    logger.info(f"Single Point Calculation arguments: {single_pt_arguments}\n")

    for xyz_file_path in xyz_file_paths:
        # We are trying to stick to Linux-style path formatting, so replace the Windows \\ with /:
        xyz_file_path = xyz_file_path.replace("\\", "/")
        xyz_filename = xyz_file_path.split("/")[-1][:-4]  # remove the .xyz as well
        mol_id = xyz_filename.split("_")[0]

        mkdir(f"{destination_path}/{mol_id}")
        mkdir(f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt")
        mkdir(f"{destination_path}/{mol_id}/{xyz_filename}_single_pt")

        # Geometry optimization:
        orca_job(
            path_to_xyz_file=xyz_file_path,
            xyz_filename=xyz_filename,
            destination_path=f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt",
            job_type="Geometry Optimization",
            functional=geom_opt_arguments["functional"],
            basis_set=geom_opt_arguments["basis_set"],
            newgto=geom_opt_arguments["newgto"],
            RI=geom_opt_arguments["RI"],
            dispersion_correction=geom_opt_arguments["dispersion_correction"],
            solvent=geom_opt_arguments["solvent"]
        )

        # Single point calculation from the geometry optimization .xyz file:
        orca_job(
            path_to_xyz_file=f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt/{xyz_filename}_geom_opt.xyz",
            xyz_filename=xyz_filename,
            destination_path=f"{destination_path}/{mol_id}/{xyz_filename}_single_pt",
            job_type="Single Point Calculation",
            functional=single_pt_arguments["functional"],
            basis_set=single_pt_arguments["basis_set"],
            newgto=single_pt_arguments["newgto"],
            RI=single_pt_arguments["RI"],
            dispersion_correction=single_pt_arguments["dispersion_correction"],
            solvent=single_pt_arguments["solvent"],
            NMR=single_pt_arguments["NMR"],
            freq=single_pt_arguments["freq"]
        )
        logger.info("Geometry optimization and single point calculation steps complete.\n\n")
