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


def make_uv_vis_plot(path_to_tddft_out: str, wavenum_min: int=10000, wavenum_max: int=50000,
                     wavenum_broadening: int=4000) -> None:
    """Use orca_mapspc to generate .abs.dat and .abs.stk files from TDDFT ORCA calculation."""
    # Load the absolute path to ORCA from config.yaml; this is necessary for calculations run in parallel.
    with open("config.yaml") as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
        orca_path = cfg['orca_path']

    orca_command = f"{orca_path}_mapspc {path_to_tddft_out} ABS -x0{wavenum_min} -x1{wavenum_max} -w{wavenum_broadening}"
    subprocess.run(orca_command, shell=True)
    # This produces two files: one that ends in .abs.dat (broadened transitions) and out that ends in .abs.stk.
    # The .abs.dat file has 5 columns which, left to right, are (I believe): frequency, intensity, electric dipole, magnetic dipole, and quadrupole.
    # I am getting that in part from: https://github.com/Ianiusha/Uselfuls_Scripts/blob/master/plot_UVVIS.py
    # The .abs.stk file includes "energies and intensities".
    # Also see: https://www.youtube.com/watch?v=CDQRIDkOO8Q


def make_inp_from_xyz(xyz_filename: str, inp_destination_path: str, job_type: str, RI: str, functional: str,
                      basis_set: str, newgto: str, dispersion_correction: str, solvent: str, grid: str, charge: int=0,
                      freq: bool=False, NMR: bool=False, orbitals: bool=False, cores=6) -> None:
    """Produces an ORCA input file for geom. opt. from xyz file."""
    # The choices of keywords here are from: https://sites.google.com/site/orcainputlibrary/geometry-optimizations
    # RI: RI-J approximation for Coulomb integrals: speed calculations at cost of small error, used for GGA calcs.
    # RIJCOSX should be used with hybrid functionals. NOTE: for functionals other than def2 family, you will need to
    # specify an auxiliary basis set, which can be done automatically: e.g. "autoaux RIJCOSX".
    # D3BJ: Grimme's D3 dispersion correction, with Becke-Johnson damping.
    # TIGHTSCF: convergence tolerance level, this is recommended for geometry optimizations.
    # NormalSCF: used for single point calculations.
    # Here we're using 6 cores by default.
    tddft = ""
    if job_type == "Geometry Optimization":
        job_keywords = "TightSCF Opt"
    elif job_type == "Single Point Calculation":
        job_keywords = "TightSCF"
    elif job_type == "TDDFT Calculation" or job_type == "Single Point Calculation":
        job_keywords = "TightSCF"
        tddft = "\n%tddft\nnroots 40\nmaxdim 5\nend\n"
    else:
        logger.error(f"Unknown job type: {job_type}")
        raise Exception(f"Unknown job type {job_type}")

    freq = "Freq" if freq else ""
    largeprint = "! LARGEPRINT\n" if orbitals else ""
    NMR = "NMR" if NMR else ""
    if grid != "":
        grid = f"! {grid}\n"
    if newgto != "":
        newgto = f"%basis\nnewgto {newgto} end\nend\n"

    header = (f"! {RI} {functional} {basis_set} {dispersion_correction} {job_keywords} {solvent} {freq} {NMR}\n"
              f"{newgto}"
              f"{grid}"
              f"{tddft}"
              f"{largeprint}"
              f"%pal\nnprocs {cores}\nend\n\n* xyz {charge} 1\n")
    with open(xyz_filename, 'r') as xyz_file:
        # Remove the initial lines of the xyz file, leaving only the atoms and their coordinates:
        inp_contents = "\n".join(xyz_file.read().splitlines()[2:])

    with open(inp_destination_path, "w") as inp_file:
        inp_file.write(header + inp_contents + "\n*")


# TODO: change it so path_to_xyz_file contains the xyz_filename and we just parse it out in the function?
def orca_job(path_to_xyz_file: str, xyz_filename_no_extension: str, destination_path: str, job_type: str,
             RI: str, functional: str, basis_set: str, newgto: str, dispersion_correction: str, solvent: str,
             grid: str, charge: int=0, freq: bool=False, NMR: bool=False, orbitals: bool=False) -> None:
    """Performs ORCA calculations based on the inputs and given xyz file."""
    # We are trying to stick to Linux-style path formatting, so replace the Windows \\ with /:
    path_to_xyz_file = path_to_xyz_file.replace("\\", "/")

    if job_type == "Geometry Optimization":
        full_filename = xyz_filename_no_extension + "_geom_opt"
        # Note that xyz_filename should not contain the extension ".xyz"!! It should be the filename part only.
    elif job_type == "Single Point Calculation":
        full_filename = xyz_filename_no_extension + "_single_pt_quick"  # TODO: change back?
    elif job_type == "TDDFT Calculation":
        full_filename = xyz_filename_no_extension + "_tddft"
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
        orbitals=orbitals,
        xyz_filename=path_to_xyz_file,
        inp_destination_path=f"{destination_path}/{full_filename}.inp"
    )

    # If an existing .out file is found in the directory, check if it seems that it was from a successful calc.
    # If so, skip doing the calculation.
    if os.path.exists(f"{destination_path}/{full_filename}.out"):
        with open(f"{destination_path}/{full_filename}.out", 'r') as out_file:

            # TODO: change this so that it uses regex instead of splitting the list?
            # I suspect that we are getting some false positives when it comes to detecting whether a calc didn't finish properly
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

    logger.info(f"Performing {job_type} on {xyz_filename_no_extension}: ")
    start = time.time()
    orca_command = f"{orca_path} {destination_path}/{full_filename}.inp > {destination_path}/{full_filename}.out"
    subprocess.run(orca_command, shell=True)
    logger.info(f"complete. Total time: {datetime.timedelta(seconds=time.time() - start)}")

    if job_type == "TDDFT Calculation":
        logger.info("Producing UV-Vis spectra from TDDFT Calculation.")
        make_uv_vis_plot(path_to_tddft_out=f"{destination_path}/{full_filename}.out")
    logger.info("\n")


def orca_job_sequence(path_to_conf_search_xyz_files: str, destination_path: str,
                      geom_opt_arguments: dict, part_2_arguments: dict,
                      geom_opt: bool=True, single_pt: bool=False, tddft: bool=False) -> None:
    """
    Perform sequential geometry optimization, single point and/or TDDFT calcs from .xyz files in a given directory.

    Mote: single point and TDDFT calculations will be based on the geometry optimization .xyz file, so make sure it
    is there if you are running those calculations. In other words, they won't be done based on the .xyz files from the
    INPUT path, which are assumed to be from a conformer search and therefore not suitable for these calculations.
    """
    # Recursively search for the paths to all xyz files anywhere under the given path:
    xyz_file_paths = glob.glob(f"{path_to_conf_search_xyz_files}/**/*.xyz", recursive=True)

    logger.info(f"List of .xyz files found: {xyz_file_paths}.")
    logger.info(f"Geometry Optimization arguments: {geom_opt_arguments}")
    if tddft:
        logger.info(f"TDDFT Calculation arguments: {part_2_arguments}\n")
    else:
        logger.info(f"Single Point Calculation arguments: {part_2_arguments}\n")

    for xyz_file_path in xyz_file_paths:
        # We are trying to stick to Linux-style path formatting, so replace the Windows \\ with /:
        xyz_file_path = xyz_file_path.replace("\\", "/")
        xyz_filename = xyz_file_path.split("/")[-1][:-4]  # remove the .xyz as well
        mol_id = xyz_filename.split("_")[0]

        mkdir(f"{destination_path}/{mol_id}")

        if geom_opt:
            mkdir(f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt")
            # Geometry optimization calculation from the conformer search .xyz file:
            orca_job(path_to_xyz_file=xyz_file_path, xyz_filename_no_extension=xyz_filename,
                     destination_path=f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt",
                     job_type="Geometry Optimization", RI=geom_opt_arguments["RI"],
                     functional=geom_opt_arguments["functional"], basis_set=geom_opt_arguments["basis_set"],
                     newgto=geom_opt_arguments["newgto"],
                     dispersion_correction=geom_opt_arguments["dispersion_correction"],
                     solvent=geom_opt_arguments["solvent"], grid=geom_opt_arguments["grid"],
                     freq=geom_opt_arguments["freq"], NMR=geom_opt_arguments["NMR"],
                     orbitals=geom_opt_arguments["orbitals"])
            logger.info(f"{mol_id} geometry optimization complete.\n")

        if tddft:
            mkdir(f"{destination_path}/{mol_id}/{xyz_filename}_tddft")
            # Time-dependent DFT calculation from the geometry optimization .xyz file:
            orca_job(
                path_to_xyz_file=f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt/{xyz_filename}_geom_opt.xyz",
                xyz_filename_no_extension=xyz_filename,
                destination_path=f"{destination_path}/{mol_id}/{xyz_filename}_tddft",
                job_type="TDDFT Calculation", RI=part_2_arguments["RI"],
                functional=part_2_arguments["functional"], basis_set=part_2_arguments["basis_set"],
                newgto=part_2_arguments["newgto"],
                dispersion_correction=part_2_arguments["dispersion_correction"],
                solvent=part_2_arguments["solvent"], grid=part_2_arguments["grid"],
                freq=part_2_arguments["freq"], NMR=part_2_arguments["NMR"], orbitals=part_2_arguments["orbitals"])
            logger.info(f"{mol_id} TDDFT calculation complete.\n")

        if single_pt:
            mkdir(f"{destination_path}/{mol_id}/{xyz_filename}_single_pt_quick")  # TODO: change back
            # Single point calculation from the geometry optimization .xyz file:
            orca_job(
                path_to_xyz_file=f"{destination_path}/{mol_id}/{xyz_filename}_geom_opt/{xyz_filename}_geom_opt.xyz",
                xyz_filename_no_extension=xyz_filename,
                destination_path=f"{destination_path}/{mol_id}/{xyz_filename}_single_pt_quick",   # TODO: change back
                job_type="Single Point Calculation", RI=part_2_arguments["RI"],
                functional=part_2_arguments["functional"], basis_set=part_2_arguments["basis_set"],
                newgto=part_2_arguments["newgto"],
                dispersion_correction=part_2_arguments["dispersion_correction"],
                solvent=part_2_arguments["solvent"], grid=part_2_arguments["grid"],
                freq=part_2_arguments["freq"], NMR=part_2_arguments["NMR"], orbitals=part_2_arguments["orbitals"])
            logger.info(f"{mol_id} single point calculation complete.\n\n")

        logger.info("\n")
