# https://sites.google.com/site/orcainputlibrary/geometry-optimizations
# Use GGA DFT functionals if they are accurate enough (depends on your system), with the RI-J approximation (default)
# as that is often the fastest useful optimization one can do. Use of the RI-J approximation leads to minimal
# geometrical errors. Often the slightly higher accuracy from hybrid functionals is not worth the effort.
import datetime
import glob
import os
import subprocess
import time

import yaml


def mkdir(dir_name: str) -> None:
    """Creates a directory."""
    try:
        os.mkdir(dir_name)
        print(f"Directory '{dir_name}' created successfully.")
    except FileExistsError:
        print(f"Warning: Directory '{dir_name}' already exists.")
    except PermissionError:
        print(f"Warning: Permission denied: Unable to create '{dir_name}'.")
    except Exception as e:
        print(f"Warning: An error occurred: {e}")


def make_inp_from_xyz(xyz_filename: str, inp_destination_path: str, job_type: str, RI: str,
                      functional: str="BP86", basis_set: str="def2-SVP", dispersion_correction:str ="D3BJ",
                      grid:str ="", charge: int=0, freq: bool=False, NMR: bool=False, cores=6) -> None:
    """Produces an ORCA input file for geom. opt. from xyz file."""
    # The choices of keywords here are from: https://sites.google.com/site/orcainputlibrary/geometry-optimizations
    # RI: RI-J approximation for Coulomb integrals: speed calculations at cost of small error, used for GGA calcs.
    # RIJCOSX should be used with hybrid functionals.
    # D3BJ: Grimme's D3 dispersion correction, with Becke-Johnson damping.
    # TIGHTSCF: convergence tolerance level, this is recommended for geometry optimizations.
    # NormalSCF: used for single point calculations.
    # Here we're using 6 cores by default.

    if job_type == "Geometry Optimization":
        job_keywords = "TIGHTSCF Opt"
    elif job_type == "Single Point Calculation":
        job_keywords = "NormalSCF"
    else:
        raise Exception(f"Unknown job type {job_type}")

    freq = "Freq" if freq else ""
    NMR = "NMR" if NMR else ""
    if grid != "":
        grid = f"\n! {grid}"

    header = (f"! {RI} {functional} {basis_set} {dispersion_correction} {job_keywords} {freq} {NMR}{grid}"
              f"\n%pal\nnprocs {cores}\nend\n\n* xyz {charge} 1\n")
    with open(xyz_filename, 'r') as xyz_file:
        # Remove the initial lines of the xyz file, leaving only the atoms and their coordinates:
        inp_contents = "\n".join(xyz_file.read().splitlines()[2:])

    with open(inp_destination_path, "w") as inp_file:
        inp_file.write(header + inp_contents + "\n*")


def orca_batch_job(path_to_xyz_files: str, destination_path: str, job_type: str, RI: str,
                   functional: str="BP86", basis_set: str="def2-SVP", dispersion_correction:str ="D3BJ", grid: str="",
                   charge: int=0, freq: bool=False, NMR: bool=False, redo_all: bool=False, log_name: str="log") -> None:
    """Performs ORCA calculations based on the xyz files in the given directory."""
    log = f"{job_type} started {datetime.datetime.now()}\n\nCreating directories:\n"

    # Recursively search for the paths to all xyz files anywhere under the given path:
    xyz_file_paths = glob.glob(f"{path_to_xyz_files}/**/*.xyz", recursive=True)
    cids = []

    for xyz_file_path in xyz_file_paths:
        # We are trying to stick to Linux-style path formatting, so replace the Windows \\ with /:
        xyz_file_path = xyz_file_path.replace("\\", "/")
        # Basically, for all intents and purposes, CID herein just means the unique label for this molecule.
        # It wouldn't have to be an actual CID.
        cid = xyz_file_path.split("/")[-1][:-4]

        # Skip trajectory .xyz files, which are produced by geometry opts. They end in _trj.xyz.
        if cid[-4:] == "_trj":
            continue
        cids.append(cid)

        try:
            mkdir(f"{destination_path}/{cid}")
        except Exception as e:
            log += f"{cid}: Directory creation in {destination_path} unsuccessful with exception {e}.\n"
            continue

        make_inp_from_xyz(
            functional=functional,
            basis_set=basis_set,
            dispersion_correction=dispersion_correction,
            RI=RI,
            job_type=job_type,
            grid=grid,
            charge=charge,
            freq=freq,
            NMR=NMR,
            xyz_filename=xyz_file_path,
            inp_destination_path=f"{destination_path}/{cid}/{cid}.inp"
        )

    log += f"\n{job_type} part:\n"
    print()
    for cid in cids:
        if not redo_all:
            # If an existing .out file is found in the directory, check if it seems that it was from a successful calc.
            # If so, skip doing the geometry optimization for that CID.
            if os.path.exists(f"{destination_path}/{cid}/{cid}.out"):
                with open(f"{destination_path}/{cid}/{cid}.out", 'r') as out_file:
                    if out_file.read().splitlines()[-2].strip() == "****ORCA TERMINATED NORMALLY****":
                        msg = (f"Existing .out file from completed calculation found in {destination_path}/{cid}/, "
                               f"skipping calculation.")
                        print(msg)
                        log += msg + "\n"
                        continue
                    else:
                        msg = (f"Existing .out file from apparently unsuccessful calculation found "
                               f"in {destination_path}/{cid}/, redoing calculation.")
                        print(msg)
                        log += msg + "\n"

        print(f"{datetime.datetime.now()} Performing {job_type} on {cid}: ", end="")

        # Load the absolute path to ORCA from config.yaml; this is necessary for calculations run in parallel.
        with open("config.yaml") as f:
            cfg = yaml.load(f, Loader=yaml.FullLoader)
            orca_path = cfg['orca_path']

        start = time.time()
        orca_command = f"{orca_path} {destination_path}/{cid}/{cid}.inp > {destination_path}/{cid}/{cid}.out"
        subprocess.run(orca_command, shell=True)
        print(f"complete. Total time: {datetime.timedelta(seconds=time.time() - start)}")
        log += f"{job_type} for {cid} completed at {datetime.datetime.now()}\n"
        with open(f"logs/{log_name}.txt", "w") as log_file:
            log_file.write(log)

    print(f"\n{job_type}s complete.")
    log += f"\n{job_type} ended {datetime.datetime.now()}\n"
    with open(f"logs/{log_name}.txt", "w") as log_file:
        log_file.write(log)
