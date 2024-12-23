import glob
import re

import pandas as pd


# TODO: work on this...
def single_pt_E_workup(path_to_single_pt_out_files: str, destination_path: str,
                       skip_atom46: bool=True, skip_tddft: bool=False) -> None:
    """
    Finds all .out files recursively in a given directory and checks if they have single point calculations.
    For the ones that do, their single point energies are compiled into an Excel spreadsheet.
    """
    single_pt_out_file_paths = glob.glob(f"{path_to_single_pt_out_files}/**/*.out", recursive=True)
    single_pt_Es = {}

    for single_pt_out_file in single_pt_out_file_paths:
        single_pt_out_file_path = single_pt_out_file.replace("\\", "/")

        calc_name = single_pt_out_file_path.split("/")[-1][:-4]  # remove .out

        # Skip those .out files that end in atom46:
        try:
            if skip_atom46 and calc_name[-6:] == "atom46":
                continue
        except IndexError:
            pass

        with open(single_pt_out_file_path, 'r') as out_file:
            outfile_contents = out_file.read()

        # If told to do so, determine if the file is a TD-DFT calculation and skip it if so:
        if skip_tddft and re.search(r"ORCA TD-DFT\/TDA CALCULATION", outfile_contents):
            continue

        # Find single point energy, skip if it happens not found:
        try:
            single_pt_Es[calc_name] = \
                re.search(r"FINAL SINGLE POINT ENERGY *(-?[\d]+[.][\d]+)\n", outfile_contents).group(1)
        except AttributeError:
            continue

    # Make Excel spreadsheet of the single point energies:
    single_pt_E_df = pd.Series(single_pt_Es).to_frame()
    single_pt_E_df = single_pt_E_df.reset_index()
    single_pt_E_df.to_excel(f"{destination_path}/single_pt_Es.xlsx", index=False)


if __name__ == "__main__":
    single_pt_E_workup(
        path_to_single_pt_out_files="data/bigjob/split1_result",
        destination_path="data/bigjob"
    )
