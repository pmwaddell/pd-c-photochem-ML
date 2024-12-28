import glob
import re
from typing import List

import pandas as pd


def single_pt_E_workup(path_to_single_pt_out_files: str, destination_path: str, skip_tddft: bool=False,
                       suffixes_to_skip: List=("atom46")) -> None:
    """
    Finds all .out files recursively in a given directory and checks if they have single point calculations.
    For the ones that do, their single point energies are compiled into an Excel spreadsheet.
    """
    single_pt_out_file_paths = glob.glob(f"{path_to_single_pt_out_files}/**/*.out", recursive=True)
    single_pt_Es = {}

    for single_pt_out_file in single_pt_out_file_paths:
        single_pt_out_file_path = single_pt_out_file.replace("\\", "/")

        calc_name = single_pt_out_file_path.split("/")[-1][:-4]  # remove .out

        skip_because_of_suffix = False
        for suffix in suffixes_to_skip:
            try:
                if calc_name[-len(suffix):] == suffix:
                    skip_because_of_suffix = True
                    break
            except IndexError:
                continue
        if skip_because_of_suffix:
            continue

        with open(single_pt_out_file_path, 'r') as out_file:
            outfile_contents = out_file.read()
        # Reverse outfile contents: we always want to find the LAST occurance of "FINAL SINGLE POINT ENERGY" in the file
        # since there can be more than one, e.g. for geometry optimizations.
        # Thus, if we reverse the file's contents line-by-line, the first match will always be the final energy.
        outfile_contents = '\n'.join(outfile_contents.split('\n')[::-1])

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
    single_pt_E_df = single_pt_E_df.sort_values(by="index")  # sort by molecule label
    single_pt_E_df.to_excel(f"{destination_path}/single_pt_Es.xlsx", index=False)


if __name__ == "__main__":
    single_pt_E_workup(
        path_to_single_pt_out_files="data/bigjob",
        destination_path="data/bigjob",
        suffixes_to_skip=["tddft", "atom46", "geom_opt"]
    )
