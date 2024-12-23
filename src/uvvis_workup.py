import glob
import os

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
from rdkit.Chem import Draw


def make_uvvis_img_from_uvvis_df(uvvis_df: pd.DataFrame, molecule_id: str, destination_path: str,
                                 path_to_smiles_csv: str, complex_picture: bool=True) -> None:
    """Creates an image with two views of the TDDFT UV-vis spectrum and an image of the complex."""
    if complex_picture:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        fig.set_figheight(5)
        fig.set_figwidth(20)
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_figheight(5)
        fig.set_figwidth(14)

    ax1.plot(uvvis_df["wavelength"], uvvis_df["intensity"])
    ax1.set_xlabel("wavelength (nm)")
    ax1.set_ylabel("intensity")
    ax1.set_title(f"{molecule_id} UV-vis")

    ax2.plot(uvvis_df["wavelength"], uvvis_df["intensity"])
    ax2.set_xlim(300, 550)

    # Set y limits to be the maximum intensity along 340-550, scale intensity appropriately:
    zoom_subrange = uvvis_df.loc[(uvvis_df["wavelength"] < 550)]
    zoom_subrange = zoom_subrange.loc[(uvvis_df["wavelength"] > 340)]
    max_y = max(zoom_subrange["intensity"])
    min_y = min(zoom_subrange["intensity"])
    ax2.set_ylim(min_y, max_y + (max_y * 0.5))
    max_intensity = max(zoom_subrange["intensity"])
    lambda_max = zoom_subrange[zoom_subrange.intensity == max_intensity]["wavelength"].iloc[0]
    ax2.annotate(f"lambda max: {lambda_max}", xy=(lambda_max, max_intensity + 100))

    ax2.set_xlabel("wavelength (nm)")
    ax2.set_ylabel("intensity")
    ax2.set_title(f"{molecule_id} UV-vis, zoom")

    if complex_picture:
        # Add the image of the complex
        PdMeCl_set = pd.read_csv(path_to_smiles_csv)
        smiles = PdMeCl_set[PdMeCl_set.molecule_id == molecule_id].iloc[0, 1]
        Draw.MolToFile(Chem.MolFromSmiles(smiles), f"{molecule_id}.png")
        img = mpimg.imread(f"{molecule_id}.png")
        ax3.imshow(img)
        ax3.axis('off')
        ax3.set_title(molecule_id)
        os.remove(f"{molecule_id}.png")

    plt.savefig(f"{destination_path}/{molecule_id}_uvvis.png")
    plt.close()


def uvvis_workup(path_to_uvvis_files: str, destination_path: str, path_to_smiles_csv: str,
                 complex_pictures: bool=True) -> None:
    """
    Finds all .out.abs.dat files recursively in a given directory and produces images of their UV-vis spectra, and
    also compiles their spectral data and lambda max data into respective excel spreadsheets.
    """
    uvvis_file_paths = glob.glob(f"{path_to_uvvis_files}/**/*.out.abs.dat", recursive=True)
    uvvis_dfs = {}
    lambda_maxes = {}

    for uvvis_file_path in uvvis_file_paths:
        uvvis_file_path = uvvis_file_path.replace("\\", "/")
        molecule_id = uvvis_file_path.split("/")[-1][:-18]  # remove the _tddft.out.abs.dat

        uvvis_df = pd.read_csv(uvvis_file_path, sep='\\s+', header=None, usecols=[0, 1])
        uvvis_df = uvvis_df.rename(columns={0: "wavenumber", 1: "intensity"})
        uvvis_df.insert(1, "wavelength", round(10000000 / uvvis_df["wavenumber"], 1))
        uvvis_dfs[molecule_id] = uvvis_df

        make_uvvis_img_from_uvvis_df(uvvis_df, molecule_id=molecule_id, destination_path=destination_path,
                                     path_to_smiles_csv=path_to_smiles_csv, complex_picture=complex_pictures)

        # Find lambda max in the range 340-550 nm:
        zoom_subrange = uvvis_df.loc[(uvvis_df["wavelength"] < 550)]
        zoom_subrange = zoom_subrange.loc[(uvvis_df["wavelength"] > 340)]
        max_intensity = max(zoom_subrange["intensity"])
        lambda_max = zoom_subrange[zoom_subrange.intensity == max_intensity]["wavelength"].iloc[0]
        lambda_maxes[molecule_id] = lambda_max

    with pd.ExcelWriter(f'{destination_path}/uvvis_data.xlsx') as writer:
        for molecule_id in uvvis_dfs.keys():
            uvvis_dfs[molecule_id].to_excel(writer, sheet_name=molecule_id, index=False)

    # Make lambda max excel spreadsheet:
    lambda_max_df = pd.Series(lambda_maxes).to_frame()
    lambda_max_df = lambda_max_df.reset_index()
    lambda_max_df.to_excel(f"{destination_path}/lambda_maxes.xlsx", index=False)


if __name__ == "__main__":
    # uvvis_workup(
    #     path_to_uvvis_files="data/bigjob/split1_result",
    #     destination_path="data/bigjob/split1_result/uvvis_workup",
    #     path_to_smiles_csv="data/complex_smiles/PdMeCl_set.csv"
    # )

    uvvis_workup(
        path_to_uvvis_files="data/electronic_effects_study/results",
        destination_path="data/electronic_effects_study/results/uvvis_workup",
        path_to_smiles_csv="data/complex_smiles/PdMeCl_set.csv",
        complex_pictures=False
    )
