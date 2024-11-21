import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs, UFFOptimizeMoleculeConfs


def lowest_conf_search_from_smiles(smiles: str, label: str, destination_path: str, force_field: str="MMFF",
                                   num_confs: int=100, use_seed: bool=False) -> None:
    """Performs a search for the lowest energy conformer of the mol from input SMILES string, saves as xyz file."""
    m = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(m)  # Necessary for reasonable conformers.
    if use_seed:
        rdDistGeom.EmbedMultipleConfs(m, numConfs=num_confs, randomSeed=0xf00d)
    else:
        rdDistGeom.EmbedMultipleConfs(m, numConfs=num_confs)

    if force_field == "MMFF":
        results = MMFFOptimizeMoleculeConfs(m)
    elif force_field == "UFF":
        results = UFFOptimizeMoleculeConfs(m)
    else:
        raise Exception("Unrecognized force_field.")

    # Find the lowest energy conformer:
    min_energy = results[0][1]
    min_energy_index = 0
    for i, conf in enumerate(results):
        if results[i][0] == 0 and results[i][1] < min_energy:
            min_energy_index = i
            min_energy = results[i][1]

    Chem.MolToXYZFile(mol=m, filename=f"{destination_path}/{label}.xyz", confId=min_energy_index)


def conformer_search(mols_filename: str, smiles_col: str, label_col: str,
                     destination_path: str, force_field: str) -> None:
    """Performs conformer search for phosphines from a DataFrame from a PubChem request."""
    for index, row in pd.read_csv(mols_filename).iterrows():
        print(f"Conformer search for label {row[f"{label_col}"]}: ", end="")
        lowest_conf_search_from_smiles(smiles=row[f"{smiles_col}"], label=row[f"{label_col}"],
                                       destination_path=destination_path, force_field=force_field)
        print("complete.")


if __name__ == "__main__":
    conformer_search(mols_filename="data/complex_smiles/pdmecl_set.csv",
                     smiles_col="smiles",
                     label_col="molecule_id",
                     destination_path="data/conf_search/pdmecl",
                     force_field="MMFF")
