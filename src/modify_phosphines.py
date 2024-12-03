from typing import Callable

from rdkit import Chem
import pandas as pd


def smiles_to_canonical_smiles(smiles) -> str:
    """Converts any given SMILES representation of a molecule to the canonical SMILES representation."""
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles)
    )


def modify_phos_smiles(phos_smiles, mod: str) -> str:
    """Modifies a phosphine and modifies it according to the input string."""
    for i, c in enumerate(phos_smiles):
        if c == "P":
            phos_smiles = phos_smiles[:i] + mod + phos_smiles[i+1:]
            break
    return smiles_to_canonical_smiles(phos_smiles)


def phos_smiles_to_pdmecl_complex_smiles(phos_smiles: str) -> str:
    """Takes SMILES string for phosphine and returns SMILES string for corresponding R3P-Pd(Me)Cl complex."""
    return modify_phos_smiles(phos_smiles, "[P+]([Pd](C)Cl)")


def phos_smiles_to_nico3_complex_smiles(phos_smiles: str) -> str:
    """Takes SMILES string for phosphine and returns SMILES string for corresponding R3P-Ni(CO)3 complex."""
    return modify_phos_smiles(phos_smiles, "[P+]([Ni]([C-]#[O+])([C-]#[O+])([C-]#[O+]))")


def phos_smiles_to_phos_oxide_smiles(phos_smiles: str) -> str:
    """Takes SMILES string for phosphine and returns SMILES string for corresponding phosphine oxide."""
    return modify_phos_smiles(phos_smiles, "P(=O)")


def make_modified_phos_df(phos_df: pd.DataFrame, mod: Callable, label_col: str, label_ending: str) -> pd.DataFrame:
    """Returns a DataFrame of modified phosphines from a phosphine dataframe."""
    result = pd.DataFrame({"smiles": [mod(smiles) for smiles in phos_df.smiles]})
    # Have to convert the labels column to strings in case they are not strings:
    result.insert(0, label_col, phos_df[label_col].astype('str') + label_ending)
    return result


if __name__ == "__main__":
    from prune_phosphines import draw_from_phos_df
    df = pd.read_csv('data/kraken_data/pruned_phos_set.csv').iloc[:, :2]

    complex_df = make_modified_phos_df(
        df, phos_smiles_to_pdmecl_complex_smiles, "molecule_id", "_PdMeCl")
    complex_df.to_csv('data/complex_smiles/PdMeCl_set.csv', index=False)
    # draw_from_phos_df(complex_df, '../misc/extra_complex_smiles/PdMeCl_set.png')
