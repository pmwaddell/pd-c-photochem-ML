from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd


def draw_from_phos_df(phos_df: pd.DataFrame, filename: str, legend: str="molecule_id") -> None:
    """Draws the set of phosphines from a phosphine DataFrame."""
    if legend:
        Draw.MolsToGridImage(
            [Chem.MolFromSmiles(smiles) for smiles in phos_df["smiles"]],
            molsPerRow=10,
            subImgSize=(400,400),
            legends=phos_df[legend].astype(str).to_list()
        ).save(filename)
    else:
        Draw.MolsToGridImage(
            [Chem.MolFromSmiles(smiles) for smiles in phos_df["smiles"]],
            molsPerRow=10,
            subImgSize=(400,400)
        ).save(filename)


p_df = pd.read_csv("data/kraken_data/ml_8_210.csv")
# print(p_df)

# We want to ensure that the phosphines would form T-shaped complexes with Pd, to the best of our ability.
# Thus, we will look for phosphines that have steric properties which indicate they are at least as large as (tBu)3P.
# In particular, I have decided to use the cone angle of the (tBu)3P-Ni(CO)3 complex and %Vbur as the limiting criteria.
# I think using the Boltzmann averages is most representative, since these shoudld control the overall thermodynamic
# viability of forming a T-shaped complex.
limiting_cone_angle = p_df[p_df.molecule_id == 8]["vbur_vbur_boltzmann_average"]


pruned = p_df[p_df.vbur_vbur_boltzmann_average >= 68]
# print(pruned)
# draw_from_phos_df(pruned)
pruned.to_csv('data/kraken_data/pruned_phos_set.csv', index=False)

potential_problems = [
    532, 835, 1136, 19872, 19874, 19973, 19981, 19991, 19998, 20061, 20101, 20065, 20136, 20137, 20140,
    20230, 20220, 20185, 20154, 20153, 20148, 20147, 20305, 20320, 20329, 20332, 192978, 183039, 183168,
    247447, 247699, 247730
]
# Maybe we want to leave them in at least at first, see if they lead to high HOMO/LUMO gaps,
# LUMOs of the wrong character, etc.