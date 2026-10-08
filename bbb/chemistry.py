"""Compute the molecular graph, descriptors, and 2D coordinates from SMILES."""
from functools import lru_cache
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors, Lipinski
from rdkit.Chem.Draw import rdMolDraw2D


@lru_cache(maxsize=256)
def molecule_data(smiles, draw=True):
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > 16000:
        raise ValueError("SMILES required (up to 16,000 characters).")
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("Invalid SMILES or atom valence.")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Multiple fragments found. Supply one molecule.")
    if mol.GetNumAtoms() > 500:
        raise ValueError("Maximum 500 atoms per molecule.")
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Undefined atoms (*) are not supported.")
    result = {
        "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "descriptors": {
            "AroR": rdMolDescriptors.CalcNumAromaticRings(mol),
            "HA": mol.GetNumHeavyAtoms(), "MW": Descriptors.MolWt(mol),
            "HBA": Lipinski.NumHAcceptors(mol), "HBD": Lipinski.NumHDonors(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
        },
        "rdkit_version": rdBase.rdkitVersion,
        "svg": None, "molblock": None,
    }
    if draw:
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(640, 400)
        drawer.drawOptions().padding = .12
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        result["svg"] = drawer.GetDrawingText()
        result["molblock"] = Chem.MolToMolBlock(mol)
    return result
