"""RDKit Mol functions."""

from rdkit.Chem.rdchem import Mol


def with_atom_numbers(mol: Mol, index_map: dict[int, int]) -> Mol:
    """Return a molecule with subscript atom numbers."""
    for atom in mol.GetAtoms():
        index = index_map[atom.GetIdx()]
        symbol = atom.GetSymbol()
        atom.SetProp("atomLabel", f"{symbol}{index}")

    return mol
