"""Core functions."""

import networkx as nx
from pydantic import BaseModel
from pydantic._internal._model_construction import ModelMetaclass
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import Mol, RWMol
from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.rdmolops import AddHs

from . import rd


class _CustomBaseModelMeta(ModelMetaclass):
    def __getattr__(self, item: str):  # noqa: ANN204
        try:
            super().__getattr__(item)  # ty:ignore[unresolved-attribute]
        except AttributeError:
            if item in self.__dict__.get("__pydantic_fields__", ()):
                return item
            raise


class CustomBaseModel(BaseModel, metaclass=_CustomBaseModelMeta):
    """A custom base model that allows accessing field names as class attributes."""


class Atom(CustomBaseModel):
    """Represents an atom in a molecule."""

    symbol: str


class Bond(CustomBaseModel):
    """Represents a bond between two atoms in a molecule."""

    order: int


def validate(graph: nx.Graph) -> nx.Graph:
    """Validate the graph structure."""
    for key, data in graph.nodes(data=True):
        if not Atom.model_validate(data):
            msg = f"Node {key} does not have a valid Atom instance."
            raise ValueError(msg)
    for key1, key2, data in graph.edges(data=True):
        if not Bond.model_validate(data):
            msg = f"Edge ({key1}, {key2}) does not have a valid Bond instance."
            raise ValueError(msg)

    return graph


def from_smiles(smiles: str) -> nx.Graph:
    """Convert a SMILES string to a graph."""
    mol = MolFromSmiles(smiles)
    mol = AddHs(mol)
    graph = nx.Graph()

    for mol_atom in mol.GetAtoms():
        atom = Atom(symbol=mol_atom.GetSymbol())
        graph.add_node(mol_atom.GetIdx(), **atom.model_dump())

    for mol_bond in mol.GetBonds():
        bond = Bond(order=mol_bond.GetBondTypeAsDouble())
        graph.add_edge(
            mol_bond.GetBeginAtomIdx(), mol_bond.GetEndAtomIdx(), **bond.model_dump()
        )

    return validate(graph)


def rdkit_mol_with_index_map(graph: nx.Graph) -> tuple[Mol, dict[int, int]]:
    """Convert a graph back to an RDKit molecule."""
    rw_mol = RWMol()
    key_map: dict[int, int] = {}

    for key in sorted(graph.nodes()):
        atom = Atom(**graph.nodes[key])
        idx = rw_mol.AddAtom(rdchem.Atom(atom.symbol))
        key_map[key] = idx

    for key1, key2 in graph.edges():
        bond = Bond(**graph.edges[key1, key2])
        rw_mol.AddBond(key_map[key1], key_map[key2], rdchem.BondType(bond.order))

    index_map = dict(map(reversed, key_map.items()))
    return rw_mol.GetMol(), index_map


def rdkit_mol(graph: nx.Graph, *, label: bool = False) -> Mol:
    """Convert a graph back to an RDKit molecule."""
    mol, index_map = rdkit_mol_with_index_map(graph)
    if label:
        mol = rd.mol.with_atom_numbers(mol, index_map=index_map)
    return mol
