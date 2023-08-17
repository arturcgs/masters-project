import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFragmentToSmiles
from rdkit.Chem.PandasTools import AddMoleculeColumnToFrame

def get_ring_systems(mol, includeSpiro=False):
    """
    Retrieve atom indices grouped into ring systems within a given molecule.
    
    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object representing the molecular structure.
        includeSpiro (bool, optional): If True, include spiro-connected atoms in the ring systems. 
            Defaults to False.

    Returns:
        list of set: A list of sets, where each set contains the indices of atoms belonging to a ring system.
    """
    ri = mol.GetRingInfo()
    systems = []

    # Scaping
    if len(ri.AtomRings()) == 0:
        return [0] # Return a carbon
    
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))

            if nInCommon and (includeSpiro or nInCommon>1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems

def get_ring_adjacent(mol):
    """
    Retrieve atom indices within ring systems along with their adjacent atoms. This functions have the get_ring_systems() as a dependency;

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule object representing the molecular structure.

    Returns:
        list of set: A list of sets, where each set contains the indices of atoms in a ring system 
        along with their adjacent atoms.
    """
    ring_systems = get_ring_systems(mol, includeSpiro=False)

    rings = []
    for ring in ring_systems:
        ring_with_adjacent = set(ring)
        for ring_atom in ring:
            neighbors = set(mol.GetAtomWithIdx(ring_atom).GetNeighbors())
            
            for neighbor_atom in neighbors:
                if (neighbor_atom.GetIdx() not in ring) and (neighbor_atom.GetIsAromatic()):
                    neighbor_atom.SetIsAromatic(False)
                    ring_with_adjacent.add(neighbor_atom.GetIdx())
                else:
                    ring_with_adjacent.add(neighbor_atom.GetIdx())
        rings.append(ring_with_adjacent)

    return rings

