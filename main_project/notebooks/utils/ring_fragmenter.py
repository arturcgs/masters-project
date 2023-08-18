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


def get_ring_fragments(list_of_smiles: list, no_rings_list=False) -> pd.DataFrame:
    """
    Recieves a list of SMILES and returnas a pd.DataFrame of parent SMILES and the respective Ring Fragment (Ring and Adjacent Atoms).
    It depends on: get_ring_systems() and get_ring_adjacent() and modules pandas and rdkit;

    Parameters:
    list_of_smiles: a list of SMILES to be extracte the ring fragments.
    no_rings_list: can also return the list of SMILES that have no ring on the list.
    """
    # empty lists
    parents_and_fragments = []
    no_rings = []

    for smiles in list_of_smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except Exception as molConversionError:
            print(f'It was not possible to convert the structure {smiles} into mol object: {molConversionError}')
            continue
        try:
            # If have no rings the get_ring_systems function returns [0]
            if get_ring_systems(mol) == [0]:
                no_rings.append(smiles)
            else:
                rings = get_ring_adjacent(mol)
                for ring in rings:
                    ring_fragment = Chem.MolToSmiles(
                                        Chem.MolFromSmiles(
                                            Chem.MolFragmentToSmiles(mol, atomsToUse=ring)
                                        )
                                    )
                    parent_smiles = smiles
                    parents_and_fragments.append(
                        {
                            'parent_smiles': parent_smiles,
                            'ring_fragment': ring_fragment
                        }
                    )
        except Exception as ringRetrievalError:
            print(f'It was not possible to get the ring structures from {smiles}: {ringRetrievalError}')
    
    if no_rings_list:
        return pd.DataFrame.from_dict(parent_smiles), no_rings
    else:
        return pd.DataFrame.from_dict(parents_and_fragments)