import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.SaltRemover import SaltRemover

def neutralize_standardize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]

    try:
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        return Chem.MolToSmiles(mol) # Isso aqui já canoniza os SMILES?
    except:
        return 'cant be neutralized'


class SmilesCleaner:
    """
    A class for cleaning and processing SMILES data in a DataFrame.
    """

    def __init__(self, dataframe:pd.DataFrame):
        """
        Constructor for SmilesCleaner.

        Parameters:
            dataframe (pandas.DataFrame): The input DataFrame containing SMILES data.
        """
        self.df =  dataframe

    def strip_salt(self, smiles_col, output_col):
        """
        Strip salts from the molecules in the specified column.

        Parameters:
            smiles_col (str): Name of the column containing SMILES strings.
            output_col (str): Name of the column to store the stripped SMILES.

        Returns:
            SmilesCleaner: An instance of SmilesCleaner with the updated DataFrame.
        """
        self.remover = SaltRemover()
        # Create mol from raw_mol column and then strip its salt
        PandasTools.AddMoleculeColumnToFrame(smilesCol=smiles_col, molCol='raw_mol', frame=self.df)
        self.df['stripped_salt_mol'] = self.df['raw_mol'].apply(lambda x: self.remover.StripMol(x, dontRemoveEverything=True))

        # Create stripped salt smiles and drop mol column
        self.df[output_col] = self.df['stripped_salt_mol'].apply(Chem.MolToSmiles)
        self.df.drop(columns=['stripped_salt_mol', 'raw_mol'], axis=1, inplace=True) # drop stripped salt mol
        
        return self

        
    def neutralize(self, smiles_col, output_col):
        """
        Neutralize the molecules in the specified column.

        Parameters:
            smiles_col (str): Name of the column containing SMILES strings.
            output_col (str): Name of the column to store the neutralized SMILES.

        Returns:
            SmilesCleaner: An instance of SmilesCleaner with the updated DataFrame.
        """
        # Recreate the mol using the stripped salt smiles
        PandasTools.AddMoleculeColumnToFrame(smilesCol=smiles_col, molCol='stripped_salt_mol', frame=self.df)
        
        # Create the neutralized smiles column and return the its smiles
        self.df[output_col] = self.df['stripped_salt_mol'].apply(neutralize_standardize_atoms)
        self.df.drop(columns=['stripped_salt_mol'], axis=1, inplace=True)
        
        return self

    
    def search_duplicate(self, smiles_col, keep_inchi=False):
        # Dependendo da funcao, ja dropar
        """
        Search duplicate molecules based on their SMILES or InChI.

        Parameters:
            smiles_col (str): Name of the column containing SMILES strings.
            keep_inchi (bool): Whether to keep the InChI column (default is False).

        Returns:
            SmilesCleaner: An instance of SmilesCleaner with the updated DataFrame.
        """
        PandasTools.AddMoleculeColumnToFrame(smilesCol=smiles_col, molCol='clean_mol', frame=self.df)
        self.df['inchi'] = self.df['clean_mol'].apply(Chem.MolToInchi)

        # now dropping the duplicates and signaling only duplicated
        self.df.drop(columns=['clean_mol'], axis=1, inplace=True)
        self.df['duplicated'] = self.df.duplicated(subset=['inchi'], keep=False)

        if not keep_inchi:
            self.df.drop(columns=['inchi'], axis=1, inplace=True)
            return self
        else:
            return self

