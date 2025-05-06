import pandas as pd
import math as m

from rdkit import Chem
from typing import Sequence
from pathlib import Path



datapath = "data/.csv"
data = pd.read_csv(datapath)


def SmilesDataFrame(file : str, sep_ : str = ";", skiprows_ : int | Sequence[int] | bool = False, skipblanklines_ : bool = True) -> pd.DataFrame:
    """
    This function will create the dataframe thanks to the file given.
    It will be then used in further functions.
    """
    
    p = Path(".")
    file_path = p / file

    try:
        file_path.exists()
        file_path.is_file()
    except FileExistsError:
        assert(f"file `{file}` does not exist, or, perhaps, is not a file.")

    try:
        df = pd.read_csv(file, sep = sep_, skiprows = skiprows_, skip_blank_lines = skipblanklines_ )
    except AssertionError:
        assert("The file does not correspond to a CSV format.")
    
    return df


def givesDataFrame(file : str, sep_ : str = ";", skiprows_ : int | Sequence[int] | bool = False, skipblanklines_ : bool = True) -> pd.DataFrame:
    """
    Take a list of smile from a file and return a dataframe 
    with proprieties of each molecule according to PubChem.

    The file should a CSV table with one column of smiles.

    The Given data frames have for entries :
    name, canonical smile, CAS, formula, mol weight,
    pKa, logP, formal charge and sterimol.

    When no pKa is found from pubchem it is predicted with 
    PypKa package.
    Sterimol is calculated thanks to wSterimol package.
    """

    df_smiles = df_smiles(file, sep_, skiprows_, skipblanklines_)
    df_smiles.iloc[0:0,:]= df_smiles.iloc[0:0,:].apply(Chem.MolToSmiles(Chem.MolFromSmiles()))

    dic_for_df = {
        "name" : [],
        "cid" : [],
        "smiles" : [],
        "pKa" : [],
        "logP" : []
    }

    for i in range(len(df_smiles)):
        props = getMoleculeInfoFromSmiles(df_smiles.iloc[0,i])
        for p in ["name", "cid", "smiles", "pKa", "logP"]:
            dic_for_df[p].append(props[p])
    
    return pd.DataFrame(dic_for_df)

            

def findpKaGaps(df : pd.DataFrame, nb : int = 1) -> dict[float, tuple[int, int]]:
    """
    A function used to find the biggest pka gaps of a dataframe. 
    It will give by default the biggest gap, but one can precise how many he wants.
    """

    sorted_df = df.sort_values(by='pKa')
    maxs : list[tuple[float, tuple[int, int]]] = [(0,(0,0)) for j in range(nb)]

    for i in range(len(sorted_df["pKa"])-1):
        for j in range(len(maxs)):
            if maxs[j] < sorted_df["pKa"].iloc(i) - sorted_df["pKa"].iloc(i+1):
                maxs[j] = (sorted_df["pKa"].iloc(i) - sorted_df["pKa"].iloc(i+1), (i, i+1))
                break
    
    maxs_dic : dict[float, tuple[int, int]]
    for i in range(len(maxs)):
        maxs_dic[maxs[i][0]] = maxs[i][1]
    
    return maxs_dic


def findLogPGaps(df : pd.DataFrame, nb : int = 1) -> dict[float, tuple[int, int]]:
    """
    A function used to find the biggest logP gaps of a dataframe. 
    It will give by default the biggest gap, but one can precise how many he wants.
    """

    sorted_df = df.sort_values(by='logP')
    maxs : list[tuple[float, tuple[int, int]]] = [(0,(0,0)) for x in range(nb)]

    for i in range(len(sorted_df["logP"])-1):
        for j in range(len(maxs)):
            if maxs[j] < (sorted_df["logP"].iloc(i) - sorted_df["logP"].iloc(i+1)):
                maxs[j] = (sorted_df["logP"].iloc(i) - sorted_df["logP"].iloc(i+1), (i, i+1))
                break
    
    maxs_dic : dict[float, tuple[int, int]]
    for i in range(len(maxs)):
        maxs_dic[maxs[i][0]] = maxs[i][1]
    
    return maxs_dic
 

def findAcidOrBase(pka : float = m.inf, smiles : str = "", logP : float = m.inf, nu : str = "") -> pd.DataFrame:
    """
    Gives a specific acid based on Smiles, pka, Nucleophilicity, ...
    It will find in our database the best match, but if one gives a specific smiles in enter,
    it can search on Pubchem if not found in the database.
    """
    
    data_sorted = data.sort_values(by="pKa")
    
    if smiles:
        try:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            data_sorted = data_sorted[data_sorted["Smiles"] == smiles]
        except Exception as e :
            assert ValueError(f"The smiles is not correct. Error : \n{e}")

        if not data_sorted.empty:
            properties : dict[str, str|float] = {
                "Names" : data_sorted["Names"].iloc[0],
                "Smiles" : data_sorted["Smiles"].iloc[0],
                "pKa" : data_sorted["pka"].iloc[0],
                "Nucleophilicity" : data_sorted["Nucleophilicity"].iloc[0]
            }
            return properties
        elif ...: # faire avec un call pubchem
            ...
        else: 
            assert KeyError("Could not find your molecule.")
    
    if nu:
        data_sorted = data_sorted[data_sorted["Nucleophilicity"] == nu]

    if logP != m.inf:
        data_sorted = data_sorted[data_sorted["LogP"] == logP]

    if ...: #mettre tous les filtres 
        ...

    if pka == m.inf:
        filter : dict[str, float|int]  = {"pKa diff" : m.inf, "index" : int(0)}
        for i in range(len(data_sorted)):
            if abs(data_sorted["pKa"].iloc[i]-pka) < filter["pKa diff"]:
                filter["pKa diff"] = abs(data_sorted["pKa"].iloc[i]-pka)
                filter["index"] = i
            
        result = data_sorted[data_sorted.iloc[filter["index"]-2:filter["index"]+2]]
        return result
    
    else: 
        return data_sorted.head()
    