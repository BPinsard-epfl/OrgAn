import pandas as pd
import math as m

from rdkit import Chem
from pathlib import Path

from pchem_rq import get_mol_info_from_smiles


datapath = "data/data.csv"
data : pd.DataFrame = pd.read_csv(datapath)



def gives_dataframe(file : str) -> pd.DataFrame:
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

    p = Path(".")
    file_path = p / file

    try:
        file_path.exists()
        file_path.is_file()
    except FileNotFoundError:
        assert(f"file `{file}` does not exist, or, perhaps, is not a file.") 

    try:
        df_smiles = pd.read_csv(file) 
    except AssertionError:
        assert("The file does not correspond to a CSV format.")

    dic_for_df : dict[str, list[str|float|int]]= {
        "name" : [],
        "cid" : [],
        "CAS" : [],
        "smiles" : [],
        "molWeight" : [],
        "molFormula" : [],
        "logP" : [],
        "is_pKa_parent_compound" : [],
        "pKa" : [],
        "charge" : [],
        "sterimol" : []
    }

    for i in range(len(df_smiles)):
        if Chem.CanonSmiles(df_smiles.iloc[i,0]) in data["smiles"]:
            data_mol = data["smiles"==Chem.CanonSmiles(df_smiles.iloc[i,0])].iloc[0]
            for key in list(data.columns.values):
                if key not in ["sterimol_L", "sterimol_B1", "sterimol_B5"]:
                    dic_for_df[key].append(data_mol[key])
            dic_for_df["is_pKa_parent_compound"].append(False)

        else:
            props = get_mol_info_from_smiles(df_smiles.iloc[i,0])
            for key in list(dic_for_df.keys()):
                dic_for_df[key].append(props[key])
    
    return pd.DataFrame(dic_for_df)

            

def find_pKa_gaps(df : pd.DataFrame, nb : int = 1) -> dict[float, tuple[int, int]]:
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


def find_logp_gaps(df : pd.DataFrame, nb : int = 1) -> dict[float, tuple[int, int]]:
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
 

def find_compounds(pka : float = m.inf, logP : float = m.inf, charge : int = 100,
                   sterimol : dict[str, float] = {}, smiles : str = "") -> pd.DataFrame:
    """
    Gives a specific acid based on Smiles, pka, Nucleophilicity, ...
    It will find in our database the best match, but if one gives a specific smiles in enter,
    it can search on Pubchem if not found in the database.
    """
    
    data_sorted = data.sort_values(by="pKa")
    
    #1 value return
    if smiles:
        try:
            smiles = Chem.CanonSmiles(smiles)
            data_sorted = data_sorted[data_sorted["smiles"] == smiles]
        except Exception as e :
            assert ValueError(f"The smiles is not correct. Error : \n{e}")

        if not data_sorted.empty:
            properties : dict[str, str|float] = {
                "Names" : data_sorted["Names"].iloc[0],
                "Smiles" : data_sorted["Smiles"].iloc[0],
                "pKa" : data_sorted["pka"].iloc[0],
                "Nucleophilicity" : data_sorted["Nucleophilicity"].iloc[0]
            }
            return pd.DataFrame(properties)
        else: 
            return pd.DataFrame(get_mol_info_from_smiles(smiles))

    
    if charge != 100:
        data_sorted = data_sorted[data_sorted["charge"] == charge]

    #Poly values return

    if sterimol:  
        data_sorted.query("`sterimol_L`<=sterimol["L"]+3 and `sterimol_L`>=sterimol["L"]-3")
        data_sorted.query("`sterimol_B1`<=sterimol["B1"]+3 and `sterimol_B1`>=sterimol["B1"]-3")
        data_sorted.query("`sterimol_B5`<=sterimol["B5"]+3 and `sterimol_B5`>=sterimol["B5"]-3")
        data_sorted.sort_values(by=["sterimol_L", "sterimol_B1", "sterimol_B5"], key = lambda col : abs(col-sterimol[col.index.split("_", 1)[1]]), inplace=True)

    if logP != m.inf:
        data_sorted.query("`logP`<=logP+1.5 and `logP`>=logP-1.5")
        data_sorted.sort_values(by="logP", key = lambda col : abs(col-logP), inplace=True)

    if pka != m.inf:
        data_sorted.sort_values(by="pKa", key = lambda col : abs(col-pka), inplace=True)
    
    if data_sorted.dropna().empty:
        print("No compound found")    
    
    return data_sorted.head()
    
    