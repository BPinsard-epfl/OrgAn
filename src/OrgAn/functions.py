import pandas as pd
import math as m

from rdkit import Chem
from pathlib import Path

from OrgAn import get_mol_info_from_smiles


datapath = "data/data.csv"
data : pd.DataFrame = pd.read_csv(datapath)



def gives_data_frame(file : str) -> pd.DataFrame: # TODO: give a more "conventional" name
    """
    Take a list of smile from a file and return a dataframe 
    with proprieties of each molecule according to PubChem.

    The file should a CSV table with one column of smiles.

    The Given data frames have for entries :
    name, canonical smile, CAS, formula, mol weight,
    pKa, logP, formal charge and sterimol.

    Sterimol is calculated using morfeus.
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
        "pKa" : [],
        "charge" : [],
        "sterimol_L": [],
        "sterimol_B1": [],
        "sterimol_B5": []
    }

    for i in range(len(df_smiles)):
        if Chem.CanonSmiles(df_smiles.iloc[i,0]) in data["smiles"].values:
            data_mol = data.loc[data["smiles"]==Chem.CanonSmiles(df_smiles.iloc[i,0])]
            for key in list(data.columns.values):
                dic_for_df[key].append(data_mol[key].tolist()[0])
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

    sorted_df = df.sort_values(by='pKa').reset_index(drop=True)
    maxs : list[tuple[float, tuple[int, int]]] = [(-1,(0,0)) for j in range(nb)]

    for i in range(len(sorted_df["pKa"])-1):
        for j in range(len(maxs)):
            if maxs[j][0] < round(float(sorted_df["pKa"].iloc[i+1] - sorted_df["pKa"].iloc[i]), 3):
                if j+1 < nb: maxs[j+1:len(maxs)] = maxs[j:len(maxs)-1]
                maxs[j] = (round(float(sorted_df["pKa"].iloc[i+1] - sorted_df["pKa"].iloc[i]), 3), (i, i+1))
                break
    
    maxs_dic : list[tuple[float, tuple[int, int]]] = []
    for i in range(len(maxs)):
        if maxs[i][1] == (0, 0): break
        maxs_dic.append(maxs[i])
    
    return maxs_dic


def find_logp_gaps(df : pd.DataFrame, nb : int = 1) -> dict[float, tuple[int, int]]: # does not handle equal values very well
    """
    A function used to find the biggest logP gaps of a dataframe. 
    It will give by default the biggest gap, but one can precise how many he wants.
    """

    sorted_df = df.sort_values(by='logP').reset_index(drop=True)
    maxs : list[tuple[float, tuple[int, int]]] = [(-1,(0,0)) for x in range(nb)]

    for i in range(len(sorted_df["logP"])-1):
        for j in range(len(maxs)):
            if maxs[j][0] <= round(float(sorted_df["logP"].iloc[i+1] - sorted_df["logP"].iloc[i]), 3):
                if j+1 < nb: maxs[j+1:len(maxs)] = maxs[j:len(maxs)-1]
                maxs[j] = (round(float(sorted_df["logP"].iloc[i+1] - sorted_df["logP"].iloc[i]), 3), (i, i+1))
                break
    
    maxs_dic : list[tuple[float, tuple[int, int]]] = []
    for i in range(len(maxs)):
        if maxs[i][1] == (0, 0): break
        maxs_dic.append(maxs[i])
    
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
                "name" : data_sorted["name"].iloc[0],
                "cid" : data_sorted["cid"].iloc[0],
                "CAS" : data_sorted["CAS"].iloc[0],
                "smiles" : data_sorted["smiles"].iloc[0],
                "molWeight" : data_sorted["molWeight"].iloc[0],
                "molFormula" : data_sorted["molFormula"].iloc[0],
                "logP" : data_sorted["logP"].iloc[0],
                "pKa" : data_sorted["pKa"].iloc[0],
                "charge" : data_sorted["charge"].iloc[0],
                "sterimol_L": data_sorted["sterimol_L"].iloc[0],
                "sterimol_B1": data_sorted["sterimol_B1"].iloc[0],
                "sterimol_B5": data_sorted["sterimol_B5"].iloc[0]
            }
            return pd.DataFrame([properties])
        else: 
            return pd.DataFrame([get_mol_info_from_smiles(smiles)])

    
    if charge != 100:
        data_sorted = data_sorted[data_sorted["charge"] == charge]

    #Poly values return

    if sterimol:  
        data_sorted.query("`sterimol_L`<=@sterimol['L']+3 and `sterimol_L`>=@sterimol['L']-3", inplace=True)
        data_sorted.query("`sterimol_B1`<=@sterimol['B1']+3 and `sterimol_B1`>=@sterimol['B1']-3", inplace=True)
        data_sorted.query("`sterimol_B5`<=@sterimol['B5']+3 and `sterimol_B5`>=@sterimol['B5']-3", inplace=True)
        data_sorted.sort_values(by=["sterimol_L", "sterimol_B1", "sterimol_B5"], key = lambda col : abs(col-sterimol[col.name.split("_", 1)[1]]), inplace=True)

    if logP != m.inf:
        data_sorted.query("`logP`<=@logP+1.5 and `logP`>=@logP-1.5", inplace=True)
        data_sorted.sort_values(by="logP", key = lambda col : abs(col-logP), inplace=True)

    if pka != m.inf:
        data_sorted.sort_values(by="pKa", key = lambda col : abs(col-pka), inplace=True)
    
    if data_sorted.dropna().empty:
        print("No compound found")    
    
    return data_sorted.head()
    
def test(col):
    print(col)
    return col