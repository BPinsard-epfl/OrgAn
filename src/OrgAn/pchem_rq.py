import requests
import json
import regex as re
import urllib.parse
from morfeus import Sterimol
from rdkit import Chem
    
def get_mol_info_from_smiles(smiles: str) -> dict:
    """
    This function takes a SMILES string as an input, and then does a series of
    requests to PubChem to get various properties of the molecule. The properties
    are then outputted 
    """
    if type(smiles) != str:
        raise TypeError("Expected str, got " + str(type(smiles)))
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON?smiles=" + urllib.parse.quote_plus(smiles))
    try:
        mol = json.loads(req.text)['PC_Compounds'][0]
    except:
        mol = json.loads(req.text)['Fault']
        raise ValueError(mol['Code'] + "\n" + mol['Message'])
    molProperties = {
        "name": None, 
        "cid": mol['id']['id']['cid'], 
        "CAS": None,
        "smiles": None, 
        "molWeight": None,
        "molFormula": None,
        "logP": None,
        "is_pKa_parent_compound": False, # TODO: implement this
        "pKa": None,
        "charge": mol["charge"],
        "sterimol_L": None,
        "sterimol_B1": None,
        "sterimol_B5": None
        }
    props = mol['props']

    i = 0
    while True:
        try:
            currentProperty = props[i]["urn"]
            val = props[i]["value"]
        except:
            break
        i += 1
        try:
            if currentProperty["label"] not in ["IUPAC Name", "Log P", "Molecular Weight", "Molecular Formula", "SMILES"]: 
                continue
            elif currentProperty["label"] == "Molecular Formula":
                molProperties["molFormula"] = val["sval"]
            elif currentProperty["label"] == "Molecular Weight":
                molProperties["molWeight"] = float(val["sval"])
            elif currentProperty["label"] == "Log P":
                molProperties["logP"] = val["fval"]
            elif currentProperty["name"] == "Preferred":
                molProperties["name"] = val["sval"]
            elif currentProperty["name"] == "Canonical":
                molProperties["smiles"] = str(Chem.CanonSmiles(val["sval"]))
        except:
            continue
    
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(mol['id']['id']['cid']) + "/JSON/?heading=CAS")
    try:
        molProperties["CAS"] = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
    except:
        molProperties["CAS"] = None
    
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(mol['id']['id']['cid']) + "/JSON/?heading=Dissociation+Constants")
    try:
        j = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        molProperties["pKa"] = float(re.search(r'\d+\.\d+', j).group())
    except:
        molProperties["pKa"] = None
    try:
        req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(mol['id']['id']['cid']) + "/JSON/?record_type=3d")
        stmol = json.loads(req.text)['PC_Compounds'][0]
    except:
        pass
    else:
        elements = stmol['atoms']['element']
        y = stmol['coords'][0]['conformers'][0]['y']
        x = stmol['coords'][0]['conformers'][0]['x']
        z = stmol['coords'][0]['conformers'][0]['z']
        coords = []
        for i in range(len(x)):
            coords.append([x[i], y[i], z[i]])
        sterimol = Sterimol(elements, coords, stmol['bonds']['aid1'][0], stmol['bonds']['aid2'][0])
        molProperties["sterimol_L"] = round(float(sterimol.L_value), 2)
        molProperties["sterimol_B1"] = round(float(sterimol.B_1_value), 2)
        molProperties["sterimol_B5"] = round(float(sterimol.B_5_value), 2)
    
    return molProperties