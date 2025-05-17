import requests
import json
import regex as re
import urllib.parse
from morfeus import Sterimol
    
def get_mol_info_from_smiles(smiles: str) -> dict:
    """
    This function takes a SMILES string as an input, and then does a series of
    requests to PubChem to get various properties of the molecule. The properties
    are then outputted 
    """
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON?smiles=" + urllib.parse.quote_plus(smiles))
    try:
        mol = json.loads(req.text)['PC_Compounds'][0]
    except:
        mol = json.loads(req.text)['Fault']
        raise Exception(mol['Code'] + "\n" + mol['Message'])
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
        "sterimol": {
            "L": None,
            "B_1": None,
            "B_5": None
        }}
    props = mol['props']

    i = 0
    while True: # This might be a bit of a barbaric approach but hey, if it works, it works. Basically I run the loop until the index goes out of range, which will raise an exception, indicating we've reached the end of the list.
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
                molProperties["smiles"] = val["sval"]
        except:
            continue
    
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(mol['id']['id']['cid']) + "/JSON/?heading=CAS")
    try:
        molProperties["CASno"] = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
    except:
        molProperties["CASno"] = None
    
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(mol['id']['id']['cid']) + "/JSON/?heading=Dissociation+Constants")
    try:
        j = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        molProperties["pKa"] = float(re.search(r'\d+\.\d+', j).group())
    except:
        molProperties["pKa"] = None
    
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(mol['id']['id']['cid']) + "/JSON/?record_type=3d")
    stmol = json.loads(req.text)['PC_Compounds'][0]
    elements = stmol['atoms']['element']
    y = stmol['coords'][0]['conformers'][0]['y']
    x = stmol['coords'][0]['conformers'][0]['x']
    z = stmol['coords'][0]['conformers'][0]['z']
    coords = []
    for i in range(len(x)):
        coords.append([x[i], y[i], z[i]])
    sterimol = Sterimol(elements, coords, stmol['bonds']['aid1'][0], stmol['bonds']['aid2'][0])
    molProperties['sterimol'] = {
        "L": round(float(sterimol.L_value), 2),
        "B_1": round(float(sterimol.B_1_value), 2),
        "B_5": round(float(sterimol.B_5_value), 2)
    }

    return molProperties
    
# test functions TODO: remove once finished
#print(getMoleculeInfoFromSmiles("O=C(O)c1c(C(O)=O)cccc1"))
#print(getMoleculeInfoFromSmiles("CCO"))