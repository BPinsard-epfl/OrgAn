import requests
import json
import regex as re
    
def getMoleculeInfoFromSmiles(smiles):
    """
    This function takes a SMILES string as an input, and then does a series of
    requests to PubChem to get various properties of the molecule. The properties
    are then outputted 
    """
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/JSON")
    try:
        mol = json.loads(req.text)['PC_Compounds'][0]
    except:
        mol = json.loads(req.text)['Fault']
        raise Exception(mol['Code'] + "\n" + mol['Message'])
    molProperties = {
        "name": None, 
        "cid": mol['id']['id']['cid'], 
        "CASno": None,
        "smiles": None, 
        "molWeight": None,
        "molFormula": None,
        "logP": None,
        "is_pKa_parent_compound": False, 
        "pKa": None, # TODO: add pKa estimation boolean
        "charge": mol["charge"],
        "sterimol": None} # TODO: add sterimol thing using Morfeus – Sterimol
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
        x = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        molProperties["pKa"] = float(re.search(r'\d+\.\d+', x).group())
    except:
        molProperties["pKa"] = None
    return molProperties
    
# test functions
print(getMoleculeInfoFromSmiles("O=C(O)c1c(C(O)=O)cccc1"))
print(getMoleculeInfoFromSmiles("CCO"))