import requests
import json
import regex as re

""" def getPugRestJsonFromSmiles(smiles):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/JSON"
    req = requests.get(url)
    return req.text """

""" def getPkaFromSmiles(smiles):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/cids/TXT"
    cid = int(requests.get(url).text)
    url1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(cid) + "/JSON/?heading=Dissociation+Constants"
    req = requests.get(url1)
    x = json.loads(req.text)
    try:
        pka = x['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        return float(re.search(r'\d+\.\d+', pka).group())
    except:
        return "No pKa was found" """
    
def getMoleculeInfoFromSmiles(smiles):
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/JSON")
    try:
        mol = json.loads(req.text)['PC_Compounds'][0]
    except:
        mol = json.loads(req.text)['Fault']
        raise Exception(mol['Code'] + "\n" + mol['Message'])
    props = mol['props']
    req = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(mol['id']['id']['cid']) + "/JSON/?heading=Dissociation+Constants")
    try:
        x = json.loads(req.text)['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        pka = float(re.search(r'\d+\.\d+', x).group())
    except:
        pka = None
    molProperties = {
        "name": props[6]['value']['sval'], 
        "cid": int(mol['id']['id']['cid']), 
        "smiles": props[20]['value']['sval'], 
        "logP": float(props[14]['value']['fval']), 
        "pKa": pka}
    return molProperties
    

print(getMoleculeInfoFromSmiles("O=C(O)c1c(C(O)=O)cccc1"))
print(getMoleculeInfoFromSmiles("CCO"))