import requests
import json
import regex as re

# r = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/165/JSON/")
# print(r.text)

def getPugRestJSONfromSmiles(smiles):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/JSON"
    req = requests.get(url)
    return req.text

def getPkaFromSmiles(smiles):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/cids/TXT"
    cid = int(requests.get(url).text)
    url1 = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/" + str(cid) + "/JSON/?heading=Dissociation+Constants"
    req = requests.get(url1)
    x = json.loads(req.text)
    try:
        pka = x['Record']['Section'][0]['Section'][0]['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
        return float(re.search(r'\d+\.\d+', pka).group())
    except:
        return "No pKa was found"
    

print(getPkaFromSmiles("CCO"))