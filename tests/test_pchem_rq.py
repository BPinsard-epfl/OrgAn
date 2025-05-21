import os
import sys
current_dir = os.getcwd()
target_dir_relative = os.path.join("src", "OrgAn")
target_dir_absolute = os.path.abspath(os.path.join(current_dir, target_dir_relative))
sys.path.append(target_dir_absolute)

from pchem_rq import get_mol_info_from_smiles

butane = {
        "name": "butane", 
        "cid": 7843, 
        "CAS": "106-97-8",
        "smiles": "CCCC", 
        "molWeight": 58.12,
        "molFormula": "C4H10",
        "logP": 2.9,
        "is_pKa_parent_compound": False,
        "pKa": None,
        "charge": 0,
        "sterimol_L": 4.68,
        "sterimol_B1": 1.99,
        "sterimol_B5": 3.26
        }

benzoic_acid = {
        "name": "benzoic acid", 
        "cid": 243, 
        "CAS": "65-85-0",
        "smiles": "O=C(O)c1ccccc1", 
        "molWeight": 122.12,
        "molFormula": "C7H6O2",
        "logP": 1.9,
        "is_pKa_parent_compound": False,
        "pKa": 4.207,
        "charge": 0,
        "sterimol_L": 6.41,
        "sterimol_B1": 1.7,
        "sterimol_B5": 6.02
        }

tetraethylammonium = {
        "name": "tetraethylazanium", 
        "cid": 5413, 
        "CAS": "66-40-0",
        "smiles": "CC[N+](CC)(CC)CC", 
        "molWeight": 130.25,
        "molFormula": "C8H20N+",
        "logP": 1.7,
        "is_pKa_parent_compound": False,
        "pKa": None,
        "charge": 1,
        "sterimol_L": 4.62,
        "sterimol_B1": 3.13,
        "sterimol_B5": 4.5
        }

sodium_tert_butoxide = {
        "name": "sodium;2-methylpropan-2-olate", 
        "cid": 23676156, 
        "CAS": "865-48-5",
        "smiles": "CC(C)(C)[O-].[Na+]", 
        "molWeight": 96.1,
        "molFormula": "C4H9NaO",
        "logP": None,
        "is_pKa_parent_compound": False,
        "pKa": None,
        "charge": 0,
        "sterimol_L": None,
        "sterimol_B1": None,
        "sterimol_B5": None
        }

def test_bad_smiles(bad_smiles : str) -> bool:
    try:
        get_mol_info_from_smiles(bad_smiles)
    except ValueError:
        return True
    except:
        return False
    else:
        return False
    
def test_bad_type(bad_variable) -> bool:
    try:
        get_mol_info_from_smiles(bad_variable)
    except TypeError:
        return True
    except:
        return False
    else:
        return False
    
assert get_mol_info_from_smiles("CCCC") == butane, "The function did not yield the correct result"
assert get_mol_info_from_smiles("C1=CC=CC=C1C(=O)O") == benzoic_acid, "The function did not yield the correct result"
assert get_mol_info_from_smiles("[N+](CC)(CC)(CC)CC") == tetraethylammonium, "The function did not yield the correct result"
assert get_mol_info_from_smiles("CC(C)(C)[O-].[Na+]") == sodium_tert_butoxide, "The function did not yield the correct result"
assert test_bad_smiles("this is a bad smiles"), "The function did not raise the correct exception"
assert not test_bad_smiles(123456789), "The function did not raise the correct exception"
assert test_bad_type(123456789), "The function did not raise the correct exception"
assert test_bad_type(True), "The function did not raise the correct exception"
assert not test_bad_type("this is a bad smiles"), "The function did not raise the correct exception"
print("All tests passed successfully")