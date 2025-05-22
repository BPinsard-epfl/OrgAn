import pandas as pd
from OrgAn import gives_data_frame, find_pKa_gaps, find_logp_gaps, find_compounds

print(find_compounds(smiles="CCO"))

df = gives_data_frame("tests/test_smiles.csv")

assert df.equals(pd.read_csv("tests/test_data.csv")), "gives_data_frame did not parse the csv correctly"
print("gives_data_frame passed all tests successfully")

assert find_pKa_gaps(df) == [(11.026, (1, 2))], "find_pKa_gaps did not find the correct gaps"
assert find_pKa_gaps(df, 2) == [(11.026, (1, 2)), (4.1, (2, 3))], "find_pKa_gaps did not find the correct gaps"
assert find_pKa_gaps(df, 3) == [(11.026, (1, 2)), (4.1, (2, 3)), (0.118, (0, 1))], "find_pKa_gaps did not find the correct gaps"
assert find_pKa_gaps(df, 5) == [(11.026, (1, 2)), (4.1, (2, 3)), (0.118, (0, 1))], "find_pKa_gaps did not find the correct gaps"
assert find_pKa_gaps(df, 30) == [(11.026, (1, 2)), (4.1, (2, 3)), (0.118, (0, 1))], "find_pKa_gaps did not find the correct gaps"
print("find_pKa_gaps passed all tests successfully")

assert find_logp_gaps(df) == [(2.8, (0, 1))], "find_logp_gaps did not find the correct gaps"
assert find_logp_gaps(df, 2) == [(2.8, (0, 1)), (1.6, (6, 7))], "find_logp_gaps did not find the correct gaps"
assert find_logp_gaps(df, 3) == [(2.8, (0, 1)), (1.6, (6, 7)), (0.9, (5, 6))], "find_logp_gaps did not find the correct gaps"
assert find_logp_gaps(df, 7) == [(2.8, (0, 1)), (1.6, (6, 7)), (0.9, (5, 6)), (0.4, (3, 4)), (0.1, (4, 5)), (0.1, (1, 2)), (0.0, (2, 3))], "find_logp_gaps did not find the correct gaps"
assert find_logp_gaps(df, 30) == [(2.8, (0, 1)), (1.6, (6, 7)), (0.9, (5, 6)), (0.4, (3, 4)), (0.1, (4, 5)), (0.1, (1, 2)), (0.0, (2, 3))], "find_logp_gaps did not find the correct gaps"
print("find_logp_gaps passed all tests successfully")

assert find_compounds(smiles="CCCC").equals(pd.DataFrame([{"name" : "butane", "cid" : 7843, "CAS" : "106-97-8", "smiles" : "CCCC", "molWeight" : 58.12, "molFormula" : "C4H10", "logP" : 2.9, "is_pKa_parent_compound": False, "pKa" : None, "charge" : 0, "sterimol_L" : 4.68, "sterimol_B1" : 1.99, "sterimol_B5" : 3.26}])), "find_compounds failed to parse information from a pubchem request"
print("find_compounds passed all tests successfully")

print("all tests were successful")