from OrgAn.chromato import get_elution_order
from OrgAn.functions import gives_data_frame

solutes = gives_data_frame("tests/test_smiles.csv")
print(get_elution_order(solutes))