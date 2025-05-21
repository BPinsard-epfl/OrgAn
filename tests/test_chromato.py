from OrgAn import get_elution_order
from OrgAn import gives_data_frame

solutes = gives_data_frame("tests/test_smiles.csv")
print(get_elution_order(solutes))