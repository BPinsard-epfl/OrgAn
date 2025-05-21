import pandas as pd
from OrgAn import gives_data_frame

def test_gives_dataframe():
    assert gives_data_frame("test_smiles.csv") == pd.read_csv("test_data.csv"), "givesDataFrame do not give the correct answer"

test_gives_dataframe()