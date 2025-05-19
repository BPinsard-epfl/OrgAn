from functions import givesDataFrame
import pandas as pd

def test_gives_dataframe():
    result = pd.DataFrame({
        "name" : [
            "acetic acid",
            "acetone",
            "acetylene"
        ],
        "cid" : [
            176,
            180,
            6326
        ],
        "CAS" :  [
            "6993-75-5",
            "67-64-1",
            "74-86-2"
        ],
        "Smiles" : [
            "CC(=O)O",
            "CC(=O)C",
            "C#C"
        ],
        "molWeight" : [
            60.05,
            58.08,
            26.04
        ],
        "molFormula" : [
            "C2H4O2",
            "C3H6O",
            "C2H2"
        ],
        "logP" : [
            -0.2,
            -0.2,
            0.4
        ],
        "pKa" : [
            4.756,
            20,
            None
        ],
        "charge" : [
            0,
            0,
            0
        ],
        "sterimol" : [

        ]
    })

    assert givesDataFrame("test_data.csv") == result, "givesDataFrame do not give the correct answer"

