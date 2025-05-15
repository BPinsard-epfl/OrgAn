import numpy as np
import matplotlib
import pandas as pd

from functions import givesDataFrame
from pchem_rq import getMoleculeInfoFromSmiles

def get_elution_order(solutes: pd.DataFrame, is_reverse_phase = False): # TODO: solutes as dataframe to get logP and thus elution order
    solutes = solutes.sort_values("logP", ascending=is_reverse_phase).reset_index(drop=True)
    #elution_order = []
    #for i in range(len(solutes)):
    #    elution_order.append(solutes.iloc[i,0])
    return solutes


def estimate_retention_factor(logP, polarity_index):
    logk = logP - 0.5 * (10.2 - polarity_index)
    return 10 ** logk

def calculate_polarity_index( # is it better to use abbreviations or the full name of the solvents for the arguments?
        cyclohex = 0, n_hex = 0, ccl4 = 0, ipr_ether = 0, 
        toluene = 0, et2o = 0, thf = 0, etoh = 0, etoac = 0, 
        dioxane = 0, meoh = 0, mecn = 0, water = 0) -> float:
    """
    Calculates the polarity index of the solvent based on the table provided here: 
    https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Instrumental_Analysis_(LibreTexts)/28%3A_High-Performance_Liquid_Chromatography/28.04%3A_Partition_Chromatography
    
    The user inputs the volume fraction of the components of the solvent, the total must be equal to 1.

    Parameters
    ----------
    cyclohex : float
        The volume fraction of cyclohexane in the solvent.
    n_hex : float
        The volume fraction of n-hexane in the solvent.
    ccl4 : float
        The volume fraction of carbon tetrachloride in the solvent.
    ipr_ether : float
        The volume fraction of isopropyl ether in the solvent.
    toluene : float
        The volume fraction of toluene in the solvent.
    et2o : float
        The volume fraction of diethyl ether in the solvent.
    thf : float
        The volume fraction of tetrahydrofuran in the solvent.
    etoh : float
        The volume fraction of ethanol in the solvent.
    etoac : float
        The volume fraction of ethyl acetate in the solvent.
    dioxane : float
        The volume fraction of dioxane in the solvent.
    meoh : float
        The volume fraction of methanol in the solvent.
    mecn : float
        The volume fraction of acetonitrile in the solvent.
    water : float
        The volume fraction of water in the solvent.

    Returns
    -------
    float
        The polarity index of the solvent.
    """
    if cyclohex + n_hex + ccl4 + ipr_ether + toluene + et2o + thf + etoh + etoac + dioxane + meoh + mecn + water != 1:
        raise ValueError("The total of the volume fractions must be equal to 1")
    return 0.04 * cyclohex + 0.1 * n_hex + 1.6 * ccl4 + 2.4 * ipr_ether + 2.4 * toluene + 2.8 * et2o + 4.0 * thf + 4.3 * etoh + 4.4 * etoac + 4.8 * dioxane + 5.1 * meoh + 5.8 * mecn + 10.2 * water

# test functions
#data = givesDataFrame("tests/test_data.csv")
#print(data)
#print(get_elution_order(data))
#print(calculate_polarity_index(water=0.6,meoh=0.4))