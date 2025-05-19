import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_elution_order(solutes: pd.DataFrame, is_reverse_phase = False) -> pd.DataFrame:
    solutes = solutes.sort_values("logP", ascending=is_reverse_phase).reset_index(drop=True)
    return solutes


def estimate_retention_factor(logP : float , polarity_index : float) -> float:
    logk = logP - 0.5 * (10.2 - polarity_index)
    return 10 ** logk


def calculate_polarity_index( # is it better to use abbreviations or the full name of the solvents for the arguments?
        cyclohex : float = 0, n_hex : float = 0, ccl4 : float  = 0,
        ipr_ether : float = 0, toluene : float = 0, et2o : float = 0,
        thf : float = 0, etoh : float = 0, etoac : float = 0, dioxane : float = 0,
        meoh : float = 0, mecn : float = 0, water : float = 0) -> float:
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


def generate_chromatogram(df : pd.DataFrame, polarity_idx : float, dead_time : float = 1, savefigas : str = "") -> None:
    
    df["retention time"] = df["logP"].apply(lambda x : (estimate_retention_factor(x, polarity_idx)+1)*dead_time)

    x_time = np.array([x/100 for x in range(round(df["retention time"].iloc[-1]*100)+100)])
    y_signal = [0 for x in range(round(df["retention time"].iloc[-1]*100)+100)]

    index_df = 0
    index_array = 0
    while index_df<len(df):
        if df[index_df]<x_time[index_array]:
            index_array +=1
        else :
            if abs(df[index_df]-x_time[index_array])<abs(x_time[index_array]-x_time[index_array-1]):
                y_signal[index_array] = 1
                index_df += 1
            else :
                y_signal[index_array-1] = 1
                index_df += 1
        index_array += 1
    
    plt.plot(x_time, y_signal, "red", linewidth= 1.25, label= "Estimation")
    plt.xlabel("Time [min]")
    plt.ylabel("Signal")
    plt.title("Estimated Chromatogram")

    if savefigas :
        plt.savefig(savefigas)

    plt.show()