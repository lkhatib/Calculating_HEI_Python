#!/usr/bin/env python
# coding: utf-8

# # THE HEI-2015 SCORING 
# 
# This code was adapted from the sample SAS code listed in: https://epi.grants.cancer.gov/hei/sas-code.html under 'Simple HEI Scoring Algorithm - Per Day' 'FFQ Example'

# In[2]:


import pandas as pd
import numpy as np


# In[3]:


def hei2015(df, kcal, vtotalleg, vdrkgrleg, f_total, fwholefrt, g_whole, d_total, 
            pfallprotleg, pfseaplantleg, monopoly, satfat, sodium, g_refined, add_sugars):    

    """
    Calculates HEI-2015 component and total scores.

    Parameters
    ----------
    df: pd.DataFrame, required
        The dataset to be used.
    kcal: pd.Series, required
        Series specifying calorie amount.
    vtotalleg: pd.Series, required
        Series specifying the intake of total vegetables plus legumes in cup equivalents.
    vdrkgrleg: pd.Series, required
        Series specifying the intake of dark green vegetables plus legumes in cup equivalents.
    f_total: pd.Series, required
        Series specifying the intake of total fruit in cup equivalents.
    fwholefrt: pd.Series, required
        Series specifying the intake of whole fruit in cup equivalents.
    g_whole: pd.Series, required
        Series specifying the intake of whole grains in ounces equivalents.
    d_total: pd.Series, required
        Series specifying the intake of total dairy in cup equivalents.
    pfallprotleg: pd.Series, required
        Series specifying the intake of total protein foods (including legumes) in ounces equivalents.
    pfseaplantleg: pd.Series, required
        Series specifying the intake of seafood, fish, and plant protein (including legumes) in ounces equivalents.
    monopoly: pd.Series, required
        Series specifying the grams of mono and polyunsaturated fat combined.
    satfat: pd.Series, required
        Series specifying the grams of saturated fat.
    sodium: pd.Series, required
        Series specifying the milligrams of sodium.
    g_refined: pd.Series, required
        Series specifying the intake of refined grains in ounces equivalents.
    add_sugars: pd.Series, required
        Series specifying the intake of added sugars in teaspoons equivalents.

    Returns
    -------
    pd.DataFrame
        A dataset with calculated HEI-2015 component scores, total score, and intermediary density values.
    """

    # Copy input dataset to prevent modification of the original
    outdat = df.copy()

    # Define density calculations
    outdat['VEGDEN'] = np.where(kcal > 0, vtotalleg / (kcal / 1000), 0)
    outdat['GRBNDEN'] = np.where(kcal > 0, vdrkgrleg / (kcal / 1000), 0)
    outdat['FRTDEN'] = np.where(kcal > 0, f_total / (kcal / 1000), 0)
    outdat['WHFRDEN'] = np.where(kcal > 0, fwholefrt / (kcal / 1000), 0)
    outdat['WGRNDEN'] = np.where(kcal > 0, g_whole / (kcal / 1000), 0)
    outdat['DAIRYDEN'] = np.where(kcal > 0, d_total / (kcal / 1000), 0)
    outdat['PROTDEN'] = np.where(kcal > 0, pfallprotleg / (kcal / 1000), 0)
    outdat['SEAPLDEN'] = np.where(kcal > 0, pfseaplantleg / (kcal / 1000), 0)
    outdat['SODDEN'] = np.where(kcal > 0, sodium / kcal, 0)
    outdat['RGDEN'] = np.where(kcal > 0, g_refined / (kcal / 1000), 0)

    # Define fatty acid ratio
    outdat['FARATIO'] = np.where(satfat > 0, monopoly / satfat, 0)

    # Define percentage calculations
    outdat['SFAT_PERC'] = np.where(kcal > 0, 100 * (satfat * 9 / kcal), 0)
    outdat['ADDSUG_PERC'] = np.where(kcal > 0, 100 * (add_sugars * 16 / kcal), 0)

    # Scoring calculations
    def score_density(density, max_score, threshold):
        score = max_score * (density / threshold)
        score = np.where(density == 0, 0, score)
        return np.minimum(score, max_score)

    outdat['HEI2015C1_TOTALVEG'] = score_density(outdat['VEGDEN'], 5, 1.1)
    outdat['HEI2015C2_GREEN_AND_BEAN'] = score_density(outdat['GRBNDEN'], 5, 0.2)
    outdat['HEI2015C3_TOTALFRUIT'] = score_density(outdat['FRTDEN'], 5, 0.8)
    outdat['HEI2015C4_WHOLEFRUIT'] = score_density(outdat['WHFRDEN'], 5, 0.4)
    outdat['HEI2015C5_WHOLEGRAIN'] = score_density(outdat['WGRNDEN'], 10, 1.5)
    outdat['HEI2015C6_TOTALDAIRY'] = score_density(outdat['DAIRYDEN'], 10, 1.3)
    outdat['HEI2015C7_TOTPROT'] = score_density(outdat['PROTDEN'], 5, 2.5)
    outdat['HEI2015C8_SEAPLANT_PROT'] = score_density(outdat['SEAPLDEN'], 5, 0.8)

    # Fatty acid ratio scoring
    FARMIN = 1.2
    FARMAX = 2.5
    outdat['HEI2015C9_FATTYACID'] = np.where(
        (satfat == 0) & (monopoly == 0), 0,
        np.where(
            (satfat == 0) & (monopoly > 0), 10,
            np.where(
                outdat['FARATIO'] >= FARMAX, 10,
                np.where(
                    outdat['FARATIO'] <= FARMIN, 0,
                    10 * (outdat['FARATIO'] - FARMIN) / (FARMAX - FARMIN)
                )
            )
        )
    )

    # Sodium scoring
    SODMIN = 1.1
    SODMAX = 2.0
    outdat['HEI2015C10_SODIUM'] = np.where(
        outdat['SODDEN'] <= SODMIN, 10,
        np.where(
            outdat['SODDEN'] >= SODMAX, 0,
            10 - (10 * (outdat['SODDEN'] - SODMIN) / (SODMAX - SODMIN))
        )
    )

    # Refined grains scoring
    RGMIN = 1.8
    RGMAX = 4.3
    outdat['HEI2015C11_REFINEDGRAIN'] = np.where(
        outdat['RGDEN'] <= RGMIN, 10,
        np.where(
            outdat['RGDEN'] >= RGMAX, 0,
            10 - (10 * (outdat['RGDEN'] - RGMIN) / (RGMAX - RGMIN))
        )
    )

    # Saturated fat scoring
    SFATMIN = 8
    SFATMAX = 16
    outdat['HEI2015C12_SFAT'] = np.where(
        outdat['SFAT_PERC'] >= SFATMAX, 0,
        np.where(
            outdat['SFAT_PERC'] <= SFATMIN, 10,
            10 - (10 * (outdat['SFAT_PERC'] - SFATMIN) / (SFATMAX - SFATMIN))
        )
    )

    # Added sugar scoring
    ADDSUGMIN = 6.5
    ADDSUGMAX = 26
    outdat['HEI2015C13_ADDSUG'] = np.where(
        outdat['ADDSUG_PERC'] >= ADDSUGMAX, 0,
        np.where(
            outdat['ADDSUG_PERC'] <= ADDSUGMIN, 10,
            10 - (10 * (outdat['ADDSUG_PERC'] - ADDSUGMIN) / (ADDSUGMAX - ADDSUGMIN))
        )
    )

    # Calculate total HEI-2015 score
    component_columns = [
        'HEI2015C1_TOTALVEG', 'HEI2015C2_GREEN_AND_BEAN', 'HEI2015C3_TOTALFRUIT',
        'HEI2015C4_WHOLEFRUIT', 'HEI2015C5_WHOLEGRAIN', 'HEI2015C6_TOTALDAIRY',
        'HEI2015C7_TOTPROT', 'HEI2015C8_SEAPLANT_PROT', 'HEI2015C9_FATTYACID',
        'HEI2015C10_SODIUM', 'HEI2015C11_REFINEDGRAIN', 'HEI2015C12_SFAT', 'HEI2015C13_ADDSUG'
    ]
    outdat['HEI2015_TOTAL_SCORE'] = outdat[component_columns].sum(axis=1)

    return outdat

