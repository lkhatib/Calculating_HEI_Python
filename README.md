# Calculating\_HEI\_Python

## Overview

This repository provides a Python script for calculating the Healthy Eating Index (HEI), an important measure of diet quality. The code is adapted from the sample SAS code available at the [NIH HEI website](https://epi.grants.cancer.gov/hei/sas-code.html) under the "Simple HEI Scoring Algorithm - Per Day" section for the FFQ example.

The calculations rely on MyPyramid Equivalents Database 2.0 (MPED) components, which are derived from Food Frequency Questionnaire (FFQ) data. Learn more about MPED components [here](https://www.whi.org/md/diet-dataset-bottom).

---

## Parameters

The HEI calculation script requires the following parameters:

| Parameter       | Type           | Description                                                                                    |
| --------------- | -------------- | ---------------------------------------------------------------------------------------------- |
| `df`            | `pd.DataFrame` | The dataset to be used.                                                                        |
| `kcal`          | `pd.Series`    | Series specifying calorie amount.                                                              |
| `vtotalleg`     | `pd.Series`    | Series specifying the intake of total vegetables plus legumes (cup eq.).                       |
| `vdrkgrleg`     | `pd.Series`    | Series specifying the intake of dark green vegetables plus legumes (cup eq.).                  |
| `f_total`       | `pd.Series`    | Series specifying the intake of total fruit (cup eq.).                                         |
| `fwholefrt`     | `pd.Series`    | Series specifying the intake of whole fruit (cup eq.).                                         |
| `g_whole`       | `pd.Series`    | Series specifying the intake of whole grains (oz. eq.).                                        |
| `d_total`       | `pd.Series`    | Series specifying the intake of total dairy (cup eq.).                                         |
| `pfallprotleg`  | `pd.Series`    | Series specifying the intake of total protein foods, including legumes (oz. eq.).              |
| `pfseaplantleg` | `pd.Series`    | Series specifying the intake of seafood, fish, and plant protein, including legumes (oz. eq.). |
| `monopoly`      | `pd.Series`    | Series specifying the grams of mono- and polyunsaturated fats combined.                        |
| `satfat`        | `pd.Series`    | Series specifying the grams of saturated fat.                                                  |
| `sodium`        | `pd.Series`    | Series specifying the milligrams of sodium.                                                    |
| `g_refined`     | `pd.Series`    | Series specifying the intake of refined grains (oz. eq.).                                      |
| `add_sugars`    | `pd.Series`    | Series specifying the intake of added sugars (teaspoons eq.).                                  |

---

## Returns

The script outputs a `pd.DataFrame` with:

- Calculated HEI-2015 component scores
- Total HEI-2015 score
- Intermediary density values for each component

---

## Files

### `calculating_HEI_Spain.ipynb`

This Jupyter notebook provides an example of calculating MPED scores from FFQ responses.

### `HEI_2015_Scoring.py`

The main script for calculating HEI-2015 scores based on the input parameters.

---

## Usage

To use the HEI calculation script, ensure that your dataset includes the required parameters as separate columns, formatted as described above. Then, run the script or adapt the example notebook to calculate HEI scores for your dataset.

For detailed guidance on the HEI methodology, refer to the official [NIH documentation](https://epi.grants.cancer.gov/hei/).

