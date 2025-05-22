# OrgAn – Molecular Analysis and Chromatographic Simulation

**OrgAn** is a Python-based tool for analyzing molecular properties and simulating chromatographic behavior. It is built to process SMILES data, fetch molecular descriptors, and generate both quantitative comparisons and visual analyses.

---

## Features

OrgAn allows you to:

- Extract molecular descriptors (logP, pKa, molecular weight, etc.) from public sources like PubChem
- Build structured molecular datasets from CSV input
- Identify gaps in descriptors across compound sets
- Suggest molecules based on custom physicochemical criteria
- Simulate chromatographic retention behavior
- Retrieve and process molecular structures for analysis

---

## Installation

OrgAn relies on Numpy, Pandas, RDKit, Morfeus, Requests, Regex, Matplotlib

A virtual environment with Python 3.10 can be installed with Anaconda as follows:

```
conda create -n organ python=3.10
conda activate organ
```

#### Installing from source

The first step is to clone the represitory in your local device:
```bash
git clone https://github.com/BPinsard-epfl/OrgAn.git
cd OrgAn
```
For simple user :
```
pip install .
```
Or for developpers :
```
pip install -e .
```

## Project Structure

```
OrgAn/
├── report.ipynb           # Documentation and usage examples
├── data/                  # Input CSV files
├── src/                   # Core logic and modules
│   ├── pchem_rq.py
│   ├── functions.py
│   └── chromato.py
├── app.py                 # Streamlit web app
├── tests/                 # Unit tests
├── pyproject.toml         # Build config
└── LICENSE
```

---

## Usage

Once installed, you can use OrgAn through Python scripts.

### Core Function Usage

Here are some of the key functions and how to use them:

#### `gives_data_frame(path: str) -> pd.DataFrame`
Loads a `.csv` file containing a column `smiles`, retrieves descriptors (logP, pKa, MolWeight, etc.) for each molecule, and returns a DataFrame.

```python
from OrgAn import gives_data_frame
df = gives_data_frame("data/example_smiles.csv")
```

#### `logp_gap(df: pd.DataFrame, threshold: float = 0.5) -> pd.DataFrame`
Identifies gaps in logP values across a compound set. Useful for optimizing compound diversity.

```python
from src.chromato import logp_gap
gaps = logp_gap(df)
```

#### `plot_gaps(df: pd.DataFrame)`
Displays a chromatographic-style plot that highlights the logP gaps between molecules.  
This helps visualize the separation potential of a compound set and supports chromatographic simulation.

```python
from src.chromato import plot_gaps
plot_gaps(df)
```

#### `find_compound(logp: float, pka: float)`
Returns compound suggestions based on a desired logP and pKa.

```python
from src.functions import find_compound
suggested = find_compound(logp=2.5, pka=4.8)
```

---

## Input Format

You must provide a `.csv` file with a `smiles` column like the example below:

```
smiles
CCO
c1ccccc1
CC(=O)O
```

---

## Authors

- Raphaël Tisseyre  
- Bastien Pinsard  
- Johan Schmidt

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
