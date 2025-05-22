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

Follow these steps to correctly set up your development environment before running the app.

### 0. Clone this GitHub repository

Start by cloning the project to your local machine:

```bash
git clone https://github.com/BPinsard-epfl/OrgAn.git
cd OrgAn
```

### 1. Install Python 3.10

This project requires **Python 3.10**. Download it here:  
[https://www.python.org/downloads/release/python-3100/](https://www.python.org/downloads/release/python-3100/)

### 2. Install Conda

If you haven’t installed Conda yet, choose one of the following:

- **Miniconda (recommended)**: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
- or **Anaconda**: [https://www.anaconda.com/products/distribution](https://www.anaconda.com/products/distribution)

### 3. Create the development environment

You have two options:

#### Option 1 — Using the `env.yml` file (recommended)

1. Run the following commands from the project directory:

```bash
conda env create -f env.yml
conda activate ppchem_project
```

> This will ensure all required packages and compatible versions are installed.

#### Option 2 — Manual setup with `pip`

1. Create a new Conda environment with Python 3.10:

```bash
conda create -n ppc_env python=3.10
conda activate ppc_env
```

2. Install required packages:

```bash
pip install streamlit streamlit_ketcher attrs==25.3.0 certifi==2025.4.26 charset-normalizer==3.4.1 exceptiongroup==1.2.2 h11==0.16.0 html5lib==1.1 idna==3.10 mechanize==0.4.10 numpy==2.2.5 outcome==1.3.0.post0 pandas==2.2.3 pillow==11.2.1 pubchempy==1.0.4 pysocks==1.7.1 python-dateutil==2.9.0.post0 pytz==2025.2 rdkit==2024.9.6 regex==2024.11.6 requests==2.32.3 selenium==4.31.0 six==1.17.0 sniffio==1.3.1 sortedcontainers==2.4.0 trio==0.30.0 trio-websocket==0.12.2 typing-extensions==4.13.2 tzdata==2025.2 urllib3==2.4.0 webencodings==0.5.1 websocket-client==1.8.0 wsproto==1.2.0
```

### Step 4 — Run the application

To start the app, run the following command in the project directory:

```bash
streamlit run app.py
```

---

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
Loads a `.csv` file containing a column `smiles`, retrieves descriptors (logP, pKa, MW, etc.) for each molecule, and returns a DataFrame.

```python
from src.functions import gives_data_frame
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

You must provide a `.csv` file with a `smiles` column:

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
