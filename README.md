# OrgAn – Molecular Analysis and Chromatographic Simulation

**OrgAn** is a Python-based framework designed to extract, analyze, and simulate molecular properties using cheminformatics tools. Developed in an academic context, it automates the retrieval of molecular data, performs structural and physicochemical analysis, and simulates chromatographic behavior.

## Features

- Extracts molecular descriptors from the PubChem database (e.g., logP, pKa, molecular weight, Sterimol parameters)
- Builds structured datasets from user-supplied SMILES strings
- Identifies key gaps in properties like logP and pKa across compound sets
- Suggests compounds based on custom physicochemical criteria
- Simulates chromatographic retention using logP and solvent polarity

## Project Structure

```
OrgAn-main/
├── report.ipynb           # Main Jupyter notebook
├── data/                  # Folder for input CSV files
├── src/                   # Python modules (PubChem requests, functions, chromatography)
│   ├── pchem_rq.py
│   ├── functions.py
│   └── chromato.py
├── tests/                 # Test suite
├── env.yml                # Conda environment specification
├── pyproject.toml         # Build configuration
└── LICENSE
```

## Installation

1. Clone this repository:

```bash
git clone https://github.com/your-username/OrgAn.git
cd OrgAn
```

2. Create the environment with Conda:

```bash
conda env create -f env.yml
conda activate organ
```

3. (Optional) Install the package locally in editable mode:

```bash
pip install -e .
```

## Input Format

The input must be a `.csv` file with a single column labeled `smiles`. Each row should contain one SMILES string:

```
smiles
CCO
c1ccccc1
CC(=O)O
```

## Usage

Open and run `report.ipynb` step-by-step. The notebook will:

- Load and process the input SMILES
- Query PubChem for descriptors
- Build a DataFrame of molecular data
- Perform analysis and chromatography simulation

All results are shown interactively within the notebook.

## Authors

- Raphaël Tisseyre  
- Bastien Pinsard  
- Johan Schmidt

## License

This project is released under the MIT License. See the `LICENSE` file for details.
