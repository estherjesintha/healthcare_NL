# Healthcare Facilities in the Netherlands

[![Python](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![Poetry](https://img.shields.io/badge/poetry-1.9.4-lightgrey.svg)](https://python-poetry.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

This repository implements a reproducible geospatial workflow to analyze healthcare accessibility across municipalities in a selected Dutch province. Using open data, APIs, spatial algorithms and visualizations, it identifies underserved areas and high-risk municipalities.

## 1. Clone the Repository

```
git clone https://github.com/estherjesintha/healthcare_NL.git
cd healthcare_NL
```

## 2. Install Poetry

Poetry is used for dependency and virtual environment management.

If Python was installed via Microsoft Store, replace `py` with `python3` in all commands below.

```
py -m pip install --user pipx
```

If you see:  
`WARNING: The script pipx.exe is installed in <USER folder>\AppData\Roaming\Python\Python3x\Scripts which is not on PATH`  
Run:

```
.\pipx.exe ensurepath
```

Restart the terminal and verify:

```
py -m pipx --version
```

Install Poetry via pipx:

```
py -m pipx install poetry
```

Restart the terminal and check:

```
py -m poetry --version
```

[Poetry Documentation](https://python-poetry.org/docs/)

## 3. Create Virtual Environment & Install Dependencies

Configure Poetry to create the virtual environment inside the project:

```
py -m poetry config virtualenvs.in-project true
```

Install dependencies:

```
py -m poetry install
```

## 4. Run Tests

Ensure the environment and code are correctly set up:

```
py -m poetry run pytest
```

## 5. Run the Main Code

Execute the healthcare accessibility analysis:

```
py -m poetry run python src/healthcare_nl/healthcare_NL.py
```


This workflow ensures a reproducible setup for running spatial analysis on healthcare accessibility in the Netherlands. All datasets are programmatically acquired and analyses are fully automated.

