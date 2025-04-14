# Protein digest

Python package using a PyO3 Rust implementation of protein digestion using user
specified digestion rules.

The Rust implementation of `get_peptide_to_protein_map` is about four times as 
fast compared to the Python version (7 vs 27 seconds on a fasta file with 190k 
target+decoy sequences).

## Build

Run `maturin develop --release` or `maturin build`.

## Test

1. install pytest: `pip install pytest`
2. run the unit tests: `python -m pytest`

## Project setup

If you do not have a Python installation, run the following first:
```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.9 python3.9-distutils python3.9-venv
python3.9 -m venv .venv
source .venv/bin/activate
```

1. install rust using the command on https://rustup.rs
2. install pipx: `pip install pipx`
3. install maturin: `pipx install maturin`
4. initialize package with maturin: `maturin new`