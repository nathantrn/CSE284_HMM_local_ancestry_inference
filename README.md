# Hidden Markov Model Local Ancestry Inference

### UCSD CSE 284 Winter 2026

### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## PROJECT OVERVIEW

This repository implements an end-to-end Hidden Markov Model (HMM) approach for local ancestry inference (LAI). The objective is to infer ancestry at each genomic position for admixed indivudals using phased genotype data and ancestry-specific reference panels.

Our implementation will replicate the basic functionality of FLARE. The model treats ancestry along a chromosome as a sequence of hidden states (AFR/EUR/EAS), with transitions driven by recombination distance and emissions based on ancestry-specific allele frequencies estimated from a reference panel.

---
## INSTALLATION

## 1. Clone this github repository
```bash
git clone https://github.com/[insert the rest of path].git
cd [name of repo]
```

### ⭐️ Recommended: Create a virtual environment ⭐️

Using pip
```python
python3 -m venv hmm_lai_env
```

OR 

Using conda 🐍
```bash
conda create --name hmm_lai_env
```

## Activate the environment
Using pip
```python
source hmm_lai_env/bin/activate
```

OR 

Using conda 🐍
```bash
conda activate hmm_lai_env
```
---

## DATASET DESCRIPTION

## Reference Set

---

## Test Set



Installing and using the tool, test data set,
