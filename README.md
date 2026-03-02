# Hidden Markov Model Local Ancestry Inference

### UCSD CSE 284 Winter 2026

### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## PROJECT OVERVIEW

This repository implements an end-to-end Hidden Markov Model (HMM) approach for local ancestry inference (LAI). The objective is to infer ancestry at each genomic position for admixed indivudals using phased genotype data and ancestry-specific reference panels.

Our implementation will replicate the basic functionality of [FLARE (Fast local ancestry estimation)](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00544-4), published by Browning et al., 2023 in The American Journal of Human Genetics. Their model treats ancestry along a chromosome as a sequence of hidden states (AFR/EUR/EAS), with transitions driven by recombination distance and emissions based on ancestry-specific allele frequencies estimated from a reference panel.

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

(note to self: may use 1000 genomes as reference set since it may be easier to use)

We downloaded high-coverage sequence data for chromosome 1 from the [Human Genome Diversity Project (HGDP)](ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/).
In line with Browning et al., 2023, we excluded variants that were not bi-allelic SNPs with <1% missingness and at least 5 copies of the minor allele in the combined data. We also assigned panels using the regional labels provided by the HGDP but omitted Oceania due to its smaller size and lack of relevance for the 1000 Genomes data. Additionally, we phased the data using [Beagle 5.2](http://faculty.washington.edu/browning/beagle/beagle.html) with the [HapMap GRCh38 map](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)

---

## Test Set

In line with Browning et al., 2023, we used FLARE with default settings to infer local ancestry in 3 populations (one for each of AFR, EUR, EAS — TODO: list the populations here) of the 1000 Genomes project, using a separate analysis for each of these populations. Ancestry proportions were obtained by averaging ancestry calls across sites and individuals.

---

Installing and using the tool, test data set,
