# Implementing a Hidden Markov Model Local Ancestry Inference (HMM LAI) Tool

### UCSD CSE 284 Winter 2026
![UCSD Logo](image-url)


### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## 📖 Project Overview

Our **HMM-LAI tool** implements an end-to-end **Hidden Markov Model (HMM)** approach for **local ancestry inference (LAI)** in admixed individuals. 

The objective is to infer ancestry at each genomic position for admixed indivudals using:
* phased genotype data
* ancestry-specific reference panels
* Recombination-aware transition modeling (**double check if we have time**)

Our implementation will replicate the basic functionality of [FLARE (Fast local ancestry estimation)](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00544-4) as published by
> *Browning et al., 2023*
> *The American Journal of Human Genetics*

In the FLARE model:
> Hidden states represent ancestry labels (AFR/EUR/EAS)
> Transitions are driven by recombination distance
> Emissions are based on ancestry-specific allele frequencies estimated from a reference panel

---
## 🚀 Installing our HMM LAI Tool

## 1. Clone this github repository
```bash
git clone https://github.com/nathantrn/CSE284_HMM_local_ancestry_inference.git
cd CSE284_HMM_local_ancestry_inference
```

### ⭐️ Recommended: Create a virtual environment ⭐️

**Using pip**
```python
python3 -m venv hmm_lai_env
```

OR 

**Using conda 🐍**
```bash
conda create --name hmm_lai_env
```

## Activate the environment
**Using pip**
```bash
source hmm_lai_env/bin/activate
```

OR 

**Using conda 🐍**
```bash
conda activate hmm_lai_env
```

## 2. Install Dependecies

We have provided the packages dependencies in requirements.txt

```bash
pip install -r requirements.txt
```

OR

you may manually install teh required Python packages:
(NEED TO EDIT THIS AFTER FINALIZING)

```bash
pip install numpy scipy pandas matplotlub scikit-learn cyvcf2
```
---

## Dataset descriptions

## Reference Set

**(note to self: may use 1000 genomes as reference set since it may be easier to use)**

We downloaded high-coverage sequence data for chromosome 1 from the [Human Genome Diversity Project (HGDP)](ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/).

In line with Browning et al., 2023, we excluded variants that were not bi-allelic SNPs with <1% missingness and at least 5 copies of the minor allele in the combined data. We also assigned panels using the regional labels provided by the HGDP but omitted Oceania due to its smaller size and lack of relevance for the 1000 Genomes data. Additionally, we phased the data using [Beagle 5.2](http://faculty.washington.edu/browning/beagle/beagle.html) with the [HapMap GRCh38 map](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)

---

## Test Set

The test set consistes of chromosome 21 (4.6 GB) **admixed individuals**. These individuals are from the 1000 Genomes Project and not includued in the reference panel.

We downloaded the test data using the following commands:
```bash
wget -c https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz

wget -c https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz.tbi
```

We inferred local ancestry in 6 populations (two for each of AFR, EUR, EAS super populations — **TODO: list the populations here**) with our HMM implementation and compared our output to FLARE for evaluation.

---

## Evaluation Strategy 

### 1. Comparison against a gold standard
Because ground-truth local ancestry labels are unavailable in real datasets, we evaluated our tool's performance against FLARE.

To establish our "gold standard", in line with Browning et al., 2023, we used FLARE with default settings to infer local ancestry in 6 populations (two for each of AFR, EUR, EAS super populations — **TODO: list the populations here**) of the 1000 Genomes project, using a separate analysis for each of these populations. Ancestry proportions were obtained by averaging ancestry calls across sites and individuals.

We compared 
* overall concordance
* per-ancestry concordance
between our model and FLARE where concordance is calculated as


$$Concordance = \frac{\\#\text{ markers where ancestry matches FLARE}}{\text{total \\# markers}}$$

### 2, Toy example comparison of ground truth

To verify correctness of our HMM LAI tool, we simulated small synthetic admixed datasets with known local ancestry labels

Simulation will allow us to validate **(need to double check this!)**
- Transition probability modeling
- Emission probability calculations
- Forward-backward inference correctness

---
## Compuatational Performance

We recorded the:
- Total runtime
- Scaling with number of markers

to assess our tools computational efficiency relative to FLARE
---

## Input Data Requirements

**TODO: UPDATE AS NEEDED**

Our tool expects
* Phased genotype data (VCF)
* 

---

## Limitations

Our HMM LAI tool and FLARE are inference methods that do not provide true ground-truth ancestry labels. Hence, our reported concordance measures agreement with established approaches rather than absolute accuracy.

Our toy simulation is provided for ground-truth validation to ensure correct implementation.

Installing and using the tool, test data set,
