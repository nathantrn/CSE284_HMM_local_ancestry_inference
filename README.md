# Implementing a Hidden Markov Model Local Ancestry Inference (HMM LAI) Tool

### UCSD CSE 284 Winter 2026

### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## 📖 Project Overview

Our **HMM-LAI tool** implements an end-to-end **Hidden Markov Model (HMM)** approach for **local ancestry inference (LAI)** in admixed individuals. 

The objective is to infer ancestry at each genomic position for admixed individuals using:
* phased genotype data
* ancestry-specific reference panels
* Recombination-aware transition modeling (**double check if we have time**)

Our implementation will replicate the basic functionality of [FLARE (Fast local ancestry estimation)](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00544-4) as published by
> *Browning et al., 2023*, Fast, accurate local ancestry inference with FLARE,
> *The American Journal of Human Genetics*

In the FLARE model:
- Hidden states represent ancestry labels (AFR/EUR/AMR)
- Transitions are driven by recombination distance
- Emissions are based on ancestry-specific allele frequencies estimated from a reference panel

---
## 🚀 Installing our HMM LAI Tool

## 1. Clone this github repository
```bash
git clone https://github.com/nathantrn/CSE284_HMM_local_ancestry_inference.git
cd CSE284_HMM_local_ancestry_inference
```

### ⭐️ Recommended: Create a virtual environment ⭐️

**Using pip🐥**
```python
python3 -m venv hmm_lai_env
```

OR 

**Using conda 🐍**
```bash
conda create --name hmm_lai_env
```

## Activate the environment
**Using pip🐥**
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

you may manually install the required Python packages:
(**NEED TO EDIT THIS AFTER FINALIZING**)

```bash
pip install numpy scipy pandas matplotlub scikit-learn cyvcf2
```
---

## 📂 Dataset descriptions

## Reference Set

**(note to self: may use 1000 genomes as reference set since it may be easier to use)**

We downloaded high-coverage sequence data for chromosome 21 from the [Human Genome Diversity Project (HGDP)](ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/) (hg38). The reference panels were made from unadmixed individuals from the AFR and EUR superpopulations and individuals from the American superpopulation in the HGDP dataset. Note that American samples are not truly unadmixed due to the history of colonization and slave trade in the Americas, but for the scope of this project, we are treating these samples as unadmixed for building our reference panel.

Prior to variant filtering and phasing, we merged the 1000 Genomes (see test set) and HGDP dataset.

### Variant Filtering
In line with Browning et al., 2023, we:
* excluded variants that were not bi-allelic SNPs
* excluded variants with >1% missingness 
* required at least 5 copies of the minor allele 
* Omitted Oceania due to its smaller size and lack of relevance for the 1000 Genomes data.

### Phasing
We phased the data using [Beagle 5.2](http://faculty.washington.edu/browning/beagle/beagle.html) with the [HapMap GRCh38 map](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip)

---

## Test Set

The main test set consistes of chromosome 21 (4.6 GB) from **admixed individuals** from the 1000 Genomes Project. These individuals **were not includued in the reference panel** and represent admixed individuals from the four populations of the AMR superpopulation (CLM & PUR) which are expected to be admixed with AMR, AFR, and EUR local ancestry. 

Data was downloaded using:
```bash
wget -c https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz

wget -c https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr21.vcf.gz.tbi
```

We inferred local ancestry in the four AMR populations (CLM, MXL, PUR, PEL) with our HMM implementation and compared our output to FLARE for evaluation. Additionally, we inferred local ancestry for one population from the EAS and SAS superpopulations in the 1000 Genomes Dataset (EAS: JPT, SAS: PJL) to evaluate performance on samples where the true ancestry is absent from our model. We expect LAI to be uncertain/unstable for these populations not represented in our model.

---

## 📊 Evaluation Strategy 

### 1. Comparison with FLARE (Reference Benchmark)
Because ground-truth local ancestry labels are unavailable in real datasets, we evaluated our tool's performance against FLARE using default settings.

We computed:
* Overall concordance
$$Concordance = \frac{\\#\text{ markers where ancestry matches FLARE}}{\text{total \\# markers}}$$
* Per-ancestry concordance
Agreement computed separately for AFR, EUR, and EAS


### 2. Toy example comparison of ground truth

To verify correctness of our HMM LAI tool, we simulated small synthetic admixed datasets with known local ancestry labels

Simulation will allow us to validate **(need to double check this!)**
- Transition probability modeling
- Emission probability calculations
- Forward-backward inference correctness

---
## ⚡ Computational Performance

We evaluated compuational efficieny by measuring:
- Total runtime
- Scaling with number of markers
to assess our tools computational efficiency relative to FLARE
---

## Input Data Requirements

**TODO: UPDATE AS NEEDED**

Our HMM LAI tool expects
* Phased genotype data (VCF format)
* Genetic map file (cM or Morgans) (**double check this**)
* Reference sample list (non-admixed individuals)
* Test sample list (admixed individuals)

---

## ⚠ Limitations of our approach

Our HMM LAI tool and FLARE are inference methods that do not provide true ground-truth ancestry labels. Hence, our reported concordance measures agreement with established approaches rather than absolute accuracy.

Our toy simulation is provided for ground-truth validation to ensure correct implementation.

---

# FOR PEER REVIEW
We would appreciate feedback on the following
* readability of the instructions for installing the tool, our goals & metrics
* We were thinking of using only 1000 Genomes data, but if we are using 1000 Genomes for our reference and test, we were thinking of using AMR individuals (who are admixed with AFR, EUR, and possibly EAS), but these individuals also have AMR specific ancestry. How could we alter our approach to accomodate this, or should we stick with using HGDP for our reference and 1000 Genomes for our test?

# REMAINING TASKS
* Implement out HMM LAI using Pomegranate forward backward algorithm
* Create our toy example and run our HMM LAI to assess correctness of implementation
* Run FLARE default parameters on our test data
* Compare FLARE and our HMM LAI tool performance
