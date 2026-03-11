# Implementing a Hidden Markov Model Local Ancestry Inference (HMM LAI) Tool

### UCSD CSE 284 Winter 2026

### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## 📖 Project Overview

Our **HMM-LAI tool** implements an end-to-end **Hidden Markov Model (HMM)** approach for **local ancestry inference (LAI)** in admixed individuals. 

The objective is to infer ancestry at each genomic position for admixed individuals using:
* phased genotype data
* ancestry-specific reference panels
* recombination-aware transition modeling

Our implementation will replicate the basic functionality of [FLARE (Fast local ancestry estimation)](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00544-4) as published by
> *Browning et al., 2023*, Fast, accurate local ancestry inference with FLARE,
> *The American Journal of Human Genetics*

In the FLARE model:
- Hidden states represent ancestry labels (AFR/EUR/AMR)
- Transitions are driven by recombination distance
- Emissions are based on ancestry-specific allele frequencies estimated from a reference panel

---
## 🚀 Installing and Running our HMM LAI Tool

## 1. Clone this Github Repository
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

First install the required Python packages:

```bash
pip install numpy scipy pandas matplotlib scikit-learn cyvcf2 scikit-allel
```

Also install FLARE, which we use for benchmarking and to generate the model parameters file. FLARE requires Java version 11 or later:

```
wget https://faculty.washington.edu/browning/flare.jar
```

## 3. Run the Tool

First, we recommend running FLARE to get the required [model parameter file](https://github.com/browning-lab/flare/tree/master?tab=readme-ov-file#model-file-format) (i.e. `flare.out.model`), as seen below:

```
java -jar flare.jar ref=flare.ref.vcf.gz gt=flare.gt.vcf.gz map=flare.chr1.map ref-panel=flare.ref.panel out=flare.out
```

Then run the following command, substituting in the appropriate paths. You can specify `-g` to run analysis without the global ancestry calculations or `-l` to run analysis without the local ancestry calculations:

```
python /path/to/CSE284_HMM_local_ancestry_inference/run_hmm.py /path/to/ref/vcf /path/to/ref/panel /path/to/admixed/vcf /path/to/genetic/map /path/to/model/params [-o /output/path] [-g global_off] [-l local_off]
```

Our tool can take a few hours to run depending on the input vcf sizes. If you are on macOS, you can precede the above code with the `caffeinate` command to make sure your computer doesn't turn off while the code is still running.

### Input Data Requirements

Our HMM LAI tool expects:
* Phased reference VCF file
* Tab-separated reference panel mapping samples in reference VCF to reference panel groups
* Phased admixed VCF file
* Genetic map file (cM), i.e. HapMap [GRCh37](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) or [GRCh38](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), provided in the hapMap directory in this repo
* Model parameters, retrieved from running the FLARE tool on the above files

### Expected Output

Running our tool will return two dataframes:

* `ancestry_haplotypes_results.txt`, which contains the per-SNP ancestry predictions for each haplotype in the admixed samples (derived from Viterbi algorithm)
* `global_ancestry_results.txt`, which contains the global ancestry proportions for each admixed samples (derived from forward-backward algorithm)

---

## 📂 Dataset Descriptions

## Reference Set

**(note to self: may use 1000 genomes as reference set since it may be easier to use)**

We downloaded high-coverage sequence data for chromosome 21 from the [Human Genome Diversity Project (HGDP)](ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/) (hg38). The reference panels were made from unadmixed individuals from the AFR and EUR superpopulations and individuals from the American superpopulation in the HGDP dataset. Note that American samples are not truly unadmixed due to the history of colonization and slave trade in the Americas, but for the scope of this project, we are treating these samples as unadmixed for building our reference panel.

Prior to variant filtering and phasing, we merged the 1000 Genomes (**see test set**) and HGDP dataset.

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

The main test set consists of chromosome 21 (4.6 GB) from **admixed individuals** from the 1000 Genomes Project. These individuals **were not includued in the reference panel** and represent admixed individuals from the four populations of the AMR superpopulation (CLM, MXL, PEL, & PUR) which are expected to be admixed with AMR, AFR, and EUR local ancestry. 

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
$$Concordance = \frac{\text{\\# markers where ancestry matches FLARE}}{\text{total \\# markers}}$$
* Per-ancestry concordance
Agreement computed separately for AFR, EUR, and AMR


### 2. Toy Example Comparison of Ground Truth

To verify correctness of our HMM LAI tool, we simulated small synthetic admixed datasets with known local ancestry labels

Simulation will allow us to validate **(need to double check this!)**
- Transition probability modeling
- Emission probability calculations
- Forward-backward inference correctness

---
## ⚡ Computational Performance

We evaluated compuational efficiency by measuring total runtime and peak memory usage.

For FLARE provided test data, FLARE took 2 seconds and had a peak memory footprint of around 120MB.
```
Wallclock Time      :  2 seconds
End Time            :  11:21 PM PDT on 10 Mar 2026
        2.84 real         6.77 user         0.82 sys
           141848576  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               54032  page reclaims
                8844  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   3  signals received
                3632  voluntary context switches
               58393  involuntary context switches
         31013993913  instructions retired
         20515550207  cycles elapsed
           119971136  peak memory footprint
```

In comparison, our tool took about 6 hours and had a peak memory footprint of around 270MB.
```
    22367.63 real     32672.35 user     17126.53 sys
           192794624  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
             2020960  page reclaims
               10299  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   0  signals received
               12138  voluntary context switches
           516968404  involuntary context switches
     478662287689567  instructions retired
     137733740758564  cycles elapsed
           271488000  peak memory footprint
```

For small subset test, FLARE took 2 seconds and had a peak memory footprint of around 140MB.
```
Wallclock Time      :  2 seconds
End Time            :  10:19 AM PDT on 09 Mar 2026
        3.19 real         8.06 user         0.87 sys
           163606528  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               60952  page reclaims
                8854  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   5  signals received
                3494  voluntary context switches
               60919  involuntary context switches
         40079284062  instructions retired
         24032757521  cycles elapsed
           140967296  peak memory footprint
```

In comparison, our tool took around 2 hours and had a peak memory footprint of around 180MB.
```
     6928.40 real     17388.04 user     16373.03 sys
           123551744  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
             1279798  page reclaims
               13496  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   0  signals received
               34571  voluntary context switches
           528245759  involuntary context switches
     204341239799490  instructions retired
      87112638393849  cycles elapsed
           184332608  peak memory footprint
```

---

## ⚠ Limitations of our approach

Our HMM LAI tool and FLARE are inference methods that do not provide true ground-truth ancestry labels. Hence, our reported concordance measures agreement with established approaches rather than absolute accuracy.

Our toy simulation is provided for ground-truth validation to ensure correct implementation.

Additionally, our reference panel contains admixed individuals for the AMR population due to the nature of European colonization and the African slave trade in the Americas. Since AMR individuals in the reference panel are admixed with AMR, AFR, and EUR ancestry, there will be lower accuracy for AMR ancestry because the panel is not 100% unadmixed American. 

---

# FOR PEER REVIEW
We would appreciate feedback on the following
* Readability of the instructions for installing the tool, our goals & metrics
* We were thinking of using only 1000 Genomes data, but if we are using 1000 Genomes for our reference and test, we were thinking of using AMR individuals (who are admixed with AFR, EUR, and possibly EAS), but these individuals also have AMR specific ancestry. How could we alter our approach to accomodate this, or should we stick with using HGDP for our reference and 1000 Genomes for our test?
* Our HMM code (in `run_hmm.py`) is very slow, so any feedback on speeding it up would be great

# REMAINING TASKS
* Create our toy example and run our HMM LAI to assess correctness of implementation
* Run FLARE default parameters on our test data
* Compare FLARE and our HMM LAI tool performance
