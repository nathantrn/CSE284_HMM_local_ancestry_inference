# HaploHMM: A Hidden Markov Model Local Ancestry Inference (HMMLAI) Tool

### UCSD CSE 284 Winter 2026

### Authors: Jonah Pacis, Safa Saeed, Nathan Tran

## 📖 Project Overview

Our **HMM-LAI tool HaploHMM** implements an end-to-end **Hidden Markov Model (HMM)** approach for **local ancestry inference (LAI)** in admixed individuals. 

The objective is to infer ancestry at each genomic position for admixed individuals using:
* phased genotype data
* ancestry-specific reference panels
* recombination-aware transition modeling

Our implementation will replicate the basic functionality of [FLARE (Fast local ancestry estimation)](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00544-4) as published by
> *Browning et al., 2023*, Fast, accurate local ancestry inference with FLARE,
> *The American Journal of Human Genetics*

In the FLARE model:
- Hidden states represent ancestry labels (ex. AFR/EUR/AMR) paired with reference haplotype labels
- Transitions are driven by recombination distance
- Emissions are based on ancestry-specific allele frequencies estimated from a reference panel

---
## 🚀 Installing and Running HaploHMM

## 1. Clone this GitHub Repository
```bash
git clone https://github.com/nathantrn/CSE284_HMM_local_ancestry_inference.git
cd CSE284_HMM_local_ancestry_inference
```

### ⭐️ Recommended: Create a virtual environment ⭐️

**Using pip🐥**
```python
python3 -m venv haploHMM
```

OR 

**Using conda 🐍**
```bash
conda create --name haploHMM
```

## Activate the environment
**Using pip🐥**
```bash
source haploHMM/bin/activate
```

OR 

**Using conda 🐍**
```bash
conda activate haploHMM
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

HaploHMM expects:
* Phased reference VCF file
* Tab-separated reference panel mapping samples in reference VCF to reference panel groups
* Phased admixed VCF file
* Genetic map file (cM), i.e. HapMap [GRCh37](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) or [GRCh38](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), provided in the hapMap directory in this repo
* Model parameters, retrieved from running the FLARE tool on the above files

**IMPORTANT:** Ensure that the reference and admixed VCF files contain the same variant positions, otherwise HaploHMM will crash. If that isn't the case, you can run `bcftools isec` beforehand to ensure that the same locations are present in both files.

### Expected Output

Running our tool will return two dataframes:

* `ancestry_haplotypes_results.txt`, which contains the per-SNP ancestry predictions for each haplotype in the admixed samples (derived from Viterbi algorithm)
* `global_ancestry_results.txt`, which contains the global ancestry proportions for each admixed samples (derived from forward-backward algorithm)

---

## 📂 Dataset Descriptions

### 1. Toy Example Comparison of Ground Truth

To verify correctness of the HaploHMM tool, we simulated small synthetic admixed datasets with known local ancestry labels (AFR: YRI, EUR: IBS, AMR: PEL).

We generated 3 test cases using haptools:

1. unadmixed individual (100% YRI)

2. 2-way admixed individual (80% YRI, 20% IBS)

3. 3-way admixed individual (40% YRI, 40% IBS, 20% PEL)

We ran FLARE and HaploHMM given the following inputs:
- A subsetted dataset of phased haplotypes from "unadmixed" YRI, IBS, and PEL individuals: `toy_example/1000G_chr21_subset_downsampled.vcf.gz`
- A subsetted reference panel from the 1000 Genomes data for YRI, IBS, and PEL individuals: `toy_example/1000genomes_sampleinfo_subset.tsv`
- A merged dataset of phased haplotypes from the 3 types of test cases:
`toy_example/toy_examples_haptools/toy_example_test.vcf.gz`
- A genetic map for chr21:
`toy_example/plink.chr21.GRCh38_renamed.map`

We compared FLARE and HaploHMM's accuracy, which was calculated as:
$$Accuracy = \frac{\text{\\# correctly assigned SNPs}}{\text{total \\# SNPs}}$$

Simulation will allow us to validate **(need to double check this!)**
- Transition probability modeling
- Emission probability calculations
- Forward-backward inference correctness

### 2. Flare Test Dataset
Our next test set comes from FLARE's GitHub repository. This includes a reference set with three ancestries, Panel.A, Panel.B, and Panel.C, each panel consisting of 100 samples, for a total of 300 reference samples. There are two admixed samples, msp_300 and msp_301. Each sample has phased genotype data for 2000 variant positions.

We ran FLARE and HaploHMM on this dataset, calculating concordance for each haplotype's local ancestry assignments and comparing global ancestry predictions between the two methods. Concordance was calculated as:
$$Concordance = \frac{\text{\\# markers where ancestry matches FLARE}}{\text{total \\# markers}}$$

### 3. 1000Genomes Dataset
Samples for both the Reference panel and the Test set are derived from the 1000 Genomes Project Phase 3 dataset for chromosome 21. Variant data can be downloaded using:

```bash
wget -c
https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

wget -c https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
```
Metadata can be found in `/1000genomes/igsr_samples.tsv`

#### Variant Filtering
In line with Browning et al., 2023, we:
* excluded variants that were not bi-allelic SNPs
* excluded variants with >1% missingness 
* required at least 5 copies of the minor allele 

#### Reference Set

The reference panel was constructed using unadmixed individuals from three superpopulations:
- AFR (African Ancestry)
- EUR (European Ancestry)
- EAS (East Asian Ancestry)


#### Test Set

The admixed test set consists of individuals from the AMR (Admixed American) superpopulation in the 1000 Genomes dataset. These individuals are expected to contain mixed ancestry components from multiple superpopulations and are thus suitable for evaluating local ancestry inference. 

#### Data Prep.ipynb
`Data Prep.ipynb` provides a step-by-step guide to creating a subsetted Reference panel and Test set. Precomputed, ready-to-go files are housed in `/1000genomes/small_subset_data_prep/`

### Global Ancestry Comparisons
We created plots to compare global ancestry proportions for each of our test datasets. These can be found in the `plot_results.ipynb` notebook of this repo.

---

## ⚡ Computational Performance

We evaluated compuational efficiency by measuring total runtime and peak memory usage, using the `time` command line tool.

For the haptools simulated toy example, FLARE took 3 seconds and had a peak memory footprint of around 190MB.
```
Wallclock Time      :  3 seconds
End Time            :  06:54 PM PDT on 11 Mar 2026
        3.80 real        10.36 user         1.05 sys
           212979712  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
               78405  page reclaims
                9180  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                  17  signals received
                3744  voluntary context switches
               82157  involuntary context switches
         52371653119  instructions retired
         30715565997  cycles elapsed
           192806656  peak memory footprint
```
In comparison, our tool took around 1 hour and had a peak memory footprint of around 120MB.
```
     3598.12 real      7267.81 user      9564.98 sys
           104079360  maximum resident set size
                   0  average shared memory size
                   0  average unshared data size
                   0  average unshared stack size
              443381  page reclaims
               10503  page faults
                   0  swaps
                   0  block input operations
                   0  block output operations
                   0  messages sent
                   0  messages received
                   0  signals received
                7913  voluntary context switches
           601936903  involuntary context switches
      98220162255567  instructions retired
      44216095399232  cycles elapsed
           119553920  peak memory footprint
```

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


