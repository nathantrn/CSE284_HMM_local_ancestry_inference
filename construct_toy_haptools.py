"""
Toy Example using haptools simgenotype
==================================================
Generates simulated admixed genotype data with ground-truth local ancestry
labels (breakpoints) using haptools simgenotype.

haptools outputs:
  - <scenario>.vcf.gz        — phased genotypes for simulated admixed individuals
  - <scenario>.bp            — ground-truth ancestry breakpoints (the "true labels")


Prior to running this script:
-------------
  pip install haptools

 This script assumes you have:
    1. A phased reference VCF (chr21 region from 1000G/HGDP data)
    2. The HapMap GRCh38 genetic map for chr21
    3. A sample_info .tsv mapping sample IDs -> population codes

  See INPUTS section below to set these paths before running.

Scenarios
---------
    Trivial         — 100% YRI (AFR) as an unadmixed control
    Two ancestries  — AFR + EUR admixture, 7 generations
    Three ancestries — AFR + EUR + AMR admixture, 7 generations

Usage
-----
  python construct_toy_haptools.py

  Outputs written to ./toy_examples_haptools/
"""

import subprocess
import textwrap
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# INPUTS — set these to match your local file paths before running
# ─────────────────────────────────────────────────────────────────────────────

# Phased reference VCF (must be tabix-indexed .vcf.gz)
# Use your chr21 reference panel — make sure it includes YRI, IBS/GBR, and PEL samples
REF_VCF = "data/reference_chr21.vcf.gz"

# Directory containing HapMap GRCh38 .map files
# File for chr21 must be named so it contains "chr21" in the filename, e.g.:
#   plink.chr21.GRCh38.map
# Download from: http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
MAP_DIR = "data/genetic_map_GRCh38/"

# Sample info TSV: two columns, no header — sample_id <tab> population_code
# Population codes must match exactly what you use in the model files below
# (YRI, IBS, PEL — or whatever populations are in your reference VCF)
# haptools provides a ready-made 1000G version:
#   wget https://raw.githubusercontent.com/CAST-genomics/haptools/main/example-files/1000genomes_sampleinfo.tsv
SAMPLE_INFO = "data/1000genomes_sampleinfo.tsv"

# Region to simulate (subset of chr21 for speed — adjust to your SNP range)
REGION = "chr21:10000000-15000000"

# Number of admixed individuals to simulate per scenario
N_SAMPLES = 20 #can also be 1?

# Output directory
OUT_DIR = Path("toy_examples_haptools")

# ─────────────────────────────────────────────────────────────────────────────
# Model file construction
# ─────────────────────────────────────────────────────────────────────────────
# The .dat format (inherited from admix-simu, used directly by haptools):
#
#   Line 1 (header):  <n_haplotypes>  Admixed  <POP1>  <POP2>  ...
#   Line 2+:          <generation>    <p_admixed>  <p_POP1>  <p_POP2>  ...
#
# - n_haplotypes = N_SAMPLES * 2  (diploid)
# - p values on each line must sum to 1.0
# - generation 1 = the founding admixture event
# - p_admixed on line 1 must be 0 (no admixed ancestors yet in gen 1)
# - subsequent lines can draw from admixed pool to model ongoing admixture
#
# Population codes must exactly match what is in SAMPLE_INFO.
# ─────────────────────────────────────────────────────────────────────────────

def write_model(path: Path, content: str):
    path.write_text(textwrap.dedent(content).strip() + "\n")
    print(f"    model file : {path}")


def build_model_files(model_dir: Path, n_haplotypes: int):
    """Write .dat model files for each scenario."""
    model_dir.mkdir(parents=True, exist_ok=True)

    # ── S1: Unadmixed control — 100% YRI for 1 generation ──
    # All haplotypes draw from YRI. p_admixed=0, p_YRI=1.0.
    # This produces a flat karyogram — everything should be called AFR.
    write_model(model_dir / "S1_trivial.dat", f"""
        {n_haplotypes}\tAdmixed\tYRI
        1\t0\t1.0
    """)

    # ── S2: Two ancestries — AFR + EUR pulse, 7 generations ──
    # Generation 1: 80% YRI (AFR), 20% IBS (EUR).
    # Generations 2-7: purely admixed (draws from admixed pool only).
    # 7 generations ≈ ~175-200 years, appropriate for Americas admixture.
    gens_two = "\n".join(
        f"        {g}\t1.0\t0\t0" for g in range(2, 8)
    )
    write_model(model_dir / "S2_two_ancestries.dat", f"""
        {n_haplotypes}\tAdmixed\tYRI\tIBS
        1\t0\t0.8\t0.2
{gens_two}
    """)

    # ── S3: Three ancestries — AFR + EUR + AMR, 7 generations ──
    # Generation 1: 40% YRI (AFR), 40% IBS (EUR), 20% PEL (AMR proxy).
    # Generations 2-7: purely admixed.
    # PEL (Peruvians) is used as the AMR reference due to high proportion of indigenous AMR ancestry
    gens_three = "\n".join(
        f"        {g}\t1.0\t0\t0\t0" for g in range(2, 8)
    )
    write_model(model_dir / "S3_three_ancestries.dat", f"""
        {n_haplotypes}\tAdmixed\tYRI\tIBS\tPEL
        1\t0\t0.4\t0.4\t0.2
{gens_three}
    """)


# ─────────────────────────────────────────────────────────────────────────────
# Run haptools simgenotype
# ─────────────────────────────────────────────────────────────────────────────

def run_simgenotype(scenario_name: str, model_file: Path, out_prefix: Path):
    """
    Run haptools simgenotype for one scenario.

    Outputs:
      <out_prefix>.vcf.gz   — simulated phased genotypes
      <out_prefix>.bp       — ground-truth ancestry breakpoints
    """
    out_vcf = out_prefix.with_suffix(".vcf.gz")

    cmd = [
        "haptools", "simgenotype",
        "--model",       str(model_file),
        "--mapdir",      MAP_DIR,
        "--ref_vcf",     REF_VCF,
        "--sample_info", SAMPLE_INFO,
        "--region",      REGION,
        "--out",         str(out_vcf),
        "--pop_field",          # annotate ancestry in VCF FORMAT field
    ]

    print(f"\n  Running {scenario_name}...")
    print(f"    cmd: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR in {scenario_name}:")
        print(result.stderr)
    else:
        bp_file = out_prefix.parent / (out_prefix.stem + ".bp")
        print(f"    genotypes  → {out_vcf}")
        print(f"    breakpoints→ {bp_file}")
        if result.stderr:
            # haptools logs to stderr by default; show last few lines
            log_lines = result.stderr.strip().splitlines()
            for line in log_lines[-5:]:
                print(f"    [log] {line}")


# ─────────────────────────────────────────────────────────────────────────────
# Parse breakpoints into a readable summary (optional inspection utility)
# ─────────────────────────────────────────────────────────────────────────────

def summarize_breakpoints(bp_file: Path):
    """
    Print a human-readable summary of the .bp breakpoints file.

    .bp format (tab-delimited, per haplotype per individual):
        sample_id  haplotype_index  population  chrom  start_bp  end_bp

    This is your ground-truth local ancestry — use it to evaluate your HMM.
    """
    if not bp_file.exists():
        print(f"  [skip] {bp_file} not found")
        return

    print(f"\n  Breakpoints summary: {bp_file.name}")
    samples_seen = set()
    n_segments = 0
    pop_counts: dict = {}

    with open(bp_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            sample, hap_idx, pop, chrom = parts[0], parts[1], parts[2], parts[3]
            samples_seen.add(sample)
            n_segments += 1
            pop_counts[pop] = pop_counts.get(pop, 0) + 1

    print(f"    individuals : {len(samples_seen)}")
    print(f"    total segments : {n_segments}")
    print(f"    ancestry segment counts: {pop_counts}")


if __name__ == "__main__":

    n_haplotypes = N_SAMPLES * 2  # diploid

    OUT_DIR.mkdir(exist_ok=True)
    model_dir = OUT_DIR / "models"

    print("\nStep 1: Writing model files...")
    build_model_files(model_dir, n_haplotypes)

    print("\nStep 2: Running haptools simgenotype for each scenario...")

    scenarios = [
        ("S1_trivial",          model_dir / "S1_trivial.dat"),
        ("S2_two_ancestries",   model_dir / "S2_two_ancestries.dat"),
        ("S3_three_ancestries", model_dir / "S3_three_ancestries.dat"),
    ]

    for name, model_file in scenarios:
        out_prefix = OUT_DIR / name
        run_simgenotype(name, model_file, out_prefix)

    print("\nStep 3: Breakpoint summaries (ground truth inspection)...")
    for name, _ in scenarios:
        bp_file = OUT_DIR / f"{name}.bp"
        summarize_breakpoints(bp_file)

    print(f"""
Done. Outputs in ./{OUT_DIR}/

For each scenario you now have:
  <name>.vcf.gz   — phased genotypes to feed into your HMM as test input
  <name>.bp       — ground-truth local ancestry breakpoints for evaluation

To load breakpoints in Python for evaluation against your HMM output:

    import pandas as pd
    bp = pd.read_csv("toy_examples_haptools/S3_three_ancestries.bp",
                     sep="\\t", comment="#",
                     names=["sample","hap","population","chrom","start","end"])

To visualise ancestry tracts (optional):
    haptools karyogram toy_examples_haptools/S3_three_ancestries.bp \\
        --output S3_karyogram.png
""")
