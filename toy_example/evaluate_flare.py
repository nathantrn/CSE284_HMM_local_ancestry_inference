"""
evaluate_flare.py
Compares FLARE ancestry calls against haptools ground truth breakpoints
at the population level (YRI, IBS, PEL) — no superpopulation grouping.
"""

import re
import gzip
import numpy as np
from cyvcf2 import VCF
from pathlib import Path

FLARE_DIR = Path("flare_output")
BP_DIR    = Path("toy_examples_haptools")

SCENARIOS = [
    "S1_trivial",
    "S2_two_ancestries",
    "S3_three_ancestries",
]

# Ground truth populations used in each scenario
SCENARIO_POPS = {
    "S1_trivial":          ["YRI"],
    "S2_two_ancestries":   ["YRI", "IBS"],
    "S3_three_ancestries": ["YRI", "IBS", "PEL"],
}

# ─────────────────────────────────────────────────────────────────────────────

def parse_ancestry_header(vcf_path: Path) -> dict:
    """Return {population: integer_index} from ##ANCESTRY header line."""
    ancestry = {}
    with gzip.open(vcf_path, "rt", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("##ANCESTRY"):
                if line.startswith("#CHROM"):
                    break
                continue
            for match in re.finditer(r"\b([A-Z][A-Z0-9]+)=(\d+)\b", line):
                pop, idx = match.group(1), int(match.group(2))
                if pop not in ancestry:
                    ancestry[pop] = idx
    return ancestry


def parse_bp(bp_file: Path) -> dict:
    """
    Returns {sample: {"hap1": [(end_pos, pop), ...], "hap2": [...]}}
    """
    result = {}
    current_key = None
    with open(bp_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                sample = "_".join(parts[0].split("_")[:2])
                hap    = "hap1" if parts[0].endswith("_1") else "hap2"
                result.setdefault(sample, {"hap1": [], "hap2": []})
                current_key = (sample, hap)
            elif len(parts) == 4 and current_key:
                pop, _chrom, end_pos, _ = parts
                result[current_key[0]][current_key[1]].append((int(end_pos), pop))
    return result


def bp_to_per_snp(segments: list, positions: np.ndarray) -> list:
    """Convert segment list [(end_pos, pop)] to per-SNP population labels."""
    calls = []
    for pos in positions:
        for end_pos, pop in segments:
            if pos <= end_pos:
                calls.append(pop)
                break
        else:
            calls.append(segments[-1][1])
    return calls


# ─────────────────────────────────────────────────────────────────────────────

def evaluate(scenario: str):
    print(f"\n{'='*60}")
    print(f"Scenario: {scenario}")
    print(f"{'='*60}")

    vcf_path = FLARE_DIR / scenario / f"{scenario}.anc.vcf.gz"
    bp_path  = BP_DIR / f"{scenario}.bp"

    # Parse ancestry index
    ancestry_idx = parse_ancestry_header(vcf_path)
    idx_to_pop   = {v: k for k, v in ancestry_idx.items()}

    # Read FLARE calls
    vcf       = VCF(str(vcf_path))
    samples   = vcf.samples
    positions = []
    an1_calls = {s: [] for s in samples}
    an2_calls = {s: [] for s in samples}

    for variant in vcf:
        positions.append(variant.POS)
        an1 = variant.format("AN1").flatten()
        an2 = variant.format("AN2").flatten()
        for i, s in enumerate(samples):
            an1_calls[s].append(idx_to_pop.get(int(an1[i]), "?"))
            an2_calls[s].append(idx_to_pop.get(int(an2[i]), "?"))

    positions = np.array(positions)
    bp_data   = parse_bp(bp_path)
    pops      = SCENARIO_POPS[scenario]

    for sample in samples:
        if sample not in bp_data:
            print(f"  [skip] {sample} not in .bp file")
            continue

        print(f"\n  {sample}")
        for hap, flare_calls in [("hap1", an1_calls[sample]),
                                  ("hap2", an2_calls[sample])]:
            truth = bp_to_per_snp(bp_data[sample][hap], positions)

            correct = sum(t == f for t, f in zip(truth, flare_calls))
            total   = len(truth)
            print(f"    {hap}: {correct}/{total} = {correct/total*100:.1f}%")

            # Per-population breakdown with misclassification info
            for pop in pops:
                idx = [i for i, t in enumerate(truth) if t == pop]
                if not idx:
                    continue
                pop_correct = sum(truth[i] == flare_calls[i] for i in idx)
                wrong = [flare_calls[i] for i in idx if truth[i] != flare_calls[i]]
                wrong_counts = {p: wrong.count(p) for p in set(wrong)}
                print(f"      {pop}: {pop_correct}/{len(idx)} = {pop_correct/len(idx)*100:.1f}%", end="")
                if wrong_counts:
                    print(f"  [misclassified as: {wrong_counts}]", end="")
                print()


if __name__ == "__main__":
    for scenario in SCENARIOS:
        evaluate(scenario)