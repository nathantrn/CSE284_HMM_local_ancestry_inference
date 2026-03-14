"""
HaploHMM gives out outputs slightly differently than FLARE, such as consolidating both samples into one output file instead of spreading them out.

evaluate_haplohmm.py
Compares haploHMM ancestry calls against haptools ground truth breakpoints
at the population level (YRI, IBS, PEL)

haploHMM output format:
    Single TSV: toy_example_ancestry_haplotypes.txt
    Columns: CHROM, POS, {scenario}_{sample_num}_h1, {scenario}_{sample_num}_h2, ...
    Values: integer ancestry index (1-based)

Ancestry index → population mapping :
    1 = IBS
    2 = YRI
    3 = PEL
"""

import numpy as np
import pandas as pd
from pathlib import Path

HAPLOHMM_DIR = Path("haploHMM_output")
BP_DIR       = Path("toy_examples_haptools")

SCENARIOS = [
    "S1_trivial",
    "S2_two_ancestries",
    "S3_three_ancestries",
]

SCENARIO_POPS = {
    "S1_trivial":          ["YRI"],
    "S2_two_ancestries":   ["YRI", "IBS"],
    "S3_three_ancestries": ["YRI", "IBS", "PEL"],
}

ANCESTRY_MAP = {
    1: "IBS",
    2: "YRI", 
    3: "PEL",
}


def parse_bp(bp_file: Path) -> dict:
    result = {}
    current_key = None
    with open(bp_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                tokens  = parts[0].split("_")  
                num     = int(tokens[1])
                hap_idx = tokens[2]
                hap     = "hap1" if hap_idx == "1" else "hap2"
                result.setdefault(num, {"hap1": [], "hap2": []})
                current_key = (num, hap)
            elif len(parts) == 4 and current_key:
                pop, _chrom, end_pos, _ = parts
                result[current_key[0]][current_key[1]].append((int(end_pos), pop))
    return result


def bp_to_per_snp(segments: list, positions: np.ndarray) -> list:
    calls = []
    for pos in positions:
        for end_pos, pop in segments:
            if pos <= end_pos:
                calls.append(pop)
                break
        else:
            calls.append(segments[-1][1])
    return calls


def load_haplohmm(txt_path: Path) -> tuple[np.ndarray, pd.DataFrame]:
    raw = pd.read_csv(txt_path, sep="\t")
    positions = raw["POS"].values
    hap_cols = [c for c in raw.columns if c not in ("CHROM", "POS")]
    df = raw[hap_cols].copy()
    for col in hap_cols:
        df[col] = df[col].map(ANCESTRY_MAP).fillna("?")
    return positions, df


def evaluate(scenario: str, positions: np.ndarray, df: pd.DataFrame,
             bp_data: dict):
    print(f"\n{'='*60}")
    print(f"Scenario: {scenario}")
    print(f"{'='*60}")

    pops = SCENARIO_POPS[scenario]

    scenario_cols = [c for c in df.columns if c.startswith(scenario + "_")]

    samples = sorted(set(
        c.rsplit("_h", 1)[0] for c in scenario_cols
    ))

    for i, col_sample in enumerate(samples, start=1):
        bp_sample = i

        if bp_sample not in bp_data:
            print(f"  [skip] {col_sample} (Sample_{i}) not found in .bp file")
            continue

        print(f"\n  {col_sample}")

        for hap_suffix, hap_key in [("h1", "hap1"), ("h2", "hap2")]:
            col_name = f"{col_sample}_{hap_suffix}"
            if col_name not in df.columns:
                print(f"    {hap_key}: column '{col_name}' not found — skipping")
                continue

            flare_calls = df[col_name].tolist()
            truth       = bp_to_per_snp(bp_data[bp_sample][hap_key], positions)

            correct = sum(t == f for t, f in zip(truth, flare_calls))
            total   = len(truth)
            print(f"    {hap_key}: {correct}/{total} = {correct/total*100:.1f}%")
            for pop in pops:
                idx = [i for i, t in enumerate(truth) if t == pop]
                if not idx:
                    continue
                pop_correct = sum(truth[i] == flare_calls[i] for i in idx)
                wrong = [flare_calls[i] for i in idx
                         if truth[i] != flare_calls[i]]
                wrong_counts = {p: wrong.count(p) for p in set(wrong)}
                print(f"      {pop}: {pop_correct}/{len(idx)} = "
                      f"{pop_correct/len(idx)*100:.1f}%", end="")
                if wrong_counts:
                    print(f"  [misclassified as: {wrong_counts}]", end="")
                print()

if __name__ == "__main__":
    txt_path = HAPLOHMM_DIR / "toy_example_ancestry_haplotypes.txt"
    print(f"Loading haploHMM calls from: {txt_path}")
    positions, df = load_haplohmm(txt_path)
    print(f"  {len(positions)} SNPs, {len(df.columns)} haplotype columns")

    for scenario in SCENARIOS:
        bp_path  = BP_DIR / f"{scenario}.bp"
        bp_data  = parse_bp(bp_path)
        evaluate(scenario, positions, df, bp_data)
