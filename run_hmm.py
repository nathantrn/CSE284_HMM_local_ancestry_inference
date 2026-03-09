import numpy as np
import pandas as pd
import allel
import sys
import logging
import argparse

"""
initial state probabilities: pi(i,h) = mu_i * q_ih = mu_i * p_ij / n_j
where mu_i and p_ij come from model parameters, n_j is number of haplotypes in reference panel j

transition probabilities: P(S_m = (i', h') | S_m-1 = (i, h))
d_m is the distance between markers m and m-1 in morgans (can be taken from hapmap genetic map)

if i = i' and h = h' (no change in ancestry or ref haplotype):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h' + e^(-d_m * T) * (1 - e^(-d_m * rho_i)) * q_i'h' + e^(-d_m * T) * e^(-d_m * rho_i)

if i = i' and h ≠ h' (no change in ancestry but ref haplotype changes):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h' + e^(-d_m * T) * (1 - e^(-d_m * rho_i)) * q_i'h'

if i ≠ i' (ancestry changes, and thus ref haplotype must also change):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h'

emission probabilities: e_m(i, h)
Let I_m(h) = 1 if the allele at marker m on haplotype h matches the observed allele on the admixed haplotype, and let I_m(h) = 0 otherwise
e_m(i, h) = (theta_ij)^(1 - I_m(h)) * (1 - theta_ij)^I_m(h)

We calculate the local ancestry proportions for each of a sample's two haplotypes separately.
For each position m, we will have hidden states S_m(i,h), where i is the ancestry at the target position and h is the reference haplotype.
First, calculate emission matrix for each reference haplotype and marker possibility, comparing the ref to the target at that marker to get I_m(h).
Then, calculate transition matrix starting from the second marker (since the first will come from initial state probabilities).
Can only transition between one marker to the consecutive one, but must consider switching between all ref haplotype possibilities (N^2 transitions per marker 2 and beyond).
"""

def processModelFile(file_path):
    mapping = {
        "list of ancestries": "ancestries",
        "list of reference panels": "panels",
        "T:": "T",
        "mu[i]": "mu",
        "p[i][j]": "p_ij",
        "theta[i][j]": "theta_ij",
        "rho[i]": "rho"
    }
    
    results = {k: [] for k in mapping.values()}
    current_key = None

    with open(file_path, 'r') as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line:
                continue
            
            if clean_line.startswith("#"):
                for keyword, key in mapping.items():
                    if keyword in clean_line:
                        current_key = key
                        break
                continue
            
            if current_key:
                parts = clean_line.split()
                try:
                    row = [float(x) for x in parts]
                except ValueError:
                    row = parts
                
                results[current_key].append(row)
    
    ancestries = results['ancestries'][0]
    panels = results['panels'][0]
    T = results['T'][0][0]
    mu = np.array(results['mu'][0])
    rho = np.array(results['rho'][0])
    
    p_ij = np.array(results['p_ij'])
    theta_ij = np.array(results['theta_ij'])

    df_p = pd.DataFrame(p_ij, index=ancestries, columns=panels)
    df_theta = pd.DataFrame(theta_ij, index=ancestries, columns=panels)

    return {
        "ancestries": ancestries,
        "panels": panels,
        "T": T,
        "mu": mu,
        "rho": rho,
        "p": df_p,
        "theta": df_theta
    }

def extractHaplotypes(vcf_path):
    fields = allel.read_vcf(vcf_path, fields=['variants/CHROM', 'variants/POS', 'calldata/GT', 'samples'])
    
    chrom = np.unique(fields['variants/CHROM'])[0]
    gt = fields['calldata/GT']
    samples = fields['samples']
    positions = fields['variants/POS']

    haplotypes = gt.reshape(gt.shape[0], -1)
    
    hap_names = []
    for s in samples:
        hap_names.extend([f"{s}_h1", f"{s}_h2"])
    
    hap_df = pd.DataFrame(haplotypes, columns = hap_names, index = positions)
        
    return hap_df, chrom

# pi(i,h) = mu_i * q_ih = mu_i * p_ij / n_j
def calculateInitialProbs(model_params, panel_map, ref_haps, q):
    pi = []
    for i in range(len(model_params['ancestries'])):
        for h in range(len(ref_haps.columns)):
            curr_panel = panel_map[ref_haps.columns[h]]
            q_ih = q.loc[model_params['ancestries'][i], curr_panel] 
            pi.append(model_params['mu'][i] * q_ih)
    return np.array(pi)

# e_m(i, h) = (theta_ij)^(1 - I_m(h)) * (1 - theta_ij)^I_m(h)
def calculateEmissionProbs(theta, ref_haplotypes, sample_haplotype, panel_map, marker):
    e_m = []
    for i in range(len(theta.index)):
        for h in range(len(ref_haplotypes.columns)):
            curr_panel = panel_map[ref_haplotypes.columns[h]]
            I_m = ref_haplotypes.loc[marker, ref_haplotypes.columns[h]] == sample_haplotype.loc[marker]
            e_m.append((theta.loc[theta.index[i], curr_panel] ** (1 - I_m)) * ((1 - theta.loc[theta.index[i], curr_panel]) ** I_m))
    return np.array(e_m)

"""
if i = i' and h = h' (no change in ancestry or ref haplotype):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h' + e^(-d_m * T) * (1 - e^(-d_m * rho_i)) * q_i'h' + e^(-d_m * T) * e^(-d_m * rho_i)

if i = i' and h ≠ h' (no change in ancestry but ref haplotype changes):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h' + e^(-d_m * T) * (1 - e^(-d_m * rho_i)) * q_i'h'

if i ≠ i' (ancestry changes, and thus ref haplotype must also change):
P = (1 - e^(-d_m * T)) * mu_i' * q_i'h'
"""
def calculateTransitionProbs(model_params, morgans_per_bp, panel_map, prev_marker, curr_marker, ancestry_states, haplotype_states, q):
    t_m = np.empty((len(ancestry_states), len(ancestry_states)))
    d_m = morgans_per_bp * (curr_marker - prev_marker)
    T_term = np.exp(-d_m * model_params['T'])
    rho_terms = np.exp(-d_m * model_params['rho'])
    for y in range(len(ancestry_states)):
        i_curr = ancestry_states[y]
        i_curr_idx = model_params['ancestries'].index(i_curr)
        h_curr = haplotype_states[y]
        curr_panel = panel_map[h_curr]
        q_ih = q.loc[i_curr, curr_panel] 
        rho_term = rho_terms[i_curr_idx]
        part1 = (1 - T_term) * model_params['mu'][i_curr_idx] * q_ih
        for x in range(len(ancestry_states)):
            i_prev = ancestry_states[x]
            h_prev = haplotype_states[x]
            if i_prev == i_curr:
                part2 = T_term * (1 - rho_term) * q_ih
                if h_prev == h_curr:
                    part3 = T_term * rho_term
                    t_m[x,y] = part1 + part2 + part3
                else:
                    t_m[x,y] = part1 + part2
            else:
                t_m[x,y] = part1
    return t_m

def globalAncestry(model_params, morgans_per_bp, ref_haps, panel_map, sample_haplotype, ancestry_states, haplotype_states, q):
    forward = np.zeros((len(ref_haps.index), len(ancestry_states)))
    backward = np.zeros((len(ref_haps.index), len(ancestry_states)))
    pi = calculateInitialProbs(model_params, panel_map, ref_haps, q)
    e_0 = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, panel_map, ref_haps.index[0])
    
    forward[0] = pi * e_0
    c_0 = forward[0].sum()
    if c_0 > 0:
        forward[0] = forward[0] / c_0

    for m in range(1, len(ref_haps.index)):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, panel_map, ref_haps.index[m-1], ref_haps.index[m], ancestry_states, haplotype_states, q)
        e_m = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, panel_map, ref_haps.index[m])

        forward[m] = e_m * (trans.T @ forward[m - 1])

        c_m = forward[m].sum()
        if c_m > 0:
            forward[m] = forward[m] / c_m
    logging.info(f'forward pass complete for {sample_haplotype.name}')
    
    backward[-1] = np.ones(len(ancestry_states))
        
    for m in range(len(ref_haps.index) - 2, -1, -1):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, panel_map, ref_haps.index[m], ref_haps.index[m + 1], ancestry_states, haplotype_states, q)
        e_next = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, panel_map, ref_haps.index[m + 1])

        backward[m] = trans @ (e_next * backward[m + 1])

        b_sum = backward[m].sum()
        if b_sum > 0:
            backward[m] = backward[m] / b_sum

    logging.info(f'backward pass complete for {sample_haplotype.name}')

    posterior = pd.DataFrame(forward * backward)
    posterior = posterior.div(posterior.sum(axis=1), axis=0) # normalize
            
    ancestry_counts = np.zeros(len(model_params['ancestries']))
    for m in range(len(ref_haps.index)):
        for i in range(len(ancestry_states)):
            curr_ancestry = ancestry_states[i]
            curr_ancestry_idx = model_params['ancestries'].index(curr_ancestry)
            ancestry_counts[curr_ancestry_idx] += posterior.iloc[m, i]

    global_ancestry = ancestry_counts / ancestry_counts.sum()
    logging.info(f'global ancestry calculated for {sample_haplotype.name}')
            
    return global_ancestry
        
def viterbi(model_params, morgans_per_bp, ref_haps, panel_map, sample_haplotype, ancestry_states, haplotype_states, q):
    viterbi = np.zeros((len(ref_haps.index), len(ancestry_states)))
    backpointer = np.zeros((len(ref_haps.index), len(ancestry_states)), dtype=int)

    pi = calculateInitialProbs(model_params, panel_map, ref_haps, q)
    e_0 = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, panel_map, ref_haps.index[0])
    viterbi[0] = np.log(pi + 1e-10) + np.log(e_0 + 1e-10)

    for m in range(1, len(ref_haps.index)):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, panel_map, ref_haps.index[m-1], ref_haps.index[m], ancestry_states, haplotype_states, q)
        e_m = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, panel_map, ref_haps.index[m])

        for curr_state_idx in range(len(ancestry_states)):
            prev_scores = viterbi[m - 1] + np.log(trans[:, curr_state_idx] + 1e-10)
            backpointer[m, curr_state_idx] = np.argmax(prev_scores)
            viterbi[m, curr_state_idx] = (np.max(prev_scores) + np.log(e_m[curr_state_idx] + 1e-10))
        if m % 100 == 99:
            logging.info(f'viterbi processed marker {m + 1} for {sample_haplotype.name}')

    path_indices = [np.argmax(viterbi[-1])]

    for m in range(len(ref_haps.index) - 2, -1, -1):
        path_indices.append(backpointer[m + 1, path_indices[-1]])

    path_indices.reverse()
    path = [model_params['ancestries'].index(ancestry_states[idx]) for idx in path_indices]

    return path

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('ref_vcf_path', help='path to phased VCF with reference samples')
    parser.add_argument('ref_panel_path', help='path to tab-separated reference panel')
    parser.add_argument('admixed_vcf_path', help='path to phased VCF with admixed samples')
    parser.add_argument('genetic_map_path', help='path to genetic map in cM')
    parser.add_argument('model_params_path', help='path to model parameters generated by FLARE')
    parser.add_argument('-o','--output', help='prefix for output files')
    parser.add_argument('-g','--global_off', action='store_false',help='exclude global ancestry calculation output')
    parser.add_argument('-l','--local_off', action='store_false',help='exclude global ancestry calculation output')
    args = parser.parse_args()
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=f'{args.output}.log',
        filemode='w'
    )

    logging.info('starting pipeline')
    model_params = processModelFile(args.model_params_path)
    ref_haps, chrom = extractHaplotypes(args.ref_vcf_path)
    admixed_haps, chrom = extractHaplotypes(args.admixed_vcf_path)
    ref_panel = pd.read_csv(args.ref_panel_path, sep='\t', header=None)
    with open(args.genetic_map_path, 'r') as f:
        lines = f.readlines()
        genetic_map = lines[-1].strip().split()
    morgans_per_bp = float(genetic_map[2]) / float(genetic_map[3]) / 100
    panel_map = {ref_haps.columns[h]: ref_panel.loc[ref_panel.iloc[:,0] == ref_haps.columns[h][:-3], 1].values[0] for h in range(len(ref_haps.columns))}
    ancestry_states = [item for x in [[model_params['ancestries'][i]] * len(ref_haps.columns) for i in range(len(model_params['ancestries']))] for item in x]
    haplotype_states = list(ref_haps.columns) * len(model_params['ancestries'])
    n = pd.Series(list(panel_map.values())).value_counts()
    q = model_params['p'].div(n, axis=1)
    if args.global_off: global_ancestry_df = pd.DataFrame(index = admixed_haps.columns, columns = model_params['ancestries'])
    if args.local_off: ancestry_haplotype_df = pd.DataFrame(index = admixed_haps.index, columns = admixed_haps.columns)
    logging.info('inputs processed')
    for i in range(len(admixed_haps.columns)):
        if args.global_off: global_ancestry_df.iloc[i] = globalAncestry(model_params, morgans_per_bp, ref_haps, panel_map, admixed_haps.iloc[:,i], ancestry_states, haplotype_states, q)
        if args.local_off: ancestry_haplotype_df.iloc[:,i] = viterbi(model_params, morgans_per_bp, ref_haps, panel_map, admixed_haps.iloc[:,i], ancestry_states, haplotype_states, q)
        logging.info(f'haplotype {i + 1} ({admixed_haps.iloc[:,i].name}) processed')
    logging.info('all haplotypes processed')
    if args.global_off:
        global_ancestry_averaged_df = global_ancestry_df.groupby(np.arange(len(global_ancestry_df)) // 2).mean()
        global_ancestry_averaged_df.index = [x[:-3] for x in global_ancestry_df.index[np.arange(int(len(global_ancestry_df) / 2)) * 2]]
        global_ancestry_averaged_df.to_csv(f'{args.output}_global_ancestry.txt', sep='\t')
    if args.local_off:
        ancestry_haplotype_df = ancestry_haplotype_df.reset_index()
        ancestry_haplotype_df.insert(loc = 0, column = 'CHROM', value = chrom)
        ancestry_haplotype_df = ancestry_haplotype_df.rename(columns={'index': 'POS'})
        ancestry_haplotype_df.to_csv(f'{args.output}_ancestry_haplotypes.txt', sep='\t', index=False)
    logging.info('results saved')


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("wrong number of arguments! please run the the command in the following format: python run_hmm.py <ref VCF> <ref panel> <admixed VCF> <genetic map> <model params> [-o -g -l]")
    
    else:
        main()
