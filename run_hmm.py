import numpy as np
import pandas as pd
import allel
import sys

"""
NOTES: 
- do we want to use 1000 genomes for both train and test?
- how to download hgdp
"""

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
For each sample, we will have an N possible hidden states, i.e. the total number of reference haplotypes.
So, we end up with N x M nodes in the HMM, where M is the total number of markers/variants represented in the VCF files.
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

    # --- Post-Processing into Clean Arrays/Scalars ---
    
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
def calculateInitialProbs(model_params, ref_panel, ref_haps):
    pi = []
    n = ref_panel.iloc[:,1].value_counts()
    for i in range(len(model_params['ancestries'])):
        for h in range(len(ref_haps.columns)):
            curr_panel = ref_panel.iloc[int(h/2),1]
            p_ij = model_params['p'].loc[model_params['ancestries'][i], curr_panel] 
            n_j = n[curr_panel]
            pi.append(model_params['mu'][i] * p_ij / n_j)
    return np.array(pi)

# e_m(i, h) = (theta_ij)^(1 - I_m(h)) * (1 - theta_ij)^I_m(h)
def calculateEmissionProbs(theta, ref_haplotypes, sample_haplotype, ref_panel, marker):
    e_m = []
    for i in range(len(theta.index)):
        for h in range(len(ref_haplotypes.columns)):
            curr_panel = ref_panel[ref_panel.iloc[:,0] == ref_haplotypes.columns[h][:-3]].iloc[:,1].values[0]
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
def calculateTransitionProbs(model_params, morgans_per_bp, ref_haps, ref_panel, prev_marker, curr_marker, states):
    t_m = np.empty((len(states), len(states)))
    n = ref_panel.iloc[:,1].value_counts()
    d_m = morgans_per_bp * (curr_marker - prev_marker)
    T_term = np.exp(-d_m * model_params['T'])
    for y in range(len(states)):
        i_curr = states[y].split(',')[0]
        i_curr_idx = model_params['ancestries'].index(i_curr)
        h_curr = states[y].split(',')[1]
        h_curr_idx = list(ref_haps.columns).index(h_curr)
        curr_panel = ref_panel.iloc[int(h_curr_idx/2),1]
        p_ij = model_params['p'].loc[model_params['ancestries'][i_curr_idx], curr_panel] 
        n_j = n[curr_panel]
        q_ih = p_ij / n_j
        rho_term = np.exp(-d_m * model_params['rho'][i_curr_idx])
        part1 = (1 - T_term) * model_params['mu'][i_curr_idx] * q_ih
        part2 = T_term * (1 - rho_term) * q_ih
        part3 = T_term * rho_term
        for x in range(len(states)):
            i_prev = states[x].split(',')[0]
            h_prev = states[x].split(',')[1]
            if i_prev == i_curr:
                if h_prev == h_curr:
                    t_m[x,y] = part1 + part2 + part3
                else:
                    t_m[x,y] = part1 + part2
            else:
                t_m[x,y] = part1
    return t_m

def globalAncestry(model_params, morgans_per_bp, ref_haps, ref_panel, sample_haplotype, states):
    forward = np.zeros((len(ref_haps.index), len(states)))
    backward = np.zeros((len(ref_haps.index), len(states)))
    log_likelihood = 0
    pi = calculateInitialProbs(model_params, ref_panel, ref_haps)
    e_0 = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, ref_panel, ref_haps.index[0])
    
    forward[0] = pi * e_0
    c_0 = forward[0].sum()
    if c_0 > 0:
        forward[0] = forward[0] / c_0
        log_likelihood += np.log(c_0)

    for m in range(1, len(ref_haps.index)):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, ref_haps, ref_panel, ref_haps.index[m-1], ref_haps.index[m], states)
        e_m = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, ref_panel, ref_haps.index[m])

        forward[m] = e_m * (trans.T @ forward[m - 1])

        c_m = forward[m].sum()
        if c_m > 0:
            forward[m] = forward[m] / c_m
            log_likelihood += np.log(c_m)
    
    backward[-1] = np.ones(len(states))
        
    for m in range(len(ref_haps.index) - 2, -1, -1):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, ref_haps, ref_panel, ref_haps.index[m], ref_haps.index[m + 1], states)
        e_next = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, ref_panel, ref_haps.index[m + 1])

        backward[m] = trans @ (e_next * backward[m + 1])

        b_sum = backward[m].sum()
        if b_sum > 0:
            backward[m] = backward[m] / b_sum

    posterior = forward * backward
    for m in range(len(ref_haps.index)): # normalize
        row_sum = posterior[m].sum()
        if row_sum > 0:
            posterior[m] = posterior[m] / row_sum
            
    ancestry_counts = np.zeros(len(model_params['ancestries']))
    for m in range(len(ref_haps.index)):
        for i in range(len(states)):
            curr_ancestry = states[i].split(',')[0]
            curr_ancestry_idx = model_params['ancestries'].index(curr_ancestry)
            ancestry_counts[curr_ancestry_idx] += posterior[m, i]

    global_ancestry = ancestry_counts / ancestry_counts.sum()
            
    return global_ancestry
        
def viterbi(model_params, morgans_per_bp, ref_haps, ref_panel, sample_haplotype, states):
    viterbi = np.zeros((len(ref_haps.index), len(states)))
    backpointer = np.zeros((len(ref_haps.index), len(states)), dtype=int)

    pi = calculateInitialProbs(model_params, ref_panel, ref_haps)
    e_0 = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, ref_panel, ref_haps.index[0])
    viterbi[0] = np.log(pi + 1e-10) + np.log(e_0 + 1e-10)

    for m in range(1, len(ref_haps.index)):
        trans = calculateTransitionProbs(model_params, morgans_per_bp, ref_haps, ref_panel, ref_haps.index[m-1], ref_haps.index[m], states)
        e_m = calculateEmissionProbs(model_params['theta'], ref_haps, sample_haplotype, ref_panel, ref_haps.index[m])

        for curr_state_idx in range(len(states)):
            prev_scores = viterbi[m - 1] + np.log(trans[:, curr_state_idx] + 1e-10)
            backpointer[m, curr_state_idx] = np.argmax(prev_scores)
            viterbi[m, curr_state_idx] = (np.max(prev_scores) + np.log(e_m[curr_state_idx] + 1e-10))

    path_indices = [np.argmax(viterbi[-1])]
    # log_prob = viterbi[-1, path_indices[0]]

    for m in range(len(ref_haps.index) - 2, -1, -1):
        path_indices.append(backpointer[m + 1, path_indices[-1]])

    path_indices.reverse()
    path = [states[idx] for idx in path_indices]

    return path

def main(ref_vcf_path, ref_panel_path, admixed_vcf_path, genetic_map_path, model_params_path):
    model_params = processModelFile(model_params_path)
    ref_haps, chrom = extractHaplotypes(ref_vcf_path)
    admixed_haps, chrom = extractHaplotypes(admixed_vcf_path)
    ref_panel = pd.read_csv(ref_panel_path, sep='\t', header=None)
    with open(genetic_map_path, 'r') as f:
        lines = f.readlines()
        genetic_map = lines[-1].strip().split()
    morgans_per_bp = float(genetic_map[2]) / float(genetic_map[3]) / 100
    states = [item for x in [list(model_params['ancestries'][i] + ',' + ref_haps.columns) for i in range(len(model_params['ancestries']))] for item in x]
    global_ancestry_df = pd.DataFrame(index = admixed_haps.columns, columns = model_params['ancestries'])
    ancestry_haplotype_df = pd.DataFrame(index = admixed_haps.index, columns = admixed_haps.columns)
    print('inputs processed')
    for i in range(len(admixed_haps.columns)):
        global_ancestry_df.iloc[i] = globalAncestry(model_params, morgans_per_bp, ref_haps, ref_panel, admixed_haps.iloc[:,i], states)
        ancestry_haplotype_df.iloc[:,i] = viterbi(model_params, morgans_per_bp, ref_haps, ref_panel, admixed_haps.iloc[:,i], states)
        if i % 100 == 0:
            print(f'haplotype {i + 1} processed')
    print('all haplotypes processed')
    global_ancestry_averaged_df = global_ancestry_df.groupby(global_ancestry_df.index // 2).mean()
    global_ancestry_averaged_df.to_csv('global_ancestry_results.txt', sep='\t')
    ancestry_haplotype_df = ancestry_haplotype_df.reset_index()
    ancestry_haplotype_df.insert(loc = 0, column = 'CHROM', value = chrom)
    ancestry_haplotype_df.to_csv('ancestry_haplotypes_results.txt', sep='\t')


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("wrong number of arguments! please run the the command in the following format: python train_hmm.py <ref VCF> <ref panel> <admixed VCF> <genetic map> <model params>")
    
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
