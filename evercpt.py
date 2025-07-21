import re
import subprocess
import tempfile
import os
import joblib
import numpy as np
import pandas as pd
from collections import deque, Counter


def is_valid_sequence(seq):
    seq = seq.upper()
    if len(seq) >= 10000 or len(seq) < 10:
        return False
    if re.fullmatch(r'[ACGTU]+', seq):
        return True
    return False

def save_fasta(seq, filepath):
    with open(filepath, 'w') as f:
        f.write(f">input\n{seq}\n")

def run_rnafold(fasta_file, output, circ=False):
    # subprocess.run(f"RNAfold --noLP --noPS {fasta_file} > {output}", shell=True)
    with open(os.path.join(output, 'seq.dbn'), 'w') as f:
        if circ:
            subprocess.run(['RNAfold', '-c', '--noLP', '--noPS', fasta_file], stdout=f, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(['RNAfold', '--noLP', '--noPS', fasta_file], stdout=f, stderr=subprocess.DEVNULL)
    with open(os.path.join(output, 'seq.dbn'), 'r') as f:
        lines = f.read().strip().split('\n')
        structure_line = lines[2] if len(lines) > 2 else ""
        structure = structure_line.split(' ')[0]
        mfe = float(re.findall(r'-?\d+\.\d+', structure_line)[-1])
    return structure, mfe

def run_bprna(path):
    try:
        subprocess.run(['perl', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bpRNA/bpRNA.pl'), os.path.join(path, 'seq.dbn')], 
                       cwd=path, stderr=subprocess.PIPE, text=True, check=True)
        # extract number of stem loops and multiloops from the output file
        with open(os.path.join(path, 'seq.st'), 'r') as f:
            lines = f.readlines()
        sl = 0
        ml = set()
        for line in lines:
            line = line.strip()
            if re.match(r'^S\d+', line):  # Hairpin loop (stem-loop)
                sl += 1
            elif re.match(r'^M\d+', line):  # Multiloop
                ml.add(line.split()[0].split('.')[0])
        return sl, len(ml)
    except subprocess.CalledProcessError as e:
        if "line 1558" in e.stderr:
            return 0, 0

def dot2bp(dot):
    return (len(dot) - dot.count("."))

def dot2pairs(dot):
    stack, pairs = deque(), []
    for i, n in enumerate(dot):
        if n == '(':
            stack.append(i)
        elif n == ')':
            pairs.append((stack.pop(), i))
    return np.array(pairs)

def dot2bp90(pairs):
    distance = []
    if pairs.size > 0:
        for i in pairs:
            distance.append(abs(i[1]-i[0]))
        return np.percentile(distance, 90)
    else:
        return 0

def dot2MLD(dot, pairs0):
    if len(pairs0) == 0:
        return 0

    MLD = 0
    shifts = [0]
    length = len(dot)
    seq = np.zeros(length, dtype=int)
    curr_mount = None

    while shifts:
        # Shift the sequence
        pairs = (pairs0 - shifts.pop()) % length

        # Find positions of opening and closing nucleotides
        opens, closes = np.min(pairs, axis=1), np.max(pairs, axis=1)

        seq.fill(0)
        seq[opens] = 1
        seq[closes] = -1

        mountain = np.cumsum(seq)

        MLD = max(MLD, int(np.max(mountain)))

        # In first iteration, calculate all required shifts
        if curr_mount is None:
            curr_mount = mountain[0]
            for i, m in enumerate(mountain):
                if m == curr_mount: 
                    shifts.append(i)
                else:
                    curr_mount = m
    return MLD

def structure_to_features(structure):
    pairs = dot2pairs(structure)
    return dot2bp(structure), dot2bp90(pairs), dot2MLD(structure, pairs)

def run_cpc2(fasta_file, output):
    subprocess.run(['python', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CPC2_standalone-1.0.1/bin/CPC2.py'), 
                              '-i', fasta_file, '-o', os.path.join(output, 'cpc')], cwd=output, check=True)
    with open(os.path.join(output, 'cpc.txt'), 'r') as f:
        lines = f.read().strip().split('\n')
    return float(lines[1].split('\t')[6])

def gc_at_percent(seq):
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    at = (seq.count('A') + seq.count('T')) / len(seq) * 100
    return gc, at

def nuc_freq(seq):
    freq = []
    for i in ['A', 'C', 'G', 'T']:
        freq.append(seq.count(i) / len(seq) * 100)
    return freq

def di_nuc_freq(seq):
    valid_di = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    total = len(seq) - 1
    # Count di-nucleotides
    counts = Counter(seq[i:i+2] for i in range(total) if seq[i:i+2] in valid_di)
    # Normalize frequencies
    freq = [counts[di] / total * 100 for di in valid_di]
    return freq

def count_motifs(seq, motifs):
    motifs = pd.read_csv(motifs) 
    counts = [1 if m in seq else 0 for m in motifs['motif'].tolist()]
    return np.array(counts).astype('float32').reshape(1, -1) # Reshape to 2D array for consistency

def run_encoder(counts, model):
    import tensorflow as tf
    encoder = tf.keras.models.load_model(model)
    return encoder.predict(counts)

def run_model(features, scaler, model_file):
    from sklearn.preprocessing import StandardScaler
    import xgboost as xgb
    scaler = joblib.load(scaler)
    features = scaler.transform(features)
    model = xgb.XGBClassifier()
    model.load_model(model_file)
    return model.predict_proba(features)[:, 1]

################## Main function ##################
def main(seq = None, type = None):
    seq = seq.strip()
    seq = ''.join(seq.split())

    if not is_valid_sequence(seq):
        print("Invalid sequence.")
        return

    seq = seq.upper().replace('U', 'T')        # Convert RNA to DNA

    files = []
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "seq.fa")
        save_fasta(seq, fasta_path)

        gc, at = gc_at_percent(seq)
        # RNAfold and bpRNA
        if type == 'mrna':
            structure, mfe = run_rnafold(fasta_path, tmpdir)
        elif type == 'circ':
            structure, mfe = run_rnafold(fasta_path, tmpdir, circ=True)

        sl, ml = run_bprna(tmpdir)

        bp, bp90, mld = structure_to_features(structure)

        # CPC2
        if type == 'circ':
            cpc = run_cpc2(fasta_path, tmpdir)

        # Combine features (example)
        if type == 'mrna':
            features = np.array([len(seq), gc, at, mfe*-1, mld, sl, ml, bp, bp90, 
                                *list(nuc_freq(seq)), *list(di_nuc_freq(seq))])
        elif type == 'circ':
            features = np.array([len(seq), gc, at, mfe*-1, mld, sl, ml, bp, bp90, cpc, 
                                *list(nuc_freq(seq)), *list(di_nuc_freq(seq))])
        
        if type == 'mrna':
            # counts = count_motifs(seq, os.path.join(os.path.dirname(__file__), 'mrna_motifs.csv'))
            counts = count_motifs(seq, os.path.join(os.path.dirname(__file__), 'model/attract_motifs.csv'))
            encoded = run_encoder(counts, os.path.join(os.path.dirname(__file__), 'model/mrna_encoder.keras'))
        elif type == 'circ':
            counts = count_motifs(seq, os.path.join(os.path.dirname(__file__), 'model/attract_motifs.csv'))
            encoded = run_encoder(counts, os.path.join(os.path.dirname(__file__), 'model/circ_encoder.keras'))

        model_input = np.concatenate((features, encoded[0])).reshape(1, -1)
        
        # run final model prediction
        if type == 'mrna':
            result = run_model(model_input, os.path.join(os.path.dirname(__file__), 'model/mrna_scaler.pkl'), os.path.join(os.path.dirname(__file__), 'model/xgb_mrna_model.json'))
        elif type == 'circ':
            result = run_model(model_input, os.path.join(os.path.dirname(__file__), 'model/circ_scaler.pkl'), os.path.join(os.path.dirname(__file__), 'model/xgb_circ_model.json'))


    if type == 'mrna':
        columns = ['Length', 'GC%', 'AT%', 'MFE', 'MLD', 'Stems', 'Multiloops', 'Base Pairs', 'BP90',
                  'A%', 'C%', 'G%', 'T%',
                  'ApA%', 'ApC%', 'ApG%', 'ApT%',
                  'CpA%', 'CpC%', 'CpG%', 'CpT%',
                  'GpA%', 'GpC%', 'GpG%', 'GpT%',
                  'TpA%', 'TpC%', 'TpG%', 'TpT%']
    elif type == 'circ':
        columns = ['Length', 'GC%', 'AT%', 'MFE', 'MLD', 'Stems', 'Multiloops', 'Base Pairs', 'BP90', 'Coding Probability',
                  'A%', 'C%', 'G%', 'T%',
                  'ApA%', 'ApC%', 'ApG%', 'ApT%',
                  'CpA%', 'CpC%', 'CpG%', 'CpT%',
                  'GpA%', 'GpC%', 'GpG%', 'GpT%',
                  'TpA%', 'TpC%', 'TpG%', 'TpT%']
        
    final = pd.DataFrame(features.reshape(1, -1), columns=columns)
    final['result'] = result
    return final

################## Motifs function ##################
def table_motifs(seq):
    seq = seq.strip()
    seq = ''.join(seq.split())

    if not is_valid_sequence(seq):
        print("Invalid sequence.")
        return

    seq = seq.upper().replace('U', 'T')        # Convert RNA to DNA

    motifs_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'model/attract_motifs.csv'))

    results = []
    for _, row in motifs_df.iterrows():
        motif = row['motif']
        matches = re.findall(motif, seq)
        count = len(matches)
        if count > 0:
            results.append({'RBP': row['gene'], 'motif': motif, 'count': count})

    results = pd.DataFrame(results)
    grouped = results.groupby('RBP').agg(
    sum_motif=('motif', 'count'),
    sum_count=('count', 'sum'),
    Motifs=('motif', lambda x: ','.join(x))
    )

    grouped['sum_all'] = grouped['sum_motif'] + grouped['sum_count']

    # Reset index if you want it as a column
    grouped = grouped.reset_index().sort_values(by='sum_all', ascending=False)

    return grouped

def plot_motifs(seq):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import linregress

    table = table_motifs(seq)

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=table, x='sum_motif', y='sum_count')

    # Fit light red line
    slope, intercept, _, _, _ = linregress(table['sum_motif'], table['sum_count'])
    x_vals = np.array(table['sum_motif'])
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, color='#ffcccc', linewidth=2)

    # Label top 5% genes
    threshold = table['sum_count'].quantile(0.95)
    top5 = table[table['sum_count'] >= threshold]
    for _, row in top5.iterrows():
        # plt.text(row['sum_motif'] + 0.3, row['sum_count'] + 0.3, f"{row['RBP']}\n({row['sum_motif']}, {row['sum_count']})", fontsize=9)
        plt.text(row['sum_motif'] + 0.3, row['sum_count'] + 0.3, row['RBP'], fontsize=9)

    plt.xlabel('Total number of motif occurrences')
    plt.ylabel('Number of distinct motifs')

    plt.margins(x=0.15, y=0.15)
    fig = plt.gcf()

    return fig


if __name__ == "__main__":
    main()
