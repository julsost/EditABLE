import math
from typing import List, Tuple
import argparse
from typing import Dict, Tuple
import pandas as pd

# Parameters for on target scoring
params: List[Tuple[int, str, float]] = [
    (1, 'G', -0.2753771), (2, 'A', -0.3238875), (2, 'C', 0.17212887), (3, 'C', -0.1006662),
    (4, 'C', -0.2018029), (4, 'G', 0.24595663), (5, 'A', 0.03644004), (5, 'C', 0.09837684),
    (6, 'C', -0.7411813), (6, 'G', -0.3932644), (11, 'A', -0.466099), (14, 'A', 0.08537695),
    (14, 'C', -0.013814), (15, 'A', 0.27262051), (15, 'C', -0.1190226), (15, 'T', -0.2859442),
    (16, 'A', 0.09745459), (16, 'G', -0.1755462), (17, 'C', -0.3457955), (17, 'G', -0.6780964),
    (18, 'A', 0.22508903), (18, 'C', -0.5077941), (19, 'G', -0.4173736), (19, 'T', -0.054307),
    (20, 'G', 0.37989937), (20, 'T', -0.0907126), (21, 'C', 0.05782332), (21, 'T', -0.5305673),
    (22, 'T', -0.8770074), (23, 'C', -0.8762358), (23, 'G', 0.27891626), (23, 'T', -0.4031022),
    (24, 'A', -0.0773007), (24, 'C', 0.28793562), (24, 'T', -0.2216372), (27, 'G', -0.6890167),
    (27, 'T', 0.11787758), (28, 'C', -0.1604453), (29, 'G', 0.38634258), (1, 'GT', -0.6257787),
    (4, 'GC', 0.30004332), (5, 'AA', -0.8348362), (5, 'TA', 0.76062777), (6, 'GG', -0.4908167),
    (11, 'GG', -1.5169074), (11, 'TA', 0.7092612), (11, 'TC', 0.49629861), (11, 'TT', -0.5868739),
    (12, 'GG', -0.3345637), (13, 'GA', 0.76384993), (13, 'GC', -0.5370252), (16, 'TG', -0.7981461),
    (18, 'GG', -0.6668087), (18, 'TC', 0.35318325), (19, 'CC', 0.74807209), (19, 'TG', -0.3672668),
    (20, 'AC', 0.56820913), (20, 'CG', 0.32907207), (20, 'GA', -0.8364568), (20, 'GG', -0.7822076),
    (21, 'TC', -1.029693), (22, 'CG', 0.85619782), (22, 'CT', -0.4632077), (23, 'AA', -0.5794924),
    (23, 'AG', 0.64907554), (24, 'AG', -0.0773007), (24, 'CG', 0.28793562), (24, 'TG', -0.2216372),
    (26, 'GT', 0.11787758), (28, 'GG', -0.69774)
]

#  mismatch scores and PAM scores dictionaries for off target scoring 
mismatch_scores: Dict[str, float] = {
    'rU:dT,12': 0.8, 'rU:dT,13': 0.692307692, 'rU:dC,5': 0.64,
    'rG:dA,14': 0.266666667, 'rG:dG,19': 0.448275862, 'rG:dG,18': 0.476190476,
    'rG:dG,15': 0.272727273, 'rG:dG,14': 0.428571429, 'rG:dG,17': 0.235294118,
    'rG:dG,16': 0.0, 'rC:dC,20': 0.058823529, 'rG:dT,20': 0.9375,
    'rG:dG,13': 0.421052632, 'rG:dG,12': 0.529411765, 'rU:dC,6': 0.571428571,
    'rU:dG,14': 0.285714286, 'rU:dT,18': 0.666666667, 'rA:dG,13': 0.210526316,
    'rA:dG,12': 0.263157895, 'rA:dG,11': 0.4, 'rA:dG,10': 0.333333333,
    'rA:dA,19': 0.538461538, 'rA:dA,18': 0.5, 'rA:dG,15': 0.272727273,
    'rA:dG,14': 0.214285714, 'rA:dA,15': 0.2, 'rA:dA,14': 0.533333333,
    'rA:dA,17': 0.133333333, 'rA:dA,16': 0.0, 'rA:dA,11': 0.307692308,
    'rA:dA,10': 0.882352941, 'rA:dA,13': 0.3, 'rA:dA,12': 0.333333333,
    'rG:dA,13': 0.3, 'rG:dA,12': 0.384615385, 'rG:dA,11': 0.384615385,
    'rG:dA,10': 0.8125, 'rG:dA,17': 0.25, 'rG:dA,16': 0.0, 'rG:dA,15': 0.142857143,
    'rG:dA,6': 0.666666667, 'rG:dG,20': 0.428571429, 'rG:dA,19': 0.666666667,
    'rG:dA,18': 0.666666667, 'rU:dC,4': 0.625, 'rG:dT,12': 0.933333333,
    'rG:dT,13': 0.923076923, 'rU:dG,11': 0.666666667, 'rC:dA,3': 0.6875,
    'rC:dA,2': 0.909090909, 'rC:dA,1': 1.0, 'rC:dA,7': 0.8125,
    'rC:dA,6': 0.928571429, 'rC:dA,5': 0.636363636, 'rC:dA,4': 0.8,
    'rC:dA,9': 0.875, 'rC:dA,8': 0.875, 'rU:dT,6': 0.866666667,
    'rA:dG,20': 0.227272727, 'rG:dT,18': 0.692307692, 'rU:dG,10': 0.533333333,
    'rG:dT,19': 0.714285714, 'rG:dA,20': 0.7, 'rC:dT,20': 0.5,
    'rU:dC,2': 0.84, 'rG:dG,10': 0.4, 'rC:dA,17': 0.466666667,
    'rC:dA,16': 0.307692308, 'rC:dA,15': 0.066666667, 'rC:dA,14': 0.733333333,
    'rC:dA,13': 0.7, 'rC:dA,12': 0.538461538, 'rC:dA,11': 0.307692308,
    'rC:dA,10': 0.941176471, 'rG:dG,11': 0.428571429, 'rU:dC,20': 0.176470588,
    'rG:dG,3': 0.384615385, 'rC:dA,19': 0.461538462, 'rC:dA,18': 0.642857143,
    'rU:dG,17': 0.705882353, 'rU:dG,16': 0.666666667, 'rU:dG,15': 0.272727273,
    'rG:dG,2': 0.692307692, 'rU:dG,13': 0.789473684, 'rU:dG,12': 0.947368421,
    'rG:dA,9': 0.533333333, 'rG:dA,8': 0.625, 'rG:dA,7': 0.571428571,
    'rG:dG,5': 0.785714286, 'rG:dA,5': 0.3, 'rG:dA,4': 0.363636364,
    'rG:dA,3': 0.5, 'rG:dA,2': 0.636363636, 'rG:dA,1': 1.0,
    'rG:dG,4': 0.529411765, 'rG:dG,1': 0.714285714, 'rA:dC,9': 0.666666667,
    'rG:dG,7': 0.6875, 'rG:dT,5': 0.866666667, 'rU:dT,20': 0.5625,
    'rC:dC,15': 0.05, 'rC:dC,14': 0.0, 'rC:dC,17': 0.058823529,
    'rC:dC,16': 0.153846154, 'rC:dC,11': 0.25, 'rC:dC,10': 0.388888889,
    'rC:dC,13': 0.136363636, 'rC:dC,12': 0.444444444, 'rC:dA,20': 0.3,
    'rC:dC,19': 0.125, 'rC:dC,18': 0.133333333, 'rA:dA,1': 1.0,
    'rA:dA,3': 0.705882353, 'rA:dA,2': 0.727272727, 'rA:dA,5': 0.363636364,
    'rA:dA,4': 0.636363636, 'rA:dA,7': 0.4375, 'rA:dA,6': 0.714285714,
    'rA:dA,9': 0.6, 'rA:dA,8': 0.428571429, 'rU:dG,20': 0.090909091,
    'rC:dC,9': 0.619047619, 'rC:dC,8': 0.642857143, 'rU:dT,10': 0.857142857,
    'rU:dT,11': 0.75, 'rU:dT,16': 0.909090909, 'rU:dT,17': 0.533333333,
    'rU:dT,14': 0.619047619, 'rU:dT,15': 0.578947368, 'rC:dC,1': 0.913043478,
    'rU:dT,3': 0.714285714, 'rC:dC,3': 0.5, 'rC:dC,2': 0.695652174,
    'rC:dC,5': 0.6, 'rC:dC,4': 0.5, 'rC:dC,7': 0.470588235,
    'rC:dC,6': 0.5, 'rU:dT,4': 0.476190476, 'rU:dT,8': 0.8,
    'rU:dT,9': 0.928571429, 'rA:dC,19': 0.375, 'rA:dC,18': 0.4,
    'rA:dC,17': 0.176470588, 'rA:dC,16': 0.192307692, 'rA:dC,15': 0.65,
    'rA:dC,14': 0.466666667, 'rA:dC,13': 0.652173913, 'rA:dC,12': 0.722222222,
    'rA:dC,11': 0.65, 'rA:dC,10': 0.555555556, 'rU:dC,7': 0.588235294,
    'rC:dT,8': 0.65, 'rC:dT,9': 0.857142857, 'rC:dT,6': 0.928571429,
    'rC:dT,7': 0.75, 'rC:dT,4': 0.842105263, 'rC:dT,5': 0.571428571,
    'rC:dT,2': 0.727272727, 'rC:dT,3': 0.866666667, 'rC:dT,1': 1.0,
    'rA:dC,8': 0.733333333, 'rU:dT,1': 1.0, 'rU:dC,3': 0.5,
    'rU:dC,1': 0.956521739, 'rU:dT,2': 0.846153846, 'rU:dG,19': 0.275862069,
    'rG:dT,14': 0.75, 'rG:dT,15': 0.941176471, 'rG:dT,16': 1.0,
    'rG:dT,17': 0.933333333, 'rG:dT,10': 0.933333333, 'rG:dT,11': 1.0,
    'rA:dG,9': 0.571428571, 'rA:dG,8': 0.428571429, 'rA:dG,7': 0.4375,
    'rA:dG,6': 0.454545455, 'rA:dG,5': 0.5, 'rA:dG,4': 0.352941176,
    'rA:dG,3': 0.428571429, 'rA:dG,2': 0.785714286, 'rA:dG,1': 0.857142857,
    'rU:dT,5': 0.5, 'rG:dT,2': 0.846153846, 'rA:dC,3': 0.611111111,
    'rA:dC,20': 0.764705882, 'rG:dT,1': 0.9, 'rG:dT,6': 1.0,
    'rG:dT,7': 1.0, 'rG:dT,4': 0.9, 'rC:dT,19': 0.428571429,
    'rG:dG,9': 0.538461538, 'rG:dG,8': 0.615384615, 'rG:dT,8': 1.0,
    'rG:dT,9': 0.642857143, 'rU:dG,18': 0.428571429, 'rU:dT,7': 0.875,
    'rG:dG,6': 0.681818182, 'rA:dA,20': 0.6, 'rU:dC,9': 0.619047619,
    'rA:dG,17': 0.176470588, 'rU:dC,8': 0.733333333, 'rA:dG,16': 0.0,
    'rA:dG,19': 0.206896552, 'rG:dT,3': 0.75, 'rU:dG,3': 0.428571429,
    'rU:dG,2': 0.857142857, 'rU:dG,1': 0.857142857, 'rA:dG,18': 0.19047619,
    'rU:dG,7': 0.6875, 'rU:dG,6': 0.909090909, 'rU:dG,5': 1.0,
    'rU:dG,4': 0.647058824, 'rU:dG,9': 0.923076923, 'rU:dG,8': 1.0,
    'rU:dC,19': 0.25, 'rU:dC,18': 0.333333333, 'rU:dC,13': 0.260869565,
    'rU:dC,12': 0.5, 'rU:dC,11': 0.4, 'rU:dC,10': 0.5,
    'rU:dC,17': 0.117647059, 'rU:dC,16': 0.346153846, 'rU:dC,15': 0.05,
    'rU:dC,14': 0.0, 'rC:dT,10': 0.866666667, 'rC:dT,11': 0.75,
    'rC:dT,12': 0.714285714, 'rC:dT,13': 0.384615385, 'rC:dT,14': 0.35,
    'rC:dT,15': 0.222222222, 'rC:dT,16': 1.0, 'rC:dT,17': 0.466666667,
    'rC:dT,18': 0.538461538, 'rA:dC,2': 0.8, 'rA:dC,1': 1.0,
    'rA:dC,7': 0.705882353, 'rA:dC,6': 0.714285714, 'rA:dC,5': 0.72,
    'rA:dC,4': 0.625, 'rU:dT,19': 0.285714286,
}

pam_scores: Dict[str, float] = {
    'AA': 0.0, 'AC': 0.0, 'GT': 0.016129032, 'AG': 0.259259259,
    'CC': 0.0, 'CA': 0.0, 'CG': 0.107142857, 'TT': 0.0,
    'GG': 1.0, 'GC': 0.022222222, 'AT': 0.0, 'GA': 0.069444444,
    'TG': 0.038961039, 'TA': 0.0, 'TC': 0.0, 'CT': 0.0,
}

intercept = 0.59763615
gcHigh = -0.1665878
gcLow = -0.2026259

def calculate_on_target_scores(seqs: List[str]) -> pd.DataFrame:
    scores = []
    for seq in seqs:
        score = intercept
        guide_seq = seq[4:24]
        gc_count = guide_seq.count("G") + guide_seq.count("C")
        if gc_count <= 10:
            gc_weight = gcLow
        else:
            gc_weight = gcHigh
        score += abs(10 - gc_count) * gc_weight

        for pos, model_seq, weight in params:
            sub_seq = seq[pos:pos + len(model_seq)]
            if sub_seq == model_seq:
                score += weight
        Score = 1.0 / (1.0 + math.exp(-score))
        scores.append((seq, Score))
    
    # Create a DataFrame for the results
    on_target_scores_df = pd.DataFrame(scores, columns=['sequence', 'score'])
    on_target_scores_df['score'] = on_target_scores_df['score'].round(3)*100
    return on_target_scores_df





#Everything below is in relation to off target scoring
def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Calculates CFD score')
    parser.add_argument('--wt', type=str, required=True, help='WT 23mer sgRNA sequence')
    parser.add_argument('--off', type=str, required=True, help='Off-target 23mer sgRNA sequence')
    return parser

# Reverse complements a given string
def revcom(s: str) -> str:
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = [basecomp[base] for base in s[::-1]]
    return ''.join(letters)


def calculate_off_target_scores(wt_sequences: List[str], sg_sequences: List[str], pam_sequences: List[str]) -> pd.DataFrame:
    scores = []
    for wt, sg, pam in zip(wt_sequences, sg_sequences, pam_sequences):
        score = 1.0
        sg = sg.replace('T', 'U')
        wt = wt.replace('T', 'U')
        s_list = list(sg)
        wt_list = list(wt)

        for i, sl in enumerate(s_list):
            if wt_list[i] == sl:
                score *= 1
            else:
                key = f'r{wt_list[i]}:d{revcom(sl)},{i+1}'
                score *= mismatch_scores.get(key, 1.0)

        score *= pam_scores.get(pam, 1.0)
        scores.append((wt, sg, pam, score))

    # Create a DataFrame for the results
    cfd_scores_df = pd.DataFrame(scores, columns=['spacer', 'protospacer', 'pam', 'score'])
    
    return cfd_scores_df


