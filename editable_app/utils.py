import collections
import pandas as pd
import regex as re
from Bio.Seq import MutableSeq
from Bio.SeqUtils import MeltingTemp as mt
from math import log
from itertools import combinations 
import numpy as np
from calculate_scores import calculate_on_target_scores, calculate_off_target_scores 


bases = {"A", "C", "G", "T"}

edit_start = 4
edit_end = 9
gRNA_size = 20
pam_length = None
PAMs = None
reverse_PAMs = None

#helper function to find be guide rnas, pams, and 30 nt sequence
def find_BE_guide_rnas(direction, seq, genomic_location):
    print(f"find_BE_guide_rnas called with direction: {direction}, genomic_location: {genomic_location}")
    if direction == "forward":
        start_pams = genomic_location + (gRNA_size - edit_end) + 1
        end_pams = genomic_location + (gRNA_size - edit_start) + pam_length + 1
        possible_pams = seq[start_pams : end_pams]
        print(f"possible_pams (forward): {possible_pams}")
        if len(re.findall(PAMs, str(possible_pams))) == 0:
            return "NO PAM"
        else:
            guideRNAs = []
            pams = []
            sequences = []
            for match in re.finditer(PAMs, str(possible_pams), overlapped=True):
                gRNA = seq[start_pams + match.start() - gRNA_size: start_pams + match.start()]
                pam = possible_pams[match.start(): match.start() + pam_length]
                guideRNAs.append(str(gRNA))
                pams.append(str(pam))
                start_idx = start_pams + match.start() - gRNA_size - 4
                end_idx = start_pams + match.start() + pam_length + 3
                sequence = seq[start_idx:end_idx]
                sequences.append(str(sequence))
                print(f"Found guideRNA (forward): {gRNA}, PAM: {pam}, sequence: {sequence}")
            return guideRNAs, pams, sequences

    elif direction == "reverse":
        start_pams = genomic_location - gRNA_size + edit_start - pam_length
        end_pams = genomic_location - gRNA_size + edit_end
        possible_pams = seq[start_pams : end_pams]
        print(f"possible_pams (reverse): {possible_pams}")
        if len(re.findall(reverse_PAMs, str(possible_pams))) == 0:
            return "NO PAM"
        else:
            guideRNAs = []
            pams = []
            sequences = []
            for match in re.finditer(reverse_PAMs, str(possible_pams), overlapped=True):
                gRNA_start = start_pams + match.start() + pam_length
                gRNA_end = gRNA_start + gRNA_size
                gRNA = seq[gRNA_start : gRNA_end]
                pam = possible_pams[match.start(): match.start() + pam_length]
                guideRNAs.append(str(gRNA.reverse_complement()))
                pams.append(str(pam.reverse_complement()))  # Correct the PAM to be the reverse complement
                start_idx = start_pams + match.start() - 3
                end_idx = start_pams + match.start() + gRNA_size + pam_length + 4
                sequence = seq[start_idx:end_idx].reverse_complement()
                sequences.append(str(sequence))
                print(f"Found guideRNA (reverse): {gRNA}, PAM: {pam}, sequence: {sequence}")
            return guideRNAs, pams, sequences
    else:
        return "ERROR"

    
 
 # Helper function to find guide RNAs for transversion (C>G and G>C) mutations on forward/reverse strands    
def find_CGB_guide_rnas(direction, seq, genomic_location):
    guideRNAs = []
    pams = []
    sequences = []

    if direction == "forward":
        pam = seq[genomic_location + 15: genomic_location + 15 + pam_length]
        if len(re.findall(PAMs, str(pam))) == 0:
            return "NO PAM"
        else:
            guideRNA = seq[genomic_location + 15 - gRNA_size: genomic_location + 15]
            guideRNAs.append(str(guideRNA))
            pams.append(str(pam))
            start_idx = genomic_location + 15 - gRNA_size - 4
            end_idx = genomic_location + 15 + pam_length + 3
            sequence = seq[start_idx:end_idx]
            sequences.append(str(sequence))

    elif direction == "reverse":
        pam = seq[genomic_location - 15 - 1: genomic_location - 15 - 1 + pam_length]
        if len(re.findall(reverse_PAMs, str(pam))) == 0:
            return "NO PAM"
        else:
            guideRNA = seq[genomic_location - 15 - 1 + pam_length: genomic_location - 15 - 1 + pam_length + gRNA_size]
            guideRNA.reverse_complement(inplace=True)
            guideRNAs.append(str(guideRNA))
            pams.append(str(pam.reverse_complement()))  # Correct the PAM to be the reverse complement
            start_idx = genomic_location - 15 - 1 + pam_length - 4
            end_idx = genomic_location - 15 - 1 + pam_length + gRNA_size + 3
            sequence = seq[start_idx:end_idx].reverse_complement()
            sequences.append(str(sequence))
            
        print(f"guide rnas: {guideRNAs}")
        print(f"pams: {pams}")

    else:
        return "ERROR"

    return guideRNAs, pams, sequences

# dictionary for guide rna cloning plasmid
pam_to_url = {
    "NNNRRT_NNGRRT_NNGRR": ( 'BPK2660', '70709',  "https://www.addgene.org/70709"),
    "NGG_NGN_NRN_NYN_NGA": ('pmCherry-U6-empty', '140580', "https://www.addgene.org/140580")
}

def get_cloning_url(pam):
    pam_patterns = {
        "NNNRRT": "NNNRRT_NNGRRT_NNGRR",
        "NNGRRT": "NNNRRT_NNGRRT_NNGRR",
        "NNGRR": "NNNRRT_NNGRRT_NNGRR",
        "NGG": "NGG_NGN_NRN_NYN_NGA",
        "NGN": "NGG_NGN_NRN_NYN_NGA",
        "NRN": "NGG_NGN_NRN_NYN_NGA",
        "NYN": "NGG_NGN_NRN_NYN_NGA",
        "NGA": "NGG_NGN_NRN_NYN_NGA"
    }
    
    cloning_group = pam_patterns.get(pam)
    url_info = pam_to_url.get(cloning_group)
    
    return url_info
 
#helper function to track position of base to edit within the guide rna and then return potential bystander edits   
def track_positions(guide, ref_sequence_input, substitution_position, orientation):
    ref_sequence_almost_rc = almost_reverse_complement(ref_sequence_input)
    base_to_edit = ref_sequence_input[substitution_position] if orientation != 'reverse' else ref_sequence_almost_rc[substitution_position]
    
    if orientation == 'reverse':
        guide = guide[::-1]
        all_guide_occurance_starts = [m.start() for m in re.finditer(guide, ref_sequence_almost_rc)]
    else:
        all_guide_occurance_starts = [m.start() for m in re.finditer(guide, ref_sequence_input)]

    true_starting_positions = []
    for start in all_guide_occurance_starts:
        end = start + len(guide) - 1
        if substitution_position >= start and substitution_position <= end:
            true_starting_positions.append(start)
    
    assert len(true_starting_positions) == 1, ("Error! Guide cannot be aligned properly to input original sequence", guide, ref_sequence_almost_rc, orientation, all_guide_occurance_starts, substitution_position)
    guide_start = true_starting_positions[0]

    if orientation == 'reverse':
        window_sequence = ref_sequence_almost_rc[guide_start + len(guide) - edit_end: guide_start + len(guide) - edit_start + 1][::-1]
        actual_base_position = guide_start + len(guide) - substitution_position 
    else:
        window_sequence = ref_sequence_input[guide_start + edit_start - 1: guide_start + edit_end]
        actual_base_position = substitution_position - guide_start + 1   # Adjust to 0-based index for the window

    all_positions = [i for i, base in enumerate(window_sequence, start=edit_start) if base == base_to_edit]
    
    if actual_base_position in all_positions:
        all_positions.remove(actual_base_position)

    if len(all_positions) > 0:
        return f"Yes (positions: {', '.join(map(str, all_positions))})"
    else:
        return "No"

#helper function for experimental validation of base editing guide rnas
def process_guide_rnas(guide_rnas):
    results = []
    for index, guide_rna in enumerate(guide_rnas):
        # Replace the first letter with G
        modified_guide_rna = 'G' + guide_rna[1:]
        
        # Generate the reverse complement
        reverse_guide = str(reverse_complement(modified_guide_rna))
        
        # Replace the last letter of the reverse complement with C
        reverse_complement_modified = reverse_guide[:-1] + 'C'
        
        # Format the output
        result = f"Guide {index + 1}\n5' - CACC {modified_guide_rna} - 3'\n5' - AAAC {reverse_complement_modified} - 3'"
        results.append(result)
    
    return results


def process_ng_rnas(oligo_top):
    results = []
    guide_rnas = [oligo[4:24] for oligo in oligo_top if oligo and oligo != 'n/a']  # Extract the guide RNA sequences
    for index, guide_rna in enumerate(guide_rnas):

        # Generate the reverse complement
        reverse_guide = str(reverse_complement(guide_rna))

        # Format the output
        result = f"5' - CACC {guide_rna} - 3'\n5' - AAAC {reverse_guide} - 3'"
        results.append(result)

    return results

def process_peg_rnas(oligo_top):
    results = []
    guide_rnas = [oligo[4:24] for oligo in oligo_top if oligo and oligo != 'n/a']  # Extract the guide RNA sequences
    for index, guide_rna in enumerate(guide_rnas):
        
        # Generate the reverse complement
        reverse_guide = str(reverse_complement(guide_rna))

        # Format the output
        result = f"5' - ACGggtctcg CACC {guide_rna} GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC {reverse_guide} CGCG ggagaccGTA - 3'"
        results.append(result)

    return results

#function to create array of protospacers with mutations
def generate_mutations_to_single_base(guide_rnas, max_mutations=3):
    bases = {'A', 'C', 'G', 'T'}
    all_mutated_guides = []

    def mutate_single_position(guide, position, new_base):
        new_guide = guide[:]
        new_guide[position] = new_base
        return ''.join(new_guide)

    for guide_rna in guide_rnas:
        guide_rna_list = list(guide_rna)
        for target_base in bases:
            mutated_guides = []
            for num_mutations in range(1, max_mutations + 1):
                for positions in combinations(range(len(guide_rna)), num_mutations):
                    new_guide = guide_rna_list[:]
                    for pos in positions:
                        new_guide[pos] = target_base
                    mutated_guides.append(''.join(new_guide))
            all_mutated_guides.extend(mutated_guides)

    return all_mutated_guides

def calculate_off_target_scores_for_guides(guide_rnas, pams):
    all_results = []

    for guide_rna, pam in zip(guide_rnas, pams):
        protospacers = generate_mutations_to_single_base([guide_rna])
        scores_df = calculate_off_target_scores([guide_rna] * len(protospacers), protospacers, [pam] * len(protospacers))
        scores_df['guide_rna'] = guide_rna
        all_results.append(scores_df)
    
    combined_scores = pd.concat(all_results)

    # this is to ensure only numeric columns are averaged
    numeric_columns = combined_scores.select_dtypes(include=[np.number]).columns.tolist()
    numeric_columns.append('guide_rna')

    average_scores = combined_scores[numeric_columns].groupby('guide_rna').mean().reset_index()
    average_scores['score'] = 100-(average_scores['score'].round(3)*100)

    return combined_scores, average_scores


#Define editor data, dictionary for type of editing enzymes
editor_data = {
    'ABE': {
        'NGN': ('ABE8e-NG enzyme', '138491', 'https://www.addgene.org/138491'),
        'NGG': ('ABE8e enzyme', '138489', 'https://www.addgene.org/138489'),
        'NGA': ('ABE8e-NG enzyme', '138491', 'https://www.addgene.org/138491'),
        'NNGRRT': ('ABE8e-SaCas enzyme', '138500', 'https://www.addgene.org/138500'),
        'NNNRRT': ('ABE8e-SaCas-KKH enzyme', '138502', 'https://www.addgene.org/138502'),
        'NRN': ('SpRY', '140003', 'https://www.addgene.org/140003/'),
        'NYN': ('SpRY', '140003', 'https://www.addgene.org/140003/')
    },
    'CBE': {
        'NGN': ('BE4max-NG', '138159', 'https://www.addgene.org/138159'),
        'NGG': ('BE4max', '112093', 'https://www.addgene.org/112093'),
        'NGA': ('BE4max-NG', '138159', 'https://www.addgene.org/138159'),
        'NNGRRT': ('BE3-SaCas', '85169', 'https://www.addgene.org/85169'),
        'NNNRRT': ('BE3-SaCas-KKH', '85170', 'https://www.addgene.org/85170'),
        'NRN' : ('SpRY', '139999', 'https://www.addgene.org/139999/'),
        'NYN' : ('SpRY', '139999', 'https://www.addgene.org/139999/')
    },
    'CGB': {
        'NGN' : ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/' ),
        'NGG': ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/' ),
        'NGA': ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/'),
        'NNGRRT': ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/'),
        'NNNRRT': ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/'),
        'NRN' : ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/'),
        'NYN' : ('UdgX-HF-nCas9','163559','https://www.addgene.org/163559/')
    },
    'PrimeEditor': {
    'default': ('PE2', '132775', 'https://www.addgene.org/132775')
    },
}    

#determine the mutation type for the base editors
def get_editor_info(ref_seq, edited_seq, pams):
    mutation_type = "PrimeEditor"  # Default mutation type is PrimeEditor
    for ref_base, edit_base in zip(ref_seq, edited_seq):
        if ref_base != edit_base:
            # Create a mutation identifier from reference and edited bases
            edit_type = f"{ref_base}>{edit_base}"
            
            # Determine the type of base editing required based on the mutation
            if edit_type in ["A>G", "T>C"]:
                mutation_type = "ABE"
            elif edit_type in ["G>A", "C>T"]:
                mutation_type = "CBE"
            elif edit_type in ["C>G", "G>C"]:
                mutation_type = "CGB"   
    return editor_data[mutation_type].get(pams, editor_data['PrimeEditor']['default'])


# Adds the base editable guide RNAs to the data table for export
def get_guide_RNAs(mutant_seq, edit_type, genomic_location):
    try:
        if edit_type == "A>G":
            result = find_BE_guide_rnas("forward", mutant_seq, genomic_location)
        elif edit_type == "T>C":
            result = find_BE_guide_rnas("reverse", mutant_seq, genomic_location)
        elif edit_type == "G>A":
            result = find_BE_guide_rnas("reverse", mutant_seq, genomic_location)
        elif edit_type == "C>T":
            result = find_BE_guide_rnas("forward", mutant_seq, genomic_location)
        elif edit_type == "G>C":
            result = find_CGB_guide_rnas("forward", mutant_seq, genomic_location)
        elif edit_type == "C>G":
            result = find_CGB_guide_rnas("reverse", mutant_seq, genomic_location)
        else:
            return (("NOT BASE EDITABLE", None, None), None)

        if result == "NO PAM":
            return (("NO PAM", None, None), None)
        return result, "forward" if edit_type in ["A>G", "C>T"] else "reverse"

    except Exception as e:
        return (("ERROR", None, None), None)


#helper function to set pam sequences: Initialize PAM sequences based on the provided PAM
def set_pam_sequences(PAM):
    global PAMs, reverse_PAMs, pam_length
    pam_length = len(PAM)
    if PAM == 'NGN':
        PAMs = re.compile("[A|T|G|C]G[A|T|G|C]")
        reverse_PAMs = re.compile("[A|T|G|C]C[A|T|G|C]")
    elif PAM == 'NGG':
        PAMs = re.compile("[A|T|G|C]GG")
        reverse_PAMs = re.compile("CC[A|T|G|C]")
    elif PAM == 'NGA':
        PAMs = re.compile("[A|T|G|C]GA")
        reverse_PAMs = re.compile("TC[A|T|G|C]")
    elif PAM == 'NNGRRT':
        PAMs = re.compile("[A|T|G|C][A|T|G|C]G[A|G][A|G]T")
        reverse_PAMs = re.compile("A[T|C][T|C]C[A|T|G|C][A|T|G|C]")
    elif PAM == 'NNNRRT':
        PAMs = re.compile("[A|T|G|C][A|T|G|C][A|T|G|C][A|G][A|G]T")
        reverse_PAMs = re.compile("A[T|C][T|C][A|T|G|C][A|T|G|C][A|T|G|C]")
    elif PAM == 'NRN':
        PAMs = re.compile("[A|T|G|C][G|A][A|T|G|C]")
        reverse_PAMs = re.compile("[A|T|G|C][C|T][A|T|G|C]")   
    elif PAM == 'NYN':
        PAMs = re.compile("[A|T|G|C][C|T][A|T|G|C]")
        reverse_PAMs = re.compile("[A|T|G|C][G|A][A|T|G|C]")   
    else:
        raise ValueError("Invalid PAM sequence")
 
 #helper function to create mutuable seq: Convert input sequences to mutable sequences  
def create_mutable_sequences(ref_sequence_original, edited_sequence_original):
    ref_sequence = MutableSeq(ref_sequence_original)
    edited_sequence = MutableSeq(edited_sequence_original)
    return ref_sequence, edited_sequence

#helper function to identify substitution position: Find the position of the substitution in the sequences
def identify_substitution_position(ref_sequence, edited_sequence):
    substitution_position = None
    for i in range(len(ref_sequence)):
        ref_base = ref_sequence[i]
        edited_base = edited_sequence[i]
        if ref_base != edited_base:
            substitution_position = i
            break
    return substitution_position

#helper function to determine the appropriate guide RNAs based on the mutation type.
def get_guide_rnas_and_orientation(ref_sequence, edit, substitution_position):
    try:
        result = get_guide_RNAs(ref_sequence, edit, substitution_position)

        if result is None or result[0] is None:
            raise ValueError("Invalid result from get_guide_RNAs")
        
        (guide_rnas, pams, sequences), orientation = result
        
        return guide_rnas, orientation, pams, sequences
    
    except ValueError as ve:
        print(f"ValueError: {ve}")
        print(f"result: {result}")
    except Exception as e:
        print(f"Exception: {e}")
        print(f"result: {result}")


#helper function to generate prime design outputs if necessary
def run_prime_design(ref_sequence, edited_sequence, substitution_position):
    primedesign_input = str(ref_sequence[:substitution_position] + f"({ref_sequence[substitution_position]}/{edited_sequence[substitution_position]})" + ref_sequence[substitution_position + 1:])
    primedesign_output = run_primedesign(str(primedesign_input))
    return primedesign_output

#helper function to handle deletions: process cases where there are insertions.
def handle_insertions(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original):
    insertion_start = ref_sequence.find('-')
    insertion_end = ref_sequence.rfind('-')
    if insertion_end - insertion_start + 1 > 44:
        add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "Use Twin Prime Editing/Integrase/HDR")
    else:
        primedesign_input = str(ref_sequence[:insertion_start] + f"(+{edited_sequence[insertion_start : insertion_end + 1]})" + ref_sequence[insertion_end + 1:])
        primedesign_output = run_primedesign(primedesign_input)
        if primedesign_output != "No PrimeDesign Recommended Guides":
            update_df_dict_with_primedesign_output(df_dict, ref_sequence_original, edited_sequence_original, primedesign_output, "Prime Editing")
        else:
            add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "Use Twin Prime Editing/Integrase/HDR")
    return df_dict


#helper function to handle deletions: process cases where there are deletions.
def handle_deletions(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original):
    deletion_start = edited_sequence.find('-')
    deletion_end = edited_sequence.rfind('-')
    if deletion_start - deletion_end + 1 > 80:
        add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "Use Twin Prime Editing/Integrase/HDR")
    else:
        primedesign_input = str(ref_sequence[:deletion_start] + f"(-{ref_sequence[deletion_start : deletion_end + 1]})" + ref_sequence[deletion_end + 1:])
        primedesign_output = run_primedesign(primedesign_input)
        if primedesign_output != "No PrimeDesign Recommended Guides":
            update_df_dict_with_primedesign_output(df_dict, ref_sequence_original, edited_sequence_original, primedesign_output, "Prime Editing")
        else:
            add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "Use Twin Prime Editing/Integrase/HDR")
    return df_dict

#helper function to add insertions and deletions
def add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, editing_technology):
    df_dict['Original Sequence'].append(ref_sequence_original)
    df_dict['Desired Sequence'].append(edited_sequence_original)
    df_dict['Editing Technology'].append(editing_technology)
    df_dict['Base Editing Guide'].append(None)
    df_dict['Off Target Score'].append(None)
    df_dict['On Target Score'].append(None)
    df_dict['Bystander Edits?'].append(None)
    df_dict['PrimeDesign pegRNA Annotation'].append(None)
    df_dict['PrimeDesign pegRNA PBS'].append(None)
    df_dict['PrimeDesign pegRNA RTT'].append(None)
    df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
    #df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
    df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
    #df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
    df_dict['PrimeDesign ngRNA Annotation'].append(None)
    df_dict['PrimeDesign ngRNA Distance'].append(None)
    df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
    #df_dict['PrimeDesign ngRNA Bottom Top'].append(None)

#helper function to update the dataframe dictionary with prime design output
def update_df_dict_with_primedesign_output(df_dict, ref_sequence_original, edited_sequence_original, primedesign_output, editing_technology):
    peg_spacer_top_recommended, peg_spacer_bottom_recommended, peg_ext_top_recommended, peg_ext_bottom_recommended, peg_annotation_recommended, peg_pbs_recommended, peg_rtt_recommended, ng_spacer_top_recommended, ng_spacer_bottom_recommended, ng_annotation_recommended, ng_distance_recommended = primedesign_output
    df_dict['Original Sequence'].append(ref_sequence_original)
    df_dict['Desired Sequence'].append(edited_sequence_original)
    df_dict['Editing Technology'].append(editing_technology)
    df_dict['PrimeDesign pegRNA Annotation'].append(peg_annotation_recommended)
    df_dict['PrimeDesign pegRNA PBS'].append(peg_pbs_recommended)
    df_dict['PrimeDesign pegRNA RTT'].append(peg_rtt_recommended)
    df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(peg_spacer_top_recommended)
    #df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(peg_spacer_bottom_recommended)
    df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(peg_ext_top_recommended)
    #df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(peg_ext_bottom_recommended)
    df_dict['PrimeDesign ngRNA Annotation'].append(ng_annotation_recommended)
    #df_dict['PrimeDesign ngRNA Distance'].append(ng_distance_recommended)
    df_dict['PrimeDesign ngRNA Oligo Top'].append(ng_spacer_top_recommended)
    #df_dict['PrimeDesign ngRNA Bottom Top'].append(ng_spacer_bottom_recommended)

def render_dataframe(df_dict):
    df_dict_render = collections.defaultdict(list)
    max_length = max(len(v) for v in df_dict.values())

    for key in df_dict:
        # Ensure each list has the same length by padding with None
        while len(df_dict[key]) < max_length:
            df_dict[key].append(None)
        df_dict_render[key] = df_dict[key]

    df_render = pd.DataFrame.from_dict(df_dict_render).dropna(how='all', axis=1)
    df_full = pd.DataFrame.from_dict(df_dict).dropna(how='all', axis=1)

    # Drop the 'Base Editing Guide Orientation' column if it exists
    if 'Base Editing Guide Orientation' in df_render.columns:
        df_render = df_render.drop(columns=['Base Editing Guide Orientation'])
        
    if 'Base Editing Guide Pam' in df_render.columns:
        df_render = df_render.drop(columns=['Base Editing Guide Pam'])
 
    # Conditionally filter 'Base Editing Guide' if it exists
    if 'Base Editing Guide' in df_render.columns:
        df_render = df_render[df_render['Base Editing Guide'].notna()]
        df_render = df_render[df_render['Base Editing Guide'] != '']
    
    return df_render, df_full

def get_guides(ref_sequence_original, edited_sequence_original, PAM):
    def process_guide_rnas_for_pam(guide_rnas, sequences, ref_sequence_original, edited_sequence_original, df_dict, ref_sequence, edited_sequence, substitution_position, orientation, pam):
        if guide_rnas == "NO PAM":
            return False  # Indicates that the next PAM should be tried
        elif guide_rnas == "NOT BASE EDITABLE":
            primedesign_output = run_prime_design(ref_sequence, edited_sequence, substitution_position)
            if primedesign_output != "No PrimeDesign Recommended Guides":
                update_df_dict_with_primedesign_output(df_dict, ref_sequence_original, edited_sequence_original, primedesign_output, "Prime Editing")
            else:
                add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "No Base or Prime Editing Guides Found")
            return True  # Indicates that further processing is complete
        else:
            protospacers = generate_mutations_to_single_base(guide_rnas)
            try:
                combined_scores, average_scores = calculate_off_target_scores_for_guides(guide_rnas, [pam] * len(guide_rnas))
            except Exception as e:
                print(f"Error in calculate_off_target_scores: {e}")
                print(f"Guide RNAs: {guide_rnas}")
                combined_scores, average_scores = pd.DataFrame(), pd.DataFrame()

            on_target_scores_df = calculate_on_target_scores(sequences)

            if not isinstance(guide_rnas, list):
                guide_rnas = [guide_rnas]
            for gRNA, cfd_score, on_score in zip(guide_rnas, average_scores['score'], on_target_scores_df['score']):
                position_info = track_positions(gRNA, ref_sequence_original, substitution_position, orientation)

                df_dict['Original Sequence'].append(ref_sequence_original)
                df_dict['Desired Sequence'].append(edited_sequence_original)
                df_dict['Editing Technology'].append("Base Editing")
                df_dict['Base Editing Guide'].append(str(gRNA))
                df_dict['Base Editing Guide Orientation'].append(orientation)
                df_dict['Base Editing Guide Pam'].append(pam)
                df_dict['Off Target Score (Click To Toggle)'].append(cfd_score)
                df_dict['On Target Score (Click To Toggle)'].append(on_score)
                df_dict['Bystander Edits?'].append(position_info)
                df_dict['PrimeDesign pegRNA Annotation'].append(None)
                df_dict['PrimeDesign pegRNA PBS'].append(None)
                df_dict['PrimeDesign pegRNA RTT'].append(None)
                df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
                #df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
                df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
                #df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
                df_dict['PrimeDesign ngRNA Annotation'].append(None)
                #df_dict['PrimeDesign ngRNA Distance'].append(None)
                df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
                #df_dict['PrimeDesign ngRNA Bottom Top'].append(None)
            return True

    def handle_base_editing(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original, substitution_position, pam_sequence_list):
        for pam in pam_sequence_list:
            set_pam_sequences(pam)
            guide_rnas, orientation, pams, sequences = get_guide_rnas_and_orientation(ref_sequence, edit, substitution_position)
            if process_guide_rnas_for_pam(guide_rnas, sequences, ref_sequence_original, edited_sequence_original, df_dict, ref_sequence, edited_sequence, substitution_position, orientation, pam):
                return df_dict
        primedesign_output = run_prime_design(ref_sequence, edited_sequence, substitution_position)
        if primedesign_output != "No PrimeDesign Recommended Guides":
            update_df_dict_with_primedesign_output(df_dict, ref_sequence_original, edited_sequence_original, primedesign_output, "Prime Editing")
        else:
            add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "No Base or Prime Editing Guides Found")
        return df_dict

    def handle_insertion_or_deletion(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original):
        if '-' in ref_sequence:
            return handle_insertions(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original)
        else:
            return handle_deletions(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original)

    try:
        set_pam_sequences(PAM)
        ref_sequence, edited_sequence = create_mutable_sequences(ref_sequence_original, edited_sequence_original)
        df_dict = collections.defaultdict(list)

        if len(set(ref_sequence) - bases) == 0:
            substitution_position = identify_substitution_position(ref_sequence, edited_sequence)
            edit = f"{ref_sequence[substitution_position]}>{edited_sequence[substitution_position]}"
            pam_sequence_list = [PAM, 'NRN', 'NYN']
            df_dict = handle_base_editing(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original, substitution_position, pam_sequence_list)
        else:
            df_dict = handle_insertion_or_deletion(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original)
    except Exception as e:
        print(f"Error in get_guides: {e}")
        add_insertion_deletion_entries(df_dict, ref_sequence_original, edited_sequence_original, "No Base or Prime Editing Guides Found")

    return render_dataframe(df_dict)



# Helper functions
def gc_content(sequence):
    sequence = sequence.upper()
    GC_count = sequence.count('G') + sequence.count('C')
    GC_content = float(GC_count)/float(len(sequence))

    return("%.2f" % GC_content)

# IUPAC code map
iupac2bases_dict = {'A':'A','T':'T','C':'C','G':'G','a':'a','t':'t','c':'c','g':'g',
'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACTG]',
'r':'[ag]','y':'[ct]','s':'[gc]','w':'[at]','k':'[gt]','m':'[ac]','b':'[cgt]','d':'[agt]','h':'[act]','v':'[acg]','n':'[actg]',
'(':'(',')':')','+':'+','-':'-','/':'/'}

def iupac2bases(iupac):

    try:
        bases = iupac2bases_dict[iupac]
    except:
        logger.error('Symbol %s is not within the IUPAC nucleotide code ...' % str(iupac))
        sys.exit(1)

    return(bases)

# Reverse complement function
def reverse_complement(sequence):
    sequence = sequence
    new_sequence = ''
    for base in sequence:
        if base == 'A':
            new_sequence += 'T'
        elif base == 'T':
            new_sequence += 'A'
        elif base == 'C':
            new_sequence += 'G'
        elif base == 'G':
            new_sequence += 'C'
        elif base == 'a':
            new_sequence += 't'
        elif base == 't':
            new_sequence += 'a'
        elif base == 'c':
            new_sequence += 'g'
        elif base == 'g':
            new_sequence += 'c'
        elif base == '[':
            new_sequence += ']'
        elif base == ']':
            new_sequence += '['
        elif base == '+':
            new_sequence += '+'
        elif base == '-':
            new_sequence += '-'
        elif base == '/':
            new_sequence += '/'
        elif base == '(':
            new_sequence += ')'
        elif base == ')':
            new_sequence += '('
    return(new_sequence[::-1])

def almost_reverse_complement(sequence):
    sequence = sequence
    new_sequence = ''
    for base in sequence:
        if base == 'A':
            new_sequence += 'T'
        elif base == 'T':
            new_sequence += 'A'
        elif base == 'C':
            new_sequence += 'G'
        elif base == 'G':
            new_sequence += 'C'
        elif base == 'a':
            new_sequence += 't'
        elif base == 't':
            new_sequence += 'a'
        elif base == 'c':
            new_sequence += 'g'
        elif base == 'g':
            new_sequence += 'c'
        elif base == '[':
            new_sequence += ']'
        elif base == ']':
            new_sequence += '['
        elif base == '+':
            new_sequence += '+'
        elif base == '-':
            new_sequence += '-'
        elif base == '/':
            new_sequence += '/'
        elif base == '(':
            new_sequence += ')'
        elif base == ')':
            new_sequence += '('
    return new_sequence

# Extract reference and edited sequence information
def process_sequence(input_sequence):

    input_sequence = ''.join(input_sequence.split())

    # Check formatting is correct
    format_check = ''
    for i in input_sequence:
        if i == '(':
            format_check += '('
        elif i == ')':
            format_check += ')'
        elif i == '/':
            format_check += '/'
        elif i == '+':
            format_check += '+'
        elif i == '-':
            format_check += '-'

    # Check composition of input sequence
    if len(input_sequence) != sum([1 if x in ['A','T','C','G','(',')','+','-','/'] else 0 for x in input_sequence.upper()]):
        assert False

    # Check formatting
    if format_check.count('(') == format_check.count(')') and format_check.count('(') > 0: # Left and right parantheses equal
        if '((' not in format_check: # Checks both directions for nested parantheses
            if '()' not in format_check: # Checks for empty annotations
                if sum([1 if x in format_check else 0 for x in ['++','--','//','+-','+/','-+','-/','/+','/-','/(','+(','-(',')/',')+',')-']]) == 0:
                    pass

    # Create mapping between input format and reference and edit sequence
    editformat2sequence = {}
    edits = re.findall('\(.*?\)', input_sequence)
    for edit in edits:
        if '/' in edit:
            editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','')]
        elif '+' in edit:
            editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','')]
        elif '-' in edit:
            editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '']

    # Create mapping between edit number and reference and edit sequence
    editformat2sequence = {}
    editnumber2sequence = {}
    edit_idxs = [[m.start(), m.end()] for m in re.finditer('\(.*?\)', input_sequence)]
    edit_counter = 1
    for edit_idx in edit_idxs:
        edit = input_sequence[edit_idx[0]:edit_idx[1]]

        # Create edit format and number to sequence map
        if '/' in edit:
            editformat2sequence[edit] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower(), edit_counter]
            editnumber2sequence[edit_counter] = [edit.split('/')[0].replace('(',''), edit.split('/')[1].replace(')','').lower()]

        elif '+' in edit:
            editformat2sequence[edit] = ['' , edit.split('+')[1].replace(')','').lower(), edit_counter]
            editnumber2sequence[edit_counter] = ['' , edit.split('+')[1].replace(')','').lower()]

        elif '-' in edit:
            editformat2sequence[edit] = [edit.split('-')[1].replace(')',''), '', edit_counter]
            editnumber2sequence[edit_counter] = [edit.split('-')[1].replace(')',''), '']

        edit_counter += 1

    edit_start = min([i.start() for i in re.finditer('\(', input_sequence)])
    edit_stop = max([i.start() for i in re.finditer('\)', input_sequence)])

    edit_span_sequence_w_ref = input_sequence[edit_start:edit_stop + 1]
    edit_span_sequence_w_edit = input_sequence[edit_start:edit_stop + 1]
    for edit in editformat2sequence:
        edit_span_sequence_w_ref = edit_span_sequence_w_ref.replace(edit, editformat2sequence[edit][0])
        edit_span_sequence_w_edit = edit_span_sequence_w_edit.replace(edit, editformat2sequence[edit][1])

    edit_start_in_ref = re.search('\(', input_sequence).start()
    edit_stop_in_ref_rev = re.search('\)', input_sequence[::-1]).start()

    edit_span_length_w_ref = len(edit_span_sequence_w_ref)
    edit_span_length_w_edit = len(edit_span_sequence_w_edit)

    reference_sequence = input_sequence
    edit_sequence = input_sequence
    editnumber_sequence = input_sequence
    for edit in editformat2sequence:
        reference_sequence = reference_sequence.replace(edit, editformat2sequence[edit][0])
        edit_sequence = edit_sequence.replace(edit, editformat2sequence[edit][1])
        editnumber_sequence = editnumber_sequence.replace(edit, str(editformat2sequence[edit][2]))

    return(editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev)

def run_primedesign(input_sequence):
    target_design = {}
    peg_design = {'pegRNA group':[],'type':[], 'spacer sequence':[],'spacer GC content':[],'PAM':[],'strand':[],'peg-to-edit distance':[],'nick-to-peg distance':[],'pegRNA extension':[], 'extension first base':[],'PBS length':[],'PBS GC content':[], 'PBS Tm':[], 'RTT length':[],'RTT GC content':[],'annotation':[],'spacer top strand oligo':[], 'spacer bottom strand oligo':[], 'pegRNA extension top strand oligo':[], 'pegRNA extension bottom strand oligo':[], 'CFD score':[]}

    input_sequence = ''.join(input_sequence.split())
    pe_format = 'NNNNNNNNNNNNNNNNN/NNN[NGG]'
    
    pbs_length_list = list(range(12, 15))
    rtt_length_list = list(range(10, 21))

    if 80 not in rtt_length_list:
        rtt_length_list.append(80)
        rtt_length_list = sorted(rtt_length_list)

    nicking_distance_minimum = 0
    nicking_distance_maximum = 100

    target_sequence = input_sequence
    target_name = 'user-input'

    target_sequence = target_sequence.upper()
    editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence, editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit, edit_start_in_ref, edit_stop_in_ref_rev = process_sequence(target_sequence)

    # Initialize dictionary for the design of pegRNA spacers for each target sequence and intended edit(s)
    target_design[target_name] = {'target_sequence':target_sequence, 'editformat2sequence': editformat2sequence, 'editnumber2sequence': editnumber2sequence, 'reference_sequence': reference_sequence, 'edit_sequence': edit_sequence, 'editnumber_sequence': editnumber_sequence, 'edit_span_length': [edit_span_length_w_ref, edit_span_length_w_edit], 'edit_start_in_ref': edit_start_in_ref, 'edit_stop_in_ref_rev': edit_stop_in_ref_rev, 'pegRNA':{'+':[], '-':[]}, 'ngRNA':{'+':[], '-':[]}}

    # Find indices but shift when removing annotations
    cut_idx = re.search('/', pe_format).start()
    pam_start_idx = re.search('\[', pe_format).start()
    pam_end_idx = re.search('\]', pe_format).start()

    # Find pam and total PE format search length
    pam_length = pam_end_idx - pam_start_idx - 1
    pe_format_length = len(pe_format) - 3

    # Check if cut site is left of PAM
    if cut_idx < pam_start_idx:

        # Shift indices with removal of annotations
        pam_start_idx = pam_start_idx - 1
        pam_end_idx = pam_end_idx - 2
        spacer_start_idx = 0
        spacer_end_idx = pam_start_idx

    else:
        pam_end_idx = pam_end_idx - 1
        cut_idx = cut_idx - 2
        spacer_start_idx = pam_end_idx
        spacer_end_idx = len(pe_format) - 3

    # Remove annotations and convert into regex
    pe_format_rm_annotation = pe_format.replace('/', '').replace('[', '').replace(']', '')

    # Create PE format and PAM search sequences
    pe_format_search_plus = ''
    for base in pe_format_rm_annotation:
        pe_format_search_plus += iupac2bases(base)
    pe_format_search_minus = reverse_complement(pe_format_search_plus)

    pam_search = ''
    pam_sequence = pe_format_rm_annotation[pam_start_idx:pam_end_idx]
    for base in pam_sequence:
        pam_search += iupac2bases(base)

    ##### Initialize data storage for output
    counter = 1
    for target_name in target_design:

        # pegRNA spacer search for (+) and (-) strands with reference sequence
        reference_sequence = target_design[target_name]['reference_sequence']
        find_guides_ref_plus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_plus, reference_sequence, re.IGNORECASE)]
        find_guides_ref_minus = [[m.start()] for m in re.finditer('(?=%s)' % pe_format_search_minus, reference_sequence, re.IGNORECASE)]

        # pegRNA spacer search for (+) and (-) strands with edit number sequence
        editnumber_sequence = target_design[target_name]['editnumber_sequence']
        find_guides_editnumber_plus = [[m.start()] for m in re.finditer('(?=%s)' % pam_search.replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]
        find_guides_editnumber_minus = [[m.start()] for m in re.finditer('(?=%s)' % reverse_complement(pam_search).replace('[', '[123456789'), editnumber_sequence, re.IGNORECASE)]

        # Find pegRNA spacers targeting (+) strand
        if find_guides_ref_plus:

            for match in find_guides_ref_plus:

                # Extract matched sequences and annotate type of prime editing
                full_search = reference_sequence[match[0]:match[0] + pe_format_length]
                spacer_sequence = full_search[spacer_start_idx:spacer_end_idx]
                extension_core_sequence = full_search[:cut_idx]
                downstream_sequence_ref = full_search[cut_idx:]
                downstream_sequence_length = len(downstream_sequence_ref)
                pam_ref = full_search[pam_start_idx:pam_end_idx]

                # Check to see if the extended non target strand is conserved in the edited strand
                try:
                    extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
                    downstream_sequence_edit = edit_sequence[extension_core_end_idx:extension_core_end_idx + downstream_sequence_length]
                    pam_edit = edit_sequence[extension_core_start_idx:extension_core_start_idx + pe_format_length][pam_start_idx:pam_end_idx]
                    
                    ## Annotate pegRNA
                    # Check if PAM is mutated relative to reference sequence
                    if pam_ref == pam_edit.upper():
                        pe_annotate = 'PAM_intact'

                    else:
                        # Check to see if mutation disrupts degenerate base positions within PAM
                        if re.search(pam_search, pam_edit.upper()):
                            pe_annotate = 'PAM_intact'

                        else:
                            pe_annotate = 'PAM_disrupted'

                    # Store pegRNA spacer
                    nick_ref_idx = match[0] + cut_idx
                    nick_edit_idx = extension_core_start_idx + cut_idx
                    target_design[target_name]['pegRNA']['+'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

                except:
                    continue

        # Find pegRNA spacers targeting (-) strand
        if find_guides_ref_minus:

            for match in find_guides_ref_minus:

                # Extract matched sequences and annotate type of prime editing
                full_search = reference_sequence[match[0]:match[0] + pe_format_length]
                spacer_sequence = full_search[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
                extension_core_sequence = full_search[pe_format_length - cut_idx:]
                downstream_sequence_ref = full_search[:pe_format_length - cut_idx]
                downstream_sequence_length = len(downstream_sequence_ref)
                pam_ref = full_search[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

                # Check to see if the extended non target strand is conserved in the edited strand
                try:
                    extension_core_start_idx, extension_core_end_idx = re.search(extension_core_sequence, edit_sequence).start(), re.search(extension_core_sequence, edit_sequence).end()
                    downstream_sequence_edit = edit_sequence[extension_core_start_idx - downstream_sequence_length:extension_core_start_idx]
                    pam_edit = edit_sequence[extension_core_end_idx - pe_format_length:extension_core_end_idx][pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]
                    
                    ## Annotate pegRNA
                    # Check if PAM is mutated relative to reference sequence
                    if pam_ref == pam_edit.upper():
                        pe_annotate = 'PAM_intact'

                    else:
                        # Check to see if mutation disrupts degenerate base positions within PAM
                        if re.search(reverse_complement(pam_search), pam_edit.upper()):
                            pe_annotate = 'PAM_intact'

                        else:
                            pe_annotate = 'PAM_disrupted'

                    # Store pegRNA spacer
                    nick_ref_idx = match[0] + (pe_format_length - cut_idx)
                    nick_edit_idx = extension_core_start_idx - downstream_sequence_length + (pe_format_length - cut_idx)
                    target_design[target_name]['pegRNA']['-'].append([nick_ref_idx, nick_edit_idx, full_search, spacer_sequence, pam_ref, pam_edit, pe_annotate])

                except:
                    continue

        # Find ngRNA spacers targeting (+) strand
        if find_guides_editnumber_plus:

            for match in find_guides_editnumber_plus:

                # Extract matched sequences and annotate type of prime editing
                full_search = editnumber_sequence[:match[0] + pam_length]
                
                full_search2ref = full_search
                full_search2edit = full_search
                for edit_number in editnumber2sequence:
                    full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
                    full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

                if len(full_search2edit[-pe_format_length:]) == pe_format_length:

                    # Identify ngRNA sequence information from edit sequence
                    full_search_edit = full_search2edit[-pe_format_length:]
                    spacer_sequence_edit = full_search_edit[spacer_start_idx:spacer_end_idx]
                    pam_edit = full_search_edit[pam_start_idx:pam_end_idx]

                    # Use reference sequence to find nick index
                    full_search_ref = full_search2ref[-pe_format_length:]
                    spacer_sequence_ref = full_search_ref[spacer_start_idx:spacer_end_idx]
                    pam_ref = full_search_ref[pam_start_idx:pam_end_idx]

                    # Annotate ngRNA
                    if spacer_sequence_edit.upper() == spacer_sequence_ref.upper():
                        ng_annotate = 'PE3'
                    else:
                        if spacer_sequence_edit.upper()[-10:] == spacer_sequence_ref.upper()[-10:]:
                            ng_annotate = 'PE3b-nonseed'
                        else:
                            ng_annotate = 'PE3b-seed'

                    # Store ngRNA spacer
                    nick_ref_idx = re.search(full_search_ref, reference_sequence).end() - (pe_format_length - cut_idx)
                    nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                    nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                    target_design[target_name]['ngRNA']['+'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

        # Find ngRNA spacers targeting (-) strand
        if find_guides_editnumber_minus:

            for match in find_guides_editnumber_minus:

                # Extract matched sequences and annotate type of prime editing
                full_search = editnumber_sequence[match[0]:]
                
                full_search2ref = full_search
                full_search2edit = full_search
                for edit_number in editnumber2sequence:
                    full_search2ref = full_search2ref.replace(str(edit_number), editnumber2sequence[edit_number][0])
                    full_search2edit = full_search2edit.replace(str(edit_number), editnumber2sequence[edit_number][1])

                if len(full_search2edit[:pe_format_length]) == pe_format_length:

                    # Identify ngRNA sequence information from edit sequence
                    full_search_edit = full_search2edit[:pe_format_length]
                    spacer_sequence_edit = full_search_edit[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
                    pam_edit = full_search_edit[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

                    # Use reference sequence to find nick index
                    full_search_ref = full_search2ref[:pe_format_length]
                    spacer_sequence_ref = full_search_ref[pe_format_length - spacer_end_idx:pe_format_length - spacer_start_idx]
                    pam_ref = full_search_ref[pe_format_length - pam_end_idx:pe_format_length - pam_start_idx]

                    # Annotate ngRNA
                    if spacer_sequence_edit.upper() == spacer_sequence_ref.upper():
                        ng_annotate = 'PE3'
                    else:
                        if spacer_sequence_edit.upper()[:10] == spacer_sequence_ref.upper()[:10]:
                            ng_annotate = 'PE3b-nonseed'
                        else:
                            ng_annotate = 'PE3b-seed'

                    # Store ngRNA spacer
                    nick_ref_idx = re.search(full_search_ref, reference_sequence).start() + (pe_format_length - cut_idx)
                    nick_edit_start_idx = re.search(spacer_sequence_edit, edit_sequence).start()
                    nick_edit_end_idx = re.search(spacer_sequence_edit, edit_sequence).end()
                    target_design[target_name]['ngRNA']['-'].append([nick_ref_idx, nick_edit_start_idx, nick_edit_end_idx, full_search_edit, spacer_sequence_edit, pam_edit, ng_annotate])

        # Grab index information of edits to introduce to target sequence
        edit_start_in_ref = int(target_design[target_name]['edit_start_in_ref'])
        edit_stop_in_ref_rev = int(target_design[target_name]['edit_stop_in_ref_rev'])
        edit_span_length_w_ref = int(target_design[target_name]['edit_span_length'][0])
        edit_span_length_w_edit = int(target_design[target_name]['edit_span_length'][1])

        # Design pegRNAs targeting the (+) strand
        counter = 1
        counted = []
        for peg_plus in target_design[target_name]['pegRNA']['+']:

            pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_plus
            pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '+']))

            pe_annotate_constant = pe_annotate

            # See if pegRNA spacer can introduce all edits
            nick2edit_length = edit_start_in_ref - pe_nick_ref_idx
            if nick2edit_length >= 0:

                # Loop through RTT lengths
                silent_mutation_edit = ''
                for rtt_length in rtt_length_list:

                    # See if RT length can reach entire edit
                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                    if nick2lastedit_length < rtt_length:

                        # Loop through PBS lengths
                        for pbs_length in pbs_length_list:
                            pe_pam_ref_silent_mutation = ''

                            # Construct pegRNA extension to encode intended edit(s)
                            if ((len(edit_sequence) - pe_nick_edit_idx) - rtt_length) < 0:
                                rtt_length = len(edit_sequence) - pe_nick_edit_idx

                            # Patch for NGG PAMs - may need to build something more generalizable in the future
                            pegRNA_ext = reverse_complement(edit_sequence[pe_nick_edit_idx - pbs_length:pe_nick_edit_idx + rtt_length])

                            # Check to see if pegRNA extension is within input sequence
                            if len(pegRNA_ext) == (pbs_length + rtt_length):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('pegRNA')
                                peg_design['spacer sequence'].append(pe_spacer_sequence)
                                peg_design['spacer GC content'].append(gc_content(pe_spacer_sequence))

                                if pe_pam_ref_silent_mutation == '':
                                    peg_design['PAM'].append(pe_pam_ref)
                                else:
                                    peg_design['PAM'].append(pe_pam_ref_silent_mutation)

                                peg_design['strand'].append('+')
                                peg_design['peg-to-edit distance'].append(nick2lastedit_length)
                                peg_design['nick-to-peg distance'].append('')
                                peg_design['pegRNA extension'].append(pegRNA_ext)
                                peg_design['extension first base'].append(pegRNA_ext[0])
                                peg_design['PBS length'].append(pbs_length)
                                peg_design['PBS GC content'].append(gc_content(pegRNA_ext[rtt_length:]))
                                peg_design['PBS Tm'].append(mt.Tm_NN(pegRNA_ext[-pbs_length:], nn_table=mt.R_DNA_NN1))
                                peg_design['RTT length'].append(rtt_length)
                                peg_design['RTT GC content'].append(gc_content(pegRNA_ext[:rtt_length]))
                                peg_design['annotation'].append(pe_annotate)
                                peg_design['CFD score'].append('')

                                if pe_spacer_sequence[0] == 'G':
                                    peg_design['spacer top strand oligo'].append('cacc' + pe_spacer_sequence + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement(pe_spacer_sequence))

                                else:
                                    peg_design['spacer top strand oligo'].append('caccG' + pe_spacer_sequence + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + pe_spacer_sequence))

                                peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                counted.append(counter)

                # Create ngRNAs targeting (-) strand for (+) pegRNAs
                if counter in counted:
                    for ng_minus in target_design[target_name]['ngRNA']['-']:
                        ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_minus
                        nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                        if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                            peg_design['pegRNA group'].append(counter)
                            peg_design['type'].append('ngRNA')
                            peg_design['spacer sequence'].append(reverse_complement(ng_spacer_sequence_edit))
                            peg_design['spacer GC content'].append(gc_content(reverse_complement(ng_spacer_sequence_edit)))
                            peg_design['PAM'].append(reverse_complement(ng_pam_edit))
                            peg_design['strand'].append('-')
                            peg_design['peg-to-edit distance'].append('')
                            peg_design['nick-to-peg distance'].append(nick_distance)
                            peg_design['pegRNA extension'].append('')
                            peg_design['extension first base'].append('')
                            peg_design['PBS length'].append('')
                            peg_design['PBS GC content'].append('')
                            peg_design['PBS Tm'].append('')
                            peg_design['RTT length'].append('')
                            peg_design['RTT GC content'].append('')
                            peg_design['annotation'].append(ng_annotate)
                            peg_design['CFD score'].append('')

                            if reverse_complement(ng_spacer_sequence_edit)[0] == 'G':
                                peg_design['spacer top strand oligo'].append('cacc' + reverse_complement(ng_spacer_sequence_edit))
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement(reverse_complement(ng_spacer_sequence_edit)))

                            else:
                                peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(ng_spacer_sequence_edit))
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + reverse_complement(ng_spacer_sequence_edit)))

                            peg_design['pegRNA extension top strand oligo'].append('')
                            peg_design['pegRNA extension bottom strand oligo'].append('')

                    counter += 1

        # Design pegRNAs targeting the (-) strand
        for peg_minus in target_design[target_name]['pegRNA']['-']:

            pe_nick_ref_idx, pe_nick_edit_idx, pe_full_search, pe_spacer_sequence, pe_pam_ref, pe_pam_edit, pe_annotate = peg_minus
            pegid = '_'.join(map(str, [pe_nick_ref_idx, pe_spacer_sequence, pe_pam_ref, pe_annotate, '-']))

            pe_annotate_constant = pe_annotate

            # See if pegRNA spacer can introduce all edits
            nick2edit_length = edit_stop_in_ref_rev - (len(reference_sequence) - pe_nick_ref_idx)
            if nick2edit_length >= 0:

                # Loop through RTT lengths
                silent_mutation_edit = ''
                for rtt_length in rtt_length_list:

                    # See if RT length can reach entire edit
                    nick2lastedit_length = nick2edit_length + edit_span_length_w_edit
                    if nick2lastedit_length < rtt_length:

                        # Loop through PBS lengths
                        for pbs_length in pbs_length_list:
                            pe_pam_ref_silent_mutation = ''

                            # Construct pegRNA extension to encode intended edit(s)
                            # pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]
                            if (pe_nick_edit_idx - rtt_length) < 0:
                                rtt_length = pe_nick_edit_idx

                            # Patch for NGG PAMs - may need to build something more generalizable in the future
                            pegRNA_ext = edit_sequence[pe_nick_edit_idx - rtt_length:pe_nick_edit_idx + pbs_length]

                            # Check to see if pegRNA extension is within input sequence
                            if len(pegRNA_ext) == (pbs_length + rtt_length):

                                peg_design['pegRNA group'].append(counter)
                                peg_design['type'].append('pegRNA')
                                peg_design['spacer sequence'].append(reverse_complement(pe_spacer_sequence))
                                peg_design['spacer GC content'].append(gc_content(reverse_complement(pe_spacer_sequence)))

                                if pe_pam_ref_silent_mutation == '':
                                    peg_design['PAM'].append(reverse_complement(pe_pam_ref))
                                else:
                                    peg_design['PAM'].append(pe_pam_ref_silent_mutation)

                                peg_design['strand'].append('-')
                                peg_design['peg-to-edit distance'].append(nick2lastedit_length)
                                peg_design['nick-to-peg distance'].append('')
                                peg_design['pegRNA extension'].append(pegRNA_ext)
                                peg_design['extension first base'].append(pegRNA_ext[0])
                                peg_design['PBS length'].append(pbs_length)
                                peg_design['PBS GC content'].append(gc_content(pegRNA_ext[rtt_length:]))
                                peg_design['PBS Tm'].append(mt.Tm_NN(pegRNA_ext[-pbs_length:], nn_table=mt.R_DNA_NN1))
                                peg_design['RTT length'].append(rtt_length)
                                peg_design['RTT GC content'].append(gc_content(pegRNA_ext[:rtt_length]))
                                peg_design['annotation'].append(pe_annotate)
                                peg_design['CFD score'].append('')

                                if reverse_complement(pe_spacer_sequence)[0] == 'G':
                                    peg_design['spacer top strand oligo'].append('cacc' + reverse_complement(pe_spacer_sequence) + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement(reverse_complement(pe_spacer_sequence)))

                                else:
                                    peg_design['spacer top strand oligo'].append('caccG' + reverse_complement(pe_spacer_sequence) + 'gtttt')
                                    peg_design['spacer bottom strand oligo'].append('ctctaaaac' + reverse_complement('G' + reverse_complement(pe_spacer_sequence)))

                                peg_design['pegRNA extension top strand oligo'].append('gtgc' + pegRNA_ext)
                                peg_design['pegRNA extension bottom strand oligo'].append('aaaa' + reverse_complement(pegRNA_ext))

                                counted.append(counter)

                # Create ngRNAs targeting (+) strand for (-) pegRNAs
                if counter in counted:
                    for ng_plus in target_design[target_name]['ngRNA']['+']:
                        ng_nick_ref_idx, ng_edit_start_idx, ng_edit_end_idx, ng_full_search_edit, ng_spacer_sequence_edit, ng_pam_edit, ng_annotate = ng_plus
                        nick_distance = ng_nick_ref_idx - pe_nick_ref_idx

                        if (abs(nick_distance) >= nicking_distance_minimum) and (abs(nick_distance) <= nicking_distance_maximum):

                            peg_design['pegRNA group'].append(counter)
                            peg_design['type'].append('ngRNA')
                            peg_design['spacer sequence'].append(ng_spacer_sequence_edit)
                            peg_design['spacer GC content'].append(gc_content(ng_spacer_sequence_edit))
                            peg_design['PAM'].append(ng_pam_edit)
                            peg_design['strand'].append('+')
                            peg_design['peg-to-edit distance'].append('')
                            peg_design['nick-to-peg distance'].append(nick_distance)
                            peg_design['pegRNA extension'].append('')
                            peg_design['extension first base'].append('')
                            peg_design['PBS length'].append('')
                            peg_design['PBS GC content'].append('')
                            peg_design['PBS Tm'].append('')
                            peg_design['RTT length'].append('')
                            peg_design['RTT GC content'].append('')
                            peg_design['annotation'].append(ng_annotate)
                            peg_design['CFD score'].append('')

                            if ng_spacer_sequence_edit[0] == 'G':
                                peg_design['spacer top strand oligo'].append('cacc' + ng_spacer_sequence_edit)
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement(ng_spacer_sequence_edit))

                            else:
                                peg_design['spacer top strand oligo'].append('caccG' + ng_spacer_sequence_edit)
                                peg_design['spacer bottom strand oligo'].append('aaac' + reverse_complement('G' + ng_spacer_sequence_edit))

                            peg_design['pegRNA extension top strand oligo'].append('')
                            peg_design['pegRNA extension bottom strand oligo'].append('')

                    counter += 1

    df = pd.DataFrame.from_dict(peg_design)

    try:
        df = df[~df['spacer sequence'].str.contains('TTTT')]
    except:
        pass

    df_pegs = df[df['type'] == 'pegRNA']

    # Find recommended pegRNA
    if len(df_pegs.sort_values(['annotation', 'peg-to-edit distance'])['pegRNA group']) > 0:

        edit_effective_length = max([edit_span_length_w_ref, edit_span_length_w_edit])
        if edit_effective_length <= 1:
            homology_downstream_recommended = 9
        elif edit_effective_length <= 5:
            homology_downstream_recommended = 14
        elif edit_effective_length <= 10:
            homology_downstream_recommended = 19
        elif edit_effective_length <= 15:
            homology_downstream_recommended = 24
        else:
            homology_downstream_recommended = 34

        pegrna_group = df_pegs.sort_values(['annotation', 'peg-to-edit distance'])['pegRNA group'].values[0]
        rtt_length_recommended = min(df_pegs[df_pegs['pegRNA group'] == pegrna_group]['peg-to-edit distance']) + homology_downstream_recommended
        rtt_max = max(df_pegs[(df_pegs['pegRNA group'] == pegrna_group)]['RTT length'].values)

        # find recommended PBS
        df_pegs['recommended_PBS_Tm'] = abs(df_pegs['PBS Tm'] - 37) # optimal PBS Tm of 37C
        pbs_length_recommended = df_pegs[df_pegs['pegRNA group'] == pegrna_group].sort_values(['recommended_PBS_Tm'], ascending = [True])['PBS length'].values[0]

        # find recommended RTT
        extension_first_base = 'C'
        while (extension_first_base == 'C') and (rtt_length_recommended < rtt_max):
            rtt_length_recommended += 1
            extension_first_base = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][rtt_max-int(rtt_length_recommended):rtt_max][0]

        if extension_first_base != 'C':

            pbs_extension_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][rtt_max:]
            rtt_extension_max = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['pegRNA extension'].values[0][:rtt_max]
            extension_recommended = rtt_extension_max[-rtt_length_recommended:] + pbs_extension_recommended

            peg_spacer_top_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['spacer top strand oligo'].values[0]
            peg_spacer_bottom_recommended = df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['spacer bottom strand oligo'].values[0]
            peg_ext_top_recommended = 'gtgc' + extension_recommended
            peg_ext_bottom_recommended = 'aaaa' + reverse_complement(extension_recommended)

            peg_annotation_recommended = ' %s' % str(df_pegs[(df_pegs['pegRNA group'] == pegrna_group) & (df_pegs['PBS length'] == pbs_length_recommended) & (df_pegs['RTT length'] == rtt_max)]['annotation'].values[0]).replace('_', ' ')
            peg_pbs_recommended = '%s nt' % str(pbs_length_recommended)
            peg_rtt_recommended = '%s nt' % str(rtt_length_recommended)

            # Find recommended ngRNA
            df_ngs = df[(df['type'] == 'ngRNA') & (df['pegRNA group'] == pegrna_group)]
            df_ngs['optimal_distance'] = abs(abs(df_ngs['nick-to-peg distance']) - 75) # optimal ngRNA +/- 75 bp away
            df_ngs = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True]) # prioritize PE3b

            if len(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer top strand oligo']) > 0:

                ng_spacer_top_recommended = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer top strand oligo'].values[0]
                ng_spacer_bottom_recommended = df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['spacer bottom strand oligo'].values[0]
                ng_annotation_recommended = ' %s' % str(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['annotation'].values[0]).replace('_', ' ')
                ng_distance_recommended = ' %s bp' % str(df_ngs.sort_values(['annotation', 'optimal_distance'], ascending = [False, True])['nick-to-peg distance'].values[0])

            else:

                ng_spacer_top_recommended = 'n/a'
                ng_spacer_bottom_recommended = 'n/a'
                ng_annotation_recommended = 'n/a'
                ng_distance_recommended = 'n/a'

        else:

            peg_spacer_top_recommended = 'n/a'
            peg_spacer_bottom_recommended = 'n/a'
            peg_ext_top_recommended = 'n/a'
            peg_ext_bottom_recommended = 'n/a'

            peg_annotation_recommended = ' n/a'
            peg_pbs_recommended = ' n/a'
            peg_rtt_recommended = ' n/a'

            ng_spacer_top_recommended = 'n/a'
            ng_spacer_bottom_recommended = 'n/a'
            ng_annotation_recommended = ' n/a'
            ng_distance_recommended = ' n/a'

    else:

        peg_spacer_top_recommended = ''
        peg_spacer_bottom_recommended = ''
        peg_ext_top_recommended = ''
        peg_ext_bottom_recommended = ''

        peg_annotation_recommended = ''
        peg_pbs_recommended = ''
        peg_rtt_recommended = ''

        ng_spacer_top_recommended = ''
        ng_spacer_bottom_recommended = ''
        ng_annotation_recommended = ''
        ng_distance_recommended = ''

    '''
    # Filter dataframes
    df = df[(df['extension first base'] != 'C')]

    pbs_range_list = list(range(12, 14 + 1)) + ['']
    rtt_range_list = list(range(10, 20 + 1)) + ['']

    pegrna_groups_to_keep = list(set(df_pegs[df_pegs['RTT length'].isin(rtt_range_list)]['pegRNA group'].values))

    df = df[df['PBS length'].isin(pbs_range_list)]
    df_pegs = df_pegs[df_pegs['PBS length'].isin(pbs_range_list)]

    df = df[df['RTT length'].isin(rtt_range_list)]
    df_pegs = df_pegs[df_pegs['RTT length'].isin(rtt_range_list)]

    df = df[df['pegRNA group'].isin(pegrna_groups_to_keep)]
    df_pegs = df_pegs[df_pegs['pegRNA group'].isin(pegrna_groups_to_keep)]

    df.reset_index(drop=True, inplace=True)

    df_pegs = df_pegs[['pegRNA group','spacer sequence','PAM','strand','peg-to-edit distance','spacer GC content','annotation']].drop_duplicates()
    df_pegs = df_pegs.sort_values('peg-to-edit distance')
    df_pegs.reset_index(drop=True, inplace=True)
    #return df_pegs, df, peg_spacer_top_recommended, peg_spacer_bottom_recommended, peg_ext_top_recommended, peg_ext_bottom_recommended, peg_annotation_recommended, peg_pbs_recommended, peg_rtt_recommended, ng_spacer_top_recommended, ng_spacer_bottom_recommended, ng_annotation_recommended, ng_distance_recommended
    '''
    if peg_spacer_top_recommended == '' and peg_spacer_bottom_recommended == '' and peg_ext_top_recommended == '' and peg_ext_bottom_recommended == '' and peg_annotation_recommended == '' and peg_pbs_recommended == '' and peg_rtt_recommended == '' and ng_spacer_top_recommended == '' and ng_spacer_bottom_recommended == '' and ng_annotation_recommended == '' and ng_distance_recommended == '':
        return "No PrimeDesign Recommended Guides"
    elif peg_spacer_top_recommended == 'n/a' and peg_spacer_bottom_recommended == 'n/a' and peg_ext_top_recommended == 'n/a' and peg_ext_bottom_recommended == 'n/a' and peg_annotation_recommended == ' n/a' and peg_pbs_recommended == ' n/a' and peg_rtt_recommended == ' n/a' and ng_spacer_top_recommended == 'n/a' and ng_spacer_bottom_recommended == 'n/a' and ng_annotation_recommended == ' n/a' and ng_distance_recommended == ' n/a':
        return "No PrimeDesign Recommended Guides"
    else:
        return (peg_spacer_top_recommended, peg_spacer_bottom_recommended, peg_ext_top_recommended, peg_ext_bottom_recommended, peg_annotation_recommended, peg_pbs_recommended, peg_rtt_recommended, ng_spacer_top_recommended, ng_spacer_bottom_recommended, ng_annotation_recommended, ng_distance_recommended)