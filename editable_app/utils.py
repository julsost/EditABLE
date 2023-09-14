import collections
import pandas as pd
import regex as re
from Bio.Seq import MutableSeq
from Bio.SeqUtils import MeltingTemp as mt
from math import log

bases = {"A", "C", "G", "T"}

PAMs = re.compile("[A|T|G|C]G")
reverse_PAMs = re.compile("C[A|T|G|C]")
edit_start = 4
edit_end = 8
gRNA_size = 20
pam_length = 2


# Helper function to find guide RNAs for ABE and CBE editable mutations on forward/reverse strands
def find_BE_guide_rnas(direction, seq, genomic_location):
    if direction == "forward":
        start_pams = genomic_location + (gRNA_size - edit_end) + 1
        end_pams = genomic_location + (gRNA_size - edit_start) + pam_length + 1
        possible_pams = seq[start_pams : end_pams]
        if len(re.findall(PAMs, str(possible_pams))) == 0:
            return "NO PAM"
        else:
            guideRNAs = []
            for match in re.finditer(PAMs, str(possible_pams), overlapped=True):
                gRNA = seq[start_pams + match.start() - gRNA_size: start_pams + match.start()]
                guideRNAs.append(str(gRNA))
            return guideRNAs

    elif direction == "reverse":
        start_pams = genomic_location - gRNA_size + edit_start - pam_length
        end_pams = genomic_location - gRNA_size + edit_end
        possible_pams = seq[start_pams : end_pams]
        if len(re.findall(reverse_PAMs, str(possible_pams))) == 0:
            return "NO PAM"
        else:
            guideRNAs = []
            for match in re.finditer(reverse_PAMs, str(possible_pams), overlapped=True):
                gRNA = seq[start_pams + match.start() + pam_length : start_pams + match.start() + pam_length + gRNA_size]
                gRNA.reverse_complement(inplace=True)
                guideRNAs.append(str(gRNA))
            return guideRNAs
    else:
        return "ERROR"


# Helper function to find guide RNAs for transversion (C>G and G>C) mutations on forward/reverse strands
def find_trans_guide_rnas(direction, seq, genomic_location):
    if direction == "forward":
        pam = seq[genomic_location + 15 : genomic_location + 15 + pam_length]
        if len(re.findall(PAMs, str(pam))) == 0:
            return "NO PAM"
        else:
            return seq[genomic_location + 15 - gRNA_size : genomic_location + 15]
    elif direction == "reverse":
        pam = seq[genomic_location - 15 - 1: genomic_location - 15 - 1 + pam_length]
        if len(re.findall(reverse_PAMs, str(pam))) == 0:
            return "NO PAM"
        else:
            gRNA = seq[genomic_location - 15 - 1 + pam_length: genomic_location - 15 - 1 + pam_length + gRNA_size]
            gRNA.reverse_complement(inplace=True)
            return gRNA
    else:
        return "ERROR"


# Adds the base editable guide RNAs to the data table for export
def get_guide_RNAs(mutant_seq, edit_type, genomic_location):
    if edit_type == "A>G":
        return find_BE_guide_rnas("forward", mutant_seq, genomic_location), "forward"
    elif edit_type == "T>C":
        return find_BE_guide_rnas("reverse", mutant_seq, genomic_location), "reverse"
    elif edit_type == "G>A":
        return find_BE_guide_rnas("forward", mutant_seq, genomic_location), "forward"
    elif edit_type == "C>T":
        return find_BE_guide_rnas("reverse", mutant_seq, genomic_location), "reverse"
    elif edit_type == "G>C":
        return find_trans_guide_rnas("forward", mutant_seq, genomic_location), "forward"
    elif edit_type == "C>G":
        return find_trans_guide_rnas("reverse", mutant_seq, genomic_location), "reverse"
    else:
        return "NOT BASE EDITABLE", None


def get_guides(ref_sequence_original, edited_sequence_original):
    ref_sequence = MutableSeq(ref_sequence_original)
    edited_sequence = MutableSeq(edited_sequence_original)

    df_dict = collections.defaultdict(list)
    
    if len(set(ref_sequence) - bases) == 0:
            substitution_position = None
            for i in range(len(ref_sequence)):
                ref_base = ref_sequence[i]
                edited_base = edited_sequence[i]
                if ref_base != edited_base:
                    substitution_position = i
                    break
            edit = f"{ref_sequence[substitution_position]}>{edited_sequence[substitution_position]}"
            guide_rnas, orientation = get_guide_RNAs(ref_sequence, edit, substitution_position)
            if guide_rnas == "NOT BASE EDITABLE" or guide_rnas == "NO PAM":                
                primedesign_input = str(ref_sequence[:substitution_position] + f"({ref_sequence[substitution_position]}/{edited_sequence[substitution_position]})" + ref_sequence[substitution_position + 1:])
                primedesign_output = run_primedesign(str(primedesign_input))
                if primedesign_output != "No PrimeDesign Recommended Guides":
                    peg_spacer_top_recommended = primedesign_output[0]
                    peg_spacer_bottom_recommended = primedesign_output[1]
                    peg_ext_top_recommended = primedesign_output[2]
                    peg_ext_bottom_recommended = primedesign_output[3]
                    peg_annotation_recommended = primedesign_output[4]
                    peg_pbs_recommended = primedesign_output[5]
                    peg_rtt_recommended = primedesign_output[6]
                    ng_spacer_top_recommended = primedesign_output[7]
                    ng_spacer_bottom_recommended = primedesign_output[8]
                    ng_annotation_recommended = primedesign_output[9]
                    ng_distance_recommended = primedesign_output[10]
                    
                    df_dict['Original Sequence'].append(ref_sequence_original)
                    df_dict['Edited Sequence'].append(edited_sequence_original)
                    df_dict['Editing Technology'].append("Prime Editing")
                    df_dict['Base Editing Guide'].append(None)
                    
                    df_dict['PrimeDesign pegRNA Annotation'].append(peg_annotation_recommended)
                    df_dict['PrimeDesign pegRNA PBS'].append(peg_pbs_recommended)
                    df_dict['PrimeDesign pegRNA RTT'].append(peg_rtt_recommended)

                    df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(peg_spacer_top_recommended)
                    df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(peg_spacer_bottom_recommended)
                    df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(peg_ext_top_recommended)
                    df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(peg_ext_bottom_recommended)
                    
                    df_dict['PrimeDesign ngRNA Annotation'].append(ng_annotation_recommended)
                    df_dict['PrimeDesign ngRNA Distance'].append(ng_distance_recommended)

                    df_dict['PrimeDesign ngRNA Oligo Top'].append(ng_spacer_top_recommended)
                    df_dict['PrimeDesign ngRNA Bottom Top'].append(ng_spacer_bottom_recommended)
                else:
                    df_dict['Original Sequence'].append(ref_sequence_original)
                    df_dict['Edited Sequence'].append(edited_sequence_original)
                    df_dict['Editing Technology'].append("No Base or Prime Editing Guides Found")
                    df_dict['Base Editing Guide'].append(None)
                    
                    df_dict['PrimeDesign pegRNA Annotation'].append(None)
                    df_dict['PrimeDesign pegRNA PBS'].append(None)
                    df_dict['PrimeDesign pegRNA RTT'].append(None)

                    df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
                    df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
                    df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
                    df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
                    
                    df_dict['PrimeDesign ngRNA Annotation'].append(None)
                    df_dict['PrimeDesign ngRNA Distance'].append(None)

                    df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
                    df_dict['PrimeDesign ngRNA Bottom Top'].append(None)
            else:
                if not isinstance(guide_rnas, list):
                    guide_rnas = [guide_rnas]
                for gRNA in guide_rnas:
                    df_dict['Original Sequence'].append(ref_sequence_original)
                    df_dict['Edited Sequence'].append(edited_sequence_original)
                    df_dict['Editing Technology'].append("Base Editing")
                    df_dict['Base Editing Guide'].append(str(gRNA))
                    df_dict['Base Editing Guide Orientation'].append(orientation)

                    df_dict['PrimeDesign pegRNA Annotation'].append(None)
                    df_dict['PrimeDesign pegRNA PBS'].append(None)
                    df_dict['PrimeDesign pegRNA RTT'].append(None)

                    df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
                    df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
                    df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
                    df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
                    
                    df_dict['PrimeDesign ngRNA Annotation'].append(None)
                    df_dict['PrimeDesign ngRNA Distance'].append(None)

                    df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
                    df_dict['PrimeDesign ngRNA Bottom Top'].append(None)

    # If the length of the reference sequence is less than the length of the edited sequence, there has been either an insertion or a duplication
    elif '-' in ref_sequence:
        insertion_start = ref_sequence.find('-')
        insertion_end = ref_sequence.rfind('-')
        if insertion_end - insertion_start + 1 > 44:
            df_dict['Original Sequence'].append(ref_sequence_original)
            df_dict['Edited Sequence'].append(edited_sequence_original)
            df_dict['Editing Technology'].append("Use Twin Prime Editing/Integrase/HDR")
            df_dict['Base Editing Guide'].append(None)
            
            df_dict['PrimeDesign pegRNA Annotation'].append(None)
            df_dict['PrimeDesign pegRNA PBS'].append(None)
            df_dict['PrimeDesign pegRNA RTT'].append(None)

            df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
            df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
            df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
            df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
            
            df_dict['PrimeDesign ngRNA Annotation'].append(None)
            df_dict['PrimeDesign ngRNA Distance'].append(None)

            df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
            df_dict['PrimeDesign ngRNA Bottom Top'].append(None)
        else:
            primedesign_input = str(ref_sequence[:insertion_start] + f"(+{edited_sequence[insertion_start : insertion_end + 1]})" + ref_sequence[insertion_end + 1:])
            primedesign_output = run_primedesign(primedesign_input)
            if primedesign_output != "No PrimeDesign Recommended Guides":
                peg_spacer_top_recommended = primedesign_output[0]
                peg_spacer_bottom_recommended = primedesign_output[1]
                peg_ext_top_recommended = primedesign_output[2]
                peg_ext_bottom_recommended = primedesign_output[3]
                peg_annotation_recommended = primedesign_output[4]
                peg_pbs_recommended = primedesign_output[5]
                peg_rtt_recommended = primedesign_output[6]
                ng_spacer_top_recommended = primedesign_output[7]
                ng_spacer_bottom_recommended = primedesign_output[8]
                ng_annotation_recommended = primedesign_output[9]
                ng_distance_recommended = primedesign_output[10]

                df_dict['Original Sequence'].append(ref_sequence_original)
                df_dict['Edited Sequence'].append(edited_sequence_original)
                df_dict['Editing Technology'].append("Prime Editing")
                df_dict['Base Editing Guide'].append(None)
                
                df_dict['PrimeDesign pegRNA Annotation'].append(peg_annotation_recommended)
                df_dict['PrimeDesign pegRNA PBS'].append(peg_pbs_recommended)
                df_dict['PrimeDesign pegRNA RTT'].append(peg_rtt_recommended)

                df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(peg_spacer_top_recommended)
                df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(peg_spacer_bottom_recommended)
                df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(peg_ext_top_recommended)
                df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(peg_ext_bottom_recommended)
                
                df_dict['PrimeDesign ngRNA Annotation'].append(ng_annotation_recommended)
                df_dict['PrimeDesign ngRNA Distance'].append(ng_distance_recommended)

                df_dict['PrimeDesign ngRNA Oligo Top'].append(ng_spacer_top_recommended)
                df_dict['PrimeDesign ngRNA Bottom Top'].append(ng_spacer_bottom_recommended)
            else:
                df_dict['Original Sequence'].append(ref_sequence_original)
                df_dict['Edited Sequence'].append(edited_sequence_original)
                df_dict['Editing Technology'].append("Use Twin Prime Editing/Integrase/HDR")
                df_dict['Base Editing Guide'].append(None)
                
                df_dict['PrimeDesign pegRNA Annotation'].append(None)
                df_dict['PrimeDesign pegRNA PBS'].append(None)
                df_dict['PrimeDesign pegRNA RTT'].append(None)

                df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
                df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
                df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
                df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
                
                df_dict['PrimeDesign ngRNA Annotation'].append(None)
                df_dict['PrimeDesign ngRNA Distance'].append(None)

                df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
                df_dict['PrimeDesign ngRNA Bottom Top'].append(None)
            
    # If the length of the reference sequence is more than the length of the edited sequence, there has been a deletion
    elif '-' in edited_sequence:
        deletion_start = edited_sequence.find('-')
        deletion_end = edited_sequence.rfind('-')
        if deletion_start - deletion_end + 1 > 80:
            df_dict['Original Sequence'].append(ref_sequence_original)
            df_dict['Edited Sequence'].append(edited_sequence_original)
            df_dict['Editing Technology'].append("Use Twin Prime Editing/Integrase/HDR")
            df_dict['Base Editing Guide'].append(None)
            
            df_dict['PrimeDesign pegRNA Annotation'].append(None)
            df_dict['PrimeDesign pegRNA PBS'].append(None)
            df_dict['PrimeDesign pegRNA RTT'].append(None)

            df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
            df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
            df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
            df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
            
            df_dict['PrimeDesign ngRNA Annotation'].append(None)
            df_dict['PrimeDesign ngRNA Distance'].append(None)

            df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
            df_dict['PrimeDesign ngRNA Bottom Top'].append(None)
        else:
            primedesign_input = str(ref_sequence[:deletion_start] + f"(-{ref_sequence[deletion_start : deletion_end + 1]})" + ref_sequence[deletion_end + 1:])
            primedesign_output = run_primedesign(primedesign_input)
            if primedesign_output != "No PrimeDesign Recommended Guides":
                peg_spacer_top_recommended = primedesign_output[0]
                peg_spacer_bottom_recommended = primedesign_output[1]
                peg_ext_top_recommended = primedesign_output[2]
                peg_ext_bottom_recommended = primedesign_output[3]
                peg_annotation_recommended = primedesign_output[4]
                peg_pbs_recommended = primedesign_output[5]
                peg_rtt_recommended = primedesign_output[6]
                ng_spacer_top_recommended = primedesign_output[7]
                ng_spacer_bottom_recommended = primedesign_output[8]
                ng_annotation_recommended = primedesign_output[9]
                ng_distance_recommended = primedesign_output[10]

                df_dict['Original Sequence'].append(ref_sequence_original)
                df_dict['Edited Sequence'].append(edited_sequence_original)
                df_dict['Editing Technology'].append("Prime Editing")
                df_dict['Base Editing Guide'].append(None)
                
                df_dict['PrimeDesign pegRNA Annotation'].append(peg_annotation_recommended)
                df_dict['PrimeDesign pegRNA PBS'].append(peg_pbs_recommended)
                df_dict['PrimeDesign pegRNA RTT'].append(peg_rtt_recommended)

                df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(peg_spacer_top_recommended)
                df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(peg_spacer_bottom_recommended)
                df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(peg_ext_top_recommended)
                df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(peg_ext_bottom_recommended)
                
                df_dict['PrimeDesign ngRNA Annotation'].append(ng_annotation_recommended)
                df_dict['PrimeDesign ngRNA Distance'].append(ng_distance_recommended)

                df_dict['PrimeDesign ngRNA Oligo Top'].append(ng_spacer_top_recommended)
                df_dict['PrimeDesign ngRNA Bottom Top'].append(ng_spacer_bottom_recommended)
            else:
                df_dict['Original Sequence'].append(ref_sequence_original)
                df_dict['Edited Sequence'].append(edited_sequence_original)
                df_dict['Editing Technology'].append("Use Twin Prime Editing/Integrase/HDR")
                df_dict['Base Editing Guide'].append(None)
                
                df_dict['PrimeDesign pegRNA Annotation'].append(None)
                df_dict['PrimeDesign pegRNA PBS'].append(None)
                df_dict['PrimeDesign pegRNA RTT'].append(None)

                df_dict['PrimeDesign pegRNA Spacer Oligo Top'].append(None)
                df_dict['PrimeDesign pegRNA Spacer Oligo Bottom'].append(None)
                df_dict['PrimeDesign pegRNA Extension Oligo Top'].append(None)
                df_dict['PrimeDesign pegRNA Extension Oligo Bottom'].append(None)
                
                df_dict['PrimeDesign ngRNA Annotation'].append(None)
                df_dict['PrimeDesign ngRNA Distance'].append(None)

                df_dict['PrimeDesign ngRNA Oligo Top'].append(None)
                df_dict['PrimeDesign ngRNA Bottom Top'].append(None)

    
    df_dict_render = collections.defaultdict(list)
    for key in df_dict:
        if key != 'Base Editing Guide Orientation':
            for value in df_dict[key]:
                if key == 'Editing Technology' or key == 'PrimeDesign pegRNA Annotation':
                    df_dict_render[key].append(value)
                elif value is None:
                    df_dict_render[key].append(None)
                else:
                    key_len = len(key)
                    #value += len(value) % (key_len * "B"
                    new_value = '\n'.join(value[i:i+key_len] for i in range(0, len(value), key_len))
                    #new_value += (len(value) % key_len) * "*"
                    #new_value += "\n"
                    df_dict_render[key].append(new_value)
    
    return pd.DataFrame.from_dict(df_dict_render).dropna(how='all', axis=1), pd.DataFrame.from_dict(df_dict).dropna(how='all', axis=1)
    #return pd.DataFrame.from_dict(df_dict).dropna(how='all', axis=1)

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
