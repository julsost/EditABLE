import os
import re
import subprocess
from typing import Tuple, Union, Dict, List
import pandas as pd


def _parse_sg_output(raw: bytes) -> List[Dict]:
    """
    Parse Perl sgRNA output into a list of records.
    """
    text = raw.decode().strip().splitlines()
    header = text[0].split("\t")
    records: List[Dict] = []
    for line in text[1:]:
        vals = line.split("\t")
        rec = dict(zip(header, vals))
        rec["position"] = int(rec["position"])
        records.append(rec)
    return records


def _parse_alignment_output(raw: str) -> Dict[str, Union[str, int]]:
    """
    Parse Needleman-Wunsch alignment output into alignment fields.
    """
    parts = raw.strip().split("\t")
    return {
        "aligned_ref": parts[0],
        "aligned_ed": parts[1],
        "edit_start": int(parts[2]),
        "edit_end": int(parts[3]),
    }


def _parse_extension_output(raw: str) -> List[Dict]:
    """
    Parse PBS/RTT extension design output.
    """
    lines = raw.strip().splitlines()
    header = lines[0].split("\t")
    recs: List[Dict] = []
    for line in lines[1:]:
        vals = line.split("\t")
        rec = dict(zip(header, vals))
        # Convert numeric fields
        rec["RTT_len"] = int(rec["RTT_len"])
        rec["PBS_len"] = int(rec["PBS_len"])
        rec["GC"] = float(rec["GC"])
        rec["Tm"] = float(rec["Tm"])
        recs.append(rec)
    return recs


def _annotate_ng(spacer: str, ref: str, ed: str) -> str:
    """
    Simple ngRNA annotation: PE3b-seed if seed mismatch, else PE3.
    """
    if spacer not in ref and any(bs != es for bs, es in zip(spacer, ed)):
        return "PE3b-seed"
    return "PE3"


def _reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in seq[::-1])


def process_sequence(input_sequence: str) -> Tuple[Dict, Dict, str, str, str, int, int, int, int]:
    """
    Parse an annotated input string (with parentheses, +, -) into reference and edited sequences.
    Returns mappings, reference_sequence, edited_sequence, and some index info.
    """
    # Remove whitespace
    seq = ''.join(input_sequence.split())

    # Find annotations
    edits = re.findall(r'\(.*?\)', seq)
    editformat2sequence: Dict[str, List[Union[str, int]]] = {}
    editnumber2sequence: Dict[int, List[str]] = {}
    edit_idxs = [[m.start(), m.end()] for m in re.finditer(r'\(.*?\)', seq)]
    counter = 1
    for start, end in edit_idxs:
        edit = seq[start:end]
        if '/' in edit:
            a, b = edit[1:-1].split('/')
            editformat2sequence[edit] = [a, b, counter]
            editnumber2sequence[counter] = [a, b]
        elif '+' in edit:
            b = edit.split('+')[1][:-1]
            editformat2sequence[edit] = ['', b, counter]
            editnumber2sequence[counter] = ['', b]
        elif '-' in edit:
            a = edit.split('-')[1][:-1]
            editformat2sequence[edit] = [a, '', counter]
            editnumber2sequence[counter] = [a, '']
        counter += 1

    # Build reference and edited strings
    reference_sequence = seq
    edited_sequence = seq
    for key, vals in editformat2sequence.items():
        reference_sequence = reference_sequence.replace(key, vals[0])
        edited_sequence = edited_sequence.replace(key, vals[1])

    # Simple index info
    edit_start_in_ref = re.search(r'\(', seq).start()
    edit_stop_in_ref_rev = re.search(r'\)', seq[::-1]).start()
    edit_span_sequence_w_ref = ''.join(
        reference_sequence[edit_start_in_ref:edit_start_in_ref + len(edits[0])]
    )
    edit_span_length_w_ref = len(edit_span_sequence_w_ref)

    return (
        editformat2sequence,
        editnumber2sequence,
        reference_sequence,
        edited_sequence,
        '',  # editnumber_sequence unused here
        edit_span_length_w_ref,
        0,  # edit_span_length_w_edit unused
        edit_start_in_ref,
        edit_stop_in_ref_rev,
    )


def _run_pegfinder_engine(input_sequence: str) -> Union[Tuple[str, ...], str]:
    """
    Call the underlying Perl scripts to generate pegRNA/ngRNA pairs.
    """
    seq = re.sub(r"\s+", "", input_sequence).upper()
    # Parse reference vs edited
    _, _, reference_sequence, edited_sequence, _, _, _, _, _ = process_sequence(seq)

    base_dir = os.environ.get("PEGFINDER_DIR", os.getcwd())
    sgrna_primary = os.path.join(base_dir, "lib", "sub-sgRNA-finder.pl")
    sgrna_secondary = os.path.join(base_dir, "lib", "sub-sgRNA-finder.general.pl")
    nw_script = os.path.join(base_dir, "lib", "sub-needleman-wunsch.pl")
    pbs_rtt_script = os.path.join(base_dir, "lib", "sub-PBS-RT-design.pl")

    # 1. Primary pegRNA candidates
    peg_raw = subprocess.check_output([
        "perl", sgrna_primary, reference_sequence, "NGG"
    ])
    peg_candidates = _parse_sg_output(peg_raw)
    if not peg_candidates:
        return "No pegFinder Recommended Guides"

    # 2. Align reference vs edited
    nw_raw = subprocess.check_output([
        "perl", nw_script, reference_sequence, edited_sequence
    ]).decode()
    alignment = _parse_alignment_output(nw_raw)

    # 3. Design PBS & RTT
    ext_records: List[Dict] = []
    for peg in peg_candidates:
        ext_raw = subprocess.check_output([
            "perl",
            pbs_rtt_script,
            peg["spacer"],
            alignment["aligned_ref"],
            alignment["aligned_ed"],
        ]).decode()
        for ext in _parse_extension_output(ext_raw):
            rec = {
                **peg,
                **ext,
                "alignment_start": alignment["edit_start"],
                "alignment_end": alignment["edit_end"],
            }
            ext_records.append(rec)

    df_ext = pd.DataFrame(ext_records)
    df_ext = df_ext[~df_ext["spacer"].str.contains("TTTT")]
    if df_ext.empty:
        return "No pegFinder Recommended Guides"

    # 4. Secondary (nicking) sgRNAs
    ng_raw = subprocess.check_output([
        "perl", sgrna_secondary, edited_sequence, "NGG"
    ]).decode()
    ng_candidates = _parse_sg_output(ng_raw)

    # 5. Pairing
    pair_records: List[Dict] = []
    for _, peg_row in df_ext.iterrows():
        for ng in ng_candidates:
            dist = abs(ng["position"] - peg_row["position"])
            if 40 <= dist <= 150:
                ann = _annotate_ng(ng["spacer"], reference_sequence, edited_sequence)
                pair_records.append({
                    **peg_row.to_dict(),
                    "ng_spacer": ng["spacer"],
                    "ng_annotation": ann,
                    "ng_distance": dist,
                })

    if not pair_records:
        return "No pegFinder Recommended Guides"

    df_pairs = pd.DataFrame(pair_records)
    df_pairs["pam_intact_flag"] = df_pairs["pam"] == df_pairs["pam"]
    df_pairs["pbs_tm_diff"] = (df_pairs["Tm"] - 37).abs()
    df_pairs = df_pairs.sort_values(
        by=["alignment_start", "pam_intact_flag", "pbs_tm_diff", "ng_distance"]
    )
    best = df_pairs.iloc[0]

    # 6. Format output tuple
    return (
        best["spacer"],
        _reverse_complement(best["spacer"]),
        best["RTT_seq"] + best["PBS_seq"],
        _reverse_complement(best["RTT_seq"] + best["PBS_seq"]),
        "PAM intact" if best["pam_intact_flag"] else "PAM disrupted",
        f"{best['PBS_len']} nt",
        f"{best['RTT_len']} nt",
        best["ng_spacer"],
        _reverse_complement(best["ng_spacer"]),
        best["ng_annotation"],
        f"{best['ng_distance']} bp",
    )


def run_pegFinder(
    ref_seq: str,
    edited_seq: str,
    substitution_position: int
) -> Union[Tuple[str, ...], str]:
    """
    Build the pegFinder input string and call the engine.

    Args:
        ref_seq: reference DNA sequence
        edited_seq: edited DNA sequence
        substitution_position: index of base edit

    Returns:
        Either a tuple of recommended oligos, or the string
        "No pegFinder Recommended Guides".
    """
    peg_input = (
        ref_seq[:substitution_position]
        + f"({ref_seq[substitution_position]}/{edited_seq[substitution_position]})"
        + ref_seq[substitution_position + 1:]
    )
    return _run_pegfinder_engine(peg_input)
