import asyncio
import os
import pandas as pd
import re
from datetime import date, datetime
from pathlib import Path
from utils import get_guides, almost_reverse_complement, get_editor_info, process_guide_rnas, get_cloning_url, \
    process_ng_rnas, process_peg_rnas, visualization_peg_rnas, visualization_ng_rnas, handle_insertions
from shiny import App, render, ui, reactive
from shiny.types import ImgData



bases = {"A", "C", "G", "T"}
accepted_bases = {"A", "C", "G", "T", "-"}
#df_dict_result, bxb1_seq, modified_seq = handle_insertions(ref_sequence, edited_sequence, df_dict, ref_sequence_original, edited_sequence_original)
# A card component wrapper.
def ui_card(title, id, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body", "id": id}, *args),
        ),
    )


app_ui = ui.page_fluid(
    {"id": "main-content"},
    ui.output_image("display_logo", inline=True),
    ui.output_image("stanford_logo", inline=True),
    ui.br(),
    ui.help_text(
        '''Welcome to EditABLE!! We have designed this tool as a starting point for determining the optimal gene editing tool for a single gene edit. 
        Our algorithm is pending publication, but at a high level, it prioritizes base editing due to greater efficiency followed by prime editing. 
        Please refer to the following papers for more information on base and prime editing:'''
    ),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a(
        'Komor AC, Kim YB, Packer MS, Zuris JA, Liu DR. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. Nature. 2016;533(7603):420-424.',
        {'href': 'https://pubmed.ncbi.nlm.nih.gov/27096365/', 'target': '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a(
        'Gaudelli NM, Komor AC, Rees HA, et al. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage. Nature. 2017;551(7681):464-471.',
        {'href': 'https://pubmed.ncbi.nlm.nih.gov/29160308/', 'target': '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a(
        'Anzalone AV, Randolph PB, Davis JR, et al. Search-and-replace genome editing without double-strand breaks or donor DNA. Nature. 2019;576(7785):149-157.',
        {'href': 'https://pubmed.ncbi.nlm.nih.gov/31634902/', 'target': '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a(
        'Doman JL, Pandey S, Neugebauer ME, et al. Phage-assisted evolution and protein engineering yield compact, efficient prime editors. Cell. 2023;186(18):3983-4002.e26.',
        {'href': 'https://pubmed.ncbi.nlm.nih.gov/37657419/', 'target': '_blank'})),
    ui.br(),
    ui.br(),
    ui_card(
        "How to use this app",
        "how_to_app",
        ui.help_text(
            '''In CRISPR editing experiments, one is trying to induce some change in a DNA sequence.
             Therefore, you have an original sequence you are trying to change and a desired sequence 
             (what you want your sequence to look like after the CRISPR edit).'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("There are two ways to use this app:"),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''1. If you want to find guides for a single CRISPR edit. For this use case, please enter 
            in your original sequence and desired sequence in their respective input boxes.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''2. If you have more than one CRISPR edit you want to make, you can upload a CSV file with 
            two columns, named "Original Sequence" and "Desired Sequence" that contain your original and 
            desired sequences, with each row representing one edit you would like to make. Then click the 
            blue "Upload File" button even if the progress bar under the file browser says "Upload complete"'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''After specifing your input, click the "Find Guides" button and editABLE will try to find CRISPR 
            guides that will induce your desired edits. Once editABLE finishes running, a table will appear 
            displaying either the guides that editABLE has found for each of your desired edits or a suggestion to 
            use an alternative CRISPR technology if base or prime editing guides can't be found. Lastly, you can 
            download a CSV of the guides found by editABLE by clicking on the "Download Results as CSV File" 
            button. Base editing reagents will be suggested first due to their higher reported editing efficiency. 
            If a base editing guide cannot be found, we will then provide suggested prime editing reagents if possible.'''
        )
    ),
    ui_card(
        "Input requirements",
        "input_recs",
        ui.help_text(
            '''Your original sequence(s) and the desired sequence(s) must be the same length. Only single edits 
            (SNV, insertion, deletion) are supported at this time. Only the following characters are allowed in the input
            ("A", "C", "G", "T", "a", "c", "g", "t", "-"). All whitespace is allowed but will be removed before running our pipeline.
            Your sequences need to be from 5' to 3'.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For single nucleotide variant (SNV) edits, the input sequences are the most straightforward. You can 
            input the original and desired sequences without modification. For example, this would be a valid set of inputs:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Original Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("A", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Desired Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("G", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For changes that result in deletions, use a string of "-" characters in the desired sequence to denote the 
       GATAGCTCAGCTAGCCTAGTCAAACCTATCGACGTCGATCGATCGATCACACCGCCTAATC     area of the deletion. For example, this would be a valid set of inputs for a deletion:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Original Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("GGG", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Desired Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("---", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For changes that result in insertions/duplications, use a string of "-" characters in the original sequence 
            to denote the area of the insertion/duplication. For example, this would be a valid set of inputs for a insertions/duplications:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Original Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("---", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Desired Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"),
            ui.tags.b("GCG", style="color: red; font-family: Courier,courier"),
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            ui.tags.span(
                "Please limit the length of the insertion or deletion to less than 40 nucleotides as prime editors have been primarily tested in this context. Refer to  "),
            ui.tags.a("Anzalone et al. 2019 Nature", href="https://pubmed.ncbi.nlm.nih.gov/31634902/", target="_blank"),
            ui.tags.span(" for more information.")
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''Then, use the "Select Desired Base Editing PAM" dropdown to select the base editing PAM that is desired.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''Lastly, we require at least 25 base pairs of sequence to the left and right of your desired ''',
            ui.tags.b("edit", style="color: red"),
            '''. So in each of the examples above, there must be 25 or more base pairs to the right and left 
            of the ''',
            ui.tags.b("red", style="color: red"),
            ''' highlighted regions. '''
        )
    ),
    ui_card(
        "Input",
        'input',
        ui_card(
            "Single Sequence",
            'single_sequence_input',
            ui.help_text(
                "If you want to find guides for a single CRISPR edit. Please enter in your original sequence and desired sequence in their respective input boxes."
            ),
            ui.br(),
            ui.br(),
            ui.input_text_area("ref_sequence_input", "Original Sequence", placeholder="Enter sequence", height="50%",
                               width="100%"),
            ui.input_text_area("edited_sequence_input", "Desired Sequence", placeholder="Enter sequence", height="50%",
                               width="100%"),
            ui.input_action_button("find_guides_text", "Find Single Guides", class_="btn-primary"),
        ),
        ui_card(
            "Batch Sequence",
            'batch_sequence_input',
            ui.help_text(
                'If you have more than one CRISPR edit you want to make, you can upload a CSV file with two columns, named "Original Sequence" and "Desired Sequence" that contain your original and desired sequences, with each row representing one edit you would like to make.'
            ),
            ui.br(),
            ui.br(),
            ui.output_ui("ui_input_file"),
            ui.output_ui("upload_status"),
            ui.input_action_button("find_guides_csv", "Find Batch Guides", class_="btn-primary"),
        ),
        ui.div(
            ui.input_action_button("advanced_settings_toggle", "Advanced Base Editing Settings",
                                   class_="btn-secondary me-2"),
            ui.input_action_button("clear", "Clear Inputs", class_="btn-danger"),
            style="display: flex; gap: 5px; justify-content: left; margin-bottom: 30px;"
        ),
        ui.output_ui("advanced_settings")
    ),

    ui.output_ui("run_with_csv_input"),
    ui.output_ui("run_with_text_input"),
    ui.help_text(
        '''For troubleshooting and suggested revisions, please contact the ''',
        ui.tags.a("Bhalla Lab", {'href': 'https://med.stanford.edu/bhallalab.html', 'target': '_blank'}),
        " at vbhalla@stanford.edu"
    ),
    ui.br(),
    ui.br(),
    ui.br(),
    ui.br(),
)


# checks if no guides found
def check_for_BE_or_PE_guides(guides_df):
    no_BE_or_PE = guides_df[guides_df['Editing Technology'].isin(['No Base or Prime Editing Guides Found'])]
    if not no_BE_or_PE.empty:
        return "No Prime Editing or Base Editing Guide RNAs available for this site. In this case, we recommend the use of CRISPR with homologous recombination. See this paper here: "


# checks if pam defaults to nrn or nyn bc no guide rnas for selected pam
def check_availible_pam(PAM, selected_PAM):
    if PAM != selected_PAM:
        return f'No Guide RNAs available with {selected_PAM} PAM. See available {PAM} PAMs below:'


# checks if base editing window was changed for cgb edit
def check_CGB_window(window, input_window):
    if window != input_window:
        return 'Alternative editing windows are not available for C>G and G>C edits. The default 4 to 9 window editing is shown'


def check_ref_edited_pair(ref_sequence, edited_sequence, edit_start=4, edit_end=9):
    if len(ref_sequence) == 0 or len(edited_sequence) == 0:
        return False, "Both the original sequence and the edited sequence must be of nonzero length."
    if len(ref_sequence) != len(edited_sequence):
        return False, "The length of the original sequence and the edited sequence must be the same."
    if ref_sequence == edited_sequence:
        return False, "The original sequence and the edited sequence are the same."
    if len(set(ref_sequence) - accepted_bases) != 0 or len(set(edited_sequence) - accepted_bases) != 0:
        return False, "You may only have the following characters in your sequences {A, C, G, T, a, c, g, t, -}."
    if edit_start == edit_end:
        return False, "The base editing window must span multiple nucleotides."
    if edit_start > edit_end:
        return False, "The base editing window start cannot be later than the base editing window end."
    if edit_end > 20:
        return False, "The base editing window end can not be greater than 20."
    if edit_start < 1:
        return False, "The base editing window start can not be less than 1."
    if len(set(ref_sequence) - bases) == 0 and len(set(edited_sequence) - bases) == 0:
        substitution_position = None
        for i in range(len(ref_sequence)):
            ref_base = ref_sequence[i]
            edited_base = edited_sequence[i]
            if ref_base != edited_base:
                if substitution_position is not None:
                    return False, "The original sequence and the edited sequence contain more than one SNV. EditABLE currently only supports single SNVs, insertions, or deletions."
                else:
                    substitution_position = i
        if substitution_position < 25:
            return False, f"There must be at least 25 base pairs of sequence before the desired edit. {substitution_position} base pairs were found before your edit."
        if len(ref_sequence) - 1 - substitution_position < 25:
            return False, f"There must be at least 25 base pairs of sequence after the desired edit. {len(ref_sequence) - 1 - substitution_position} base pairs were found after your edit."
    else:
        if '-' not in ref_sequence and '-' not in edited_sequence:
            return False, 'The lengths of the original sequence and the edited sequence are not the same but neither has a "-" in it.'
        if '-' in ref_sequence and '-' in edited_sequence:
            return False, 'You cannot have a "-" in both the original and edited sequences.'
        elif '-' in ref_sequence:
            start_dash_position = None
            current_dash_position = None
            for i in range(len(ref_sequence)):
                if ref_sequence[i] == '-':
                    if start_dash_position is None:
                        start_dash_position = i
                    if current_dash_position is not None and i - current_dash_position != 1:
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple insertions. EditABLE currently only supports single SNVs, insertions, or deletions.'
                    else:
                        current_dash_position = i
            if start_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence before the desired edit. {start_dash_position} base pairs were found before your edit."
            if len(ref_sequence) - 1 - current_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence after the desired edit. {len(ref_sequence) - 1 - current_dash_position} base pairs were found after your edit."
        else:
            start_dash_position = None
            current_dash_position = None
            for i in range(len(edited_sequence)):
                if edited_sequence[i] == '-':
                    if start_dash_position is None:
                        start_dash_position = i
                    if current_dash_position is not None and i - current_dash_position != 1:
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple deletions. EditABLE currently only supports single SNVs, insertions, or deletions.'
                    else:
                        current_dash_position = i
            if start_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence before the desired edit. {start_dash_position} base pairs were found before your edit."
            if len(ref_sequence) - 1 - current_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence after the desired edit. {len(ref_sequence) - 1 - current_dash_position} base pairs were found after your edit."
    return True, "Inputs verified. Proceed to get guides."

#PLASMID CARDS

def create_prime_del_plasmid_card(guides_df, editor_info, editor_url):
    return ui_card(
        "Suggested Addgene Plasmids",
        "prime_del_plasmids",
        ui.help_text("Recommended pegRNA Plasmid A: Prime-Del Backbone (Addgene: 172657)"),
        ui.br(),
        ui.help_text("Recommended pegRNA Plasmid B: Prime-Del Backbone (Addgene: 172658)"),
        ui.br(),
        ui.help_text("Recommended Prime Editor Plasmid: PE2 (Addgene: 199267)"),
        ui.br(),
        ui.br(),
        ui.div(
            {"class": "d-flex"},
            ui.tags.a(
                "View Backbone A",
                href="http://n2t.net/addgene:172657",
                target="_blank",
                class_="btn btn-primary me-2"
            ),
            ui.tags.a(
                "View Backbone B",
                href="http://n2t.net/addgene:172658",
                target="_blank",
                class_="btn btn-primary me-2"
            ),
            ui.tags.a(
                "View Prime Editor Plasmid",
                href="http://n2t.net/addgene:199267",
                target="_blank",
                class_="btn btn-primary"
            ),
        ),
    )
def create_prime_editing_plasmid_card(guides_df, editor_info, editor_url):
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()

    processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligo_top,
                                               pegRNA_oligo_extension_bottom) if pegRNA_oligo_top and any(
        pegRNA_oligo_top) else None
    processed_ngRNA_oligos = process_ng_rnas(ngRNA_oligos) if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[
        0] != 'n/a' else None

    prime_card_elements = []

    # Add the recommended pegRNA plasmid information if pegRNA oligos are not empty
    if processed_pegRNA_oligos:
        prime_card_elements.append(
            ui.help_text("Recommended pegRNA Plasmid: pU6-tevopreq1-GG-acceptor (Addgene: 174038)"))
        prime_card_elements.append(ui.br())

    # Add the recommended ngRNA plasmid information if ngRNA oligos are not empty
    if processed_ngRNA_oligos:
        prime_card_elements.append(ui.help_text("Recommended ngRNA Plasmid: pmCherry-U6-empty (Addgene: 140580)"))
        prime_card_elements.append(ui.br())

    # Add the recommended prime editing plasmid information
    prime_card_elements.append(ui.help_text(f"Recommended Prime Editor Plasmid: {editor_info}"))
    prime_card_elements.append(ui.br())
    prime_card_elements.append(ui.br())

    # Add buttons at the bottom next to each other
    prime_card_elements.append(
        ui.div(
            {"class": "d-flex"},
            ui.tags.a("View pegRNA Plasmid", href="https://www.addgene.org/174038/", target="_blank",
                      class_="btn btn-primary me-2") if processed_pegRNA_oligos else None,
            ui.tags.a("View ngRNA Plasmid", href="https://www.addgene.org/65777/", target="_blank",
                      class_="btn btn-primary me-2") if processed_ngRNA_oligos else None,
            ui.tags.a("View Prime Editor Plasmid", href=editor_url, target="_blank", class_="btn btn-primary"),
        )
    )

    return ui_card(
        "Suggested Addgene Plasmids",
        "recommended_prime_editor",
        *[element for element in prime_card_elements if element is not None]  # Filter out None elements
    )

def create_twin_prime_editing_plasmid_card(guides_df, editor_info, editor_url):
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()

    processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligo_top,
                                               pegRNA_oligo_extension_bottom) if pegRNA_oligo_top and any(
        pegRNA_oligo_top) else None
    processed_ngRNA_oligos = process_ng_rnas(ngRNA_oligos) if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[
        0] != 'n/a' else None

    prime_card_elements = []

    # Add the recommended pegRNA plasmid information if pegRNA oligos are not empty
    if processed_pegRNA_oligos:
        prime_card_elements.append(
            ui.help_text("Recommended pegRNA Plasmid: pU6-tevopreq1-GG-acceptor (Addgene: 174038)"))
        prime_card_elements.append(ui.br())

    # Add the recommended ngRNA plasmid information if ngRNA oligos are not empty
    if processed_ngRNA_oligos:
        prime_card_elements.append(ui.help_text("Recommended ngRNA Plasmid: pmCherry-U6-empty (Addgene: 140580)"))
        prime_card_elements.append(ui.br())

    # Add the recommended prime editing plasmid information
    prime_card_elements.append(ui.help_text(f"Recommended Prime Editor Plasmid for creating Bxb1 integration site: {editor_info}"))
    prime_card_elements.append(ui.br())
    prime_card_elements.append(ui.br())

    # Add buttons at the bottom next to each other
    prime_card_elements.append(
        ui.div(
            {"class": "d-flex"},
            ui.tags.a("View pegRNA Plasmid", href="https://www.addgene.org/174038/", target="_blank",
                      class_="btn btn-primary me-2") if processed_pegRNA_oligos else None,
            ui.tags.a("View ngRNA Plasmid", href="https://www.addgene.org/65777/", target="_blank",
                      class_="btn btn-primary me-2") if processed_ngRNA_oligos else None,
            ui.tags.a("View Prime Editor Plasmid", href=editor_url, target="_blank", class_="btn btn-primary"),
        )
    )

    return ui_card(
        "Suggested Addgene Plasmids - Bxb1 integration step",
        "recommended_prime_editor",
        *[element for element in prime_card_elements if element is not None]  # Filter out None elements
    )

def create_integrase_plasmid_card(guides_df, editor_info, editor_url):
    return ui_card(
        "Suggested Addgene Plasmids for Bxb1 Integration",
        "integrase_only_card",
        ui.help_text("You’ll need both of these plasmids to perform Bxb1-mediated site-specific integration."),
        ui.br(),
        ui.br(),
        ui.div(
            {"class": "d-flex gap-2"},
            ui.tags.a(
                "View Bxb1 Integrase Plasmid",
                href="https://www.addgene.org/51271/",
                target="_blank",
                class_="btn btn-primary"
            ),
            ui.tags.a(
                "View attP Donor Vector (Clone Your Insert)",
                href="https://www.addgene.org/60561/",
                target="_blank",
                class_="btn btn-primary"
            ),
        )
    )
#-------------------------------------------------------------------------------------

#PROTOCOL/VALIDATION SECTIONS
def generate_prime_protocols_section(guides_df):
    # Extract the oligos from the DataFrame
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()

    # Initialize the prime section UI elements
    prime_section = [
        ui.help_text(
            ui.tags.b("pegRNA:", style="text-decoration: bold"),
            ui.br(),
            ui.tags.span("1. Order the following pegRNA as a geneBlock from "),
            ui.tags.a("IDT", href="https://www.idtdna.com", target="_blank"),
            ui.tags.span(" (or other preferred vendor):")
        ),
        ui.br(),
        ui.br(),
    ]

    # Process and add pegRNA oligos if they are present and valid
    if pegRNA_oligo_top and any(pegRNA_oligo_top):
        processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligo_top, pegRNA_oligo_extension_bottom)
        if processed_pegRNA_oligos:
            for guide_rna in processed_pegRNA_oligos:
                guide_rna_parts = guide_rna.split('\n')
                for part in guide_rna_parts:
                    prime_section.append(ui.help_text(part))
                    prime_section.append(ui.br())
                prime_section.append(ui.br())

    prime_section.extend([
        ui.help_text(
            ui.tags.span("2. Digest geneBlock using "),
            ui.tags.a("BsaI enzyme", href="https://www.neb.com/en-us/products/r3733-bsai-hf-v2", target="_blank"),
            ui.tags.span(" (37oC for 1 hr, 80 oC for 20 min) and purify using gel extraction kit "),
            ui.tags.a("(Zymo: D4007).",
                      href="https://www.zymoresearch.com/collections/zymoclean-gel-dna-recovery-kits/products/zymoclean-gel-dna-recovery-kit",
                      target="_blank"),
        ),
    ])
    prime_section.extend([
        ui.help_text(
            ui.br(),
            ui.br(),
            ui.tags.span("3. Follow Cloning Protocol here for "),
            ui.tags.a("pegRNA plasmid.", href="https://drive.google.com/file/d/1kKD7EVwS7nZMbiZ6UH4LftK7GepoS8si/view",
                      target="_blank"),
        ),
    ])
    # Process and add ngRNA oligos if they are present and valid
    if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[0] != 'n/a':
        processed_ngRNA_oligos = process_ng_rnas(ngRNA_oligos)
        if processed_ngRNA_oligos:
            prime_section.append(ui.help_text(
                ui.br(),
                ui.br(),
                ui.tags.b("ngRNA:", style="text-decoration: bold"),
                ui.br(),
                ui.tags.span("1. Order the following paired ngRNA Oligos from "),
                ui.tags.a("IDT", href="https://www.idtdna.com", target="_blank"),
                ui.tags.span(" (or other preferred vendor):")
            ))
            # displaying the ng rna oligo top and bottom
            prime_section.append(ui.br())
            prime_section.append(ui.br())
            for guide_rna in processed_ngRNA_oligos:
                guide_rna_parts = guide_rna.split('\n')
                for part in guide_rna_parts:
                    prime_section.append(ui.help_text(part))
                    prime_section.append(ui.br())
                prime_section.append(ui.br())

                # adding cloning protocal for ng rna
                prime_section.extend([
                    ui.help_text(
                        ui.tags.span("2. Follow Cloning Protocol here for "),
                        ui.tags.a("ngRNA plasmid.",
                                  href="https://drive.google.com/file/d/1RBkiXhl7CCKxHuDNkaogvjkEykxkH4ws/view?usp=sharing",
                                  target="_blank"),
                    ),
                ])

    prime_section.append(ui.br())
    prime_section.append(ui.br())

    return ui_card("Experimental Validation of Prime Editing", 'prime_section', *prime_section)
def generate_twin_prime_protocols_section(guides_df):
    # Extract the oligos from the DataFrame
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()

    # Initialize the prime section UI elements
    prime_section = [
        ui.help_text(
            ui.tags.b("pegRNA:", style="text-decoration: bold"),
            ui.br(),
            ui.tags.span("1. Order the following pegRNA as a geneBlock from "),
            ui.tags.a("IDT", href="https://www.idtdna.com", target="_blank"),
            ui.tags.span(" (or other preferred vendor):")
        ),
        ui.br(),
        ui.br(),
    ]

    # Process and add pegRNA oligos if they are present and valid
    if pegRNA_oligo_top and any(pegRNA_oligo_top):
        processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligo_top, pegRNA_oligo_extension_bottom)
        print("processed_pegRNA_oligos")
        if processed_pegRNA_oligos:
            for guide_rna in processed_pegRNA_oligos:
                guide_rna_parts = guide_rna.split('\n')
                for part in guide_rna_parts:
                    prime_section.append(ui.help_text(part))
                    prime_section.append(ui.br())
                prime_section.append(ui.br())
    prime_section.extend([
        ui.help_text(
            ui.tags.span("2. Digest geneBlock using "),
            ui.tags.a("BsaI enzyme", href="https://www.neb.com/en-us/products/r3733-bsai-hf-v2", target="_blank"),
            ui.tags.span(" (37oC for 1 hr, 80 oC for 20 min) and purify using gel extraction kit "),
            ui.tags.a("(Zymo: D4007).",
                      href="https://www.zymoresearch.com/collections/zymoclean-gel-dna-recovery-kits/products/zymoclean-gel-dna-recovery-kit",
                      target="_blank"),
        ),
    ])
    prime_section.extend([
        ui.help_text(
            ui.br(),
            ui.br(),
            ui.tags.span("3. Follow Cloning Protocol here for "),
            ui.tags.a("pegRNA plasmid.", href="https://drive.google.com/file/d/1kKD7EVwS7nZMbiZ6UH4LftK7GepoS8si/view",
                      target="_blank"),
        ),
    ])
    # Process and add ngRNA oligos if they are present and valid
    if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[0] != 'n/a':
        processed_ngRNA_oligos = process_ng_rnas(ngRNA_oligos)
        if processed_ngRNA_oligos:
            prime_section.append(ui.help_text(
                ui.br(),
                ui.br(),
                ui.tags.b("ngRNA:", style="text-decoration: bold"),
                ui.br(),
                ui.tags.span("1. Order the following paired ngRNA Oligos from "),
                ui.tags.a("IDT", href="https://www.idtdna.com", target="_blank"),
                ui.tags.span(" (or other preferred vendor):")
            ))
            # displaying the ng rna oligo top and bottom
            prime_section.append(ui.br())
            prime_section.append(ui.br())
            for guide_rna in processed_ngRNA_oligos:
                guide_rna_parts = guide_rna.split('\n')
                for part in guide_rna_parts:
                    prime_section.append(ui.help_text(part))
                    prime_section.append(ui.br())
                prime_section.append(ui.br())

                # adding cloning protocal for ng rna
                prime_section.extend([
                    ui.help_text(
                        ui.tags.span("2. Follow Cloning Protocol here for "),
                        ui.tags.a("ngRNA plasmid.",
                                  href="https://drive.google.com/file/d/1RBkiXhl7CCKxHuDNkaogvjkEykxkH4ws/view?usp=sharing",
                                  target="_blank"),
                    ),
                ])

    prime_section.append(ui.br())
    prime_section.append(ui.br())

    return ui_card("Experimental Validation of Prime Editing-Bxb1 integration step", 'prime_section', *prime_section)
def generate_integrase_protocols_section(guides_df):
    integrase_section = [
        ui.help_text(
            ui.tags.b("Integrase-Mediated Insertion Workflow:", style="text-decoration: bold"),
            ui.br(),
            ui.tags.span("Step 1: Prepare your plasmids"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Donor plasmid containing your insert flanked by the attB site"),
                ui.tags.li("Acceptor plasmid containing the attP site"),
                ui.tags.li("Bxb1 integrase expression plasmid (e.g., pCAG-NLS-HA-Bxb1, Addgene #51271)"),
            ),
            ui.br(),
            ui.tags.span("Step 2: Transfection into Mammalian Cells"),
            ui.br(),
            ui.tags.span("Seed mammalian cells (e.g., HEK293T) in 6-well plates at ~60-70% confluency one day before transfection."),
            ui.br(),
            ui.tags.span("Prepare DNA transfection mix per well:"),
            ui.tags.table(
                ui.tags.tr(
                    ui.tags.th("Plasmid"), ui.tags.th("Amount (ng)")
                ),
                ui.tags.tr(
                    ui.tags.td("Donor plasmid (insert + attB)"), ui.tags.td("500")
                ),
                ui.tags.tr(
                    ui.tags.td("Acceptor plasmid (attP)"), ui.tags.td("500")
                ),
                ui.tags.tr(
                    ui.tags.td("pCAG-NLS-HA-Bxb1 (integrase)"), ui.tags.td("500")
                ),
            ),
            ui.br(),
            ui.tags.span("Use a transfection reagent (e.g., Lipofectamine 3000):"),
            ui.tags.ul(
                ui.tags.li("Mix DNA with reagent according to manufacturer instructions."),
                ui.tags.li("Add complexes to cells."),
                ui.tags.li("Incubate cells 24–48 hours to allow integrase expression and recombination."),
            ),
            ui.br(),
            ui.tags.span("Step 3: Selection (if applicable)"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("If donor or acceptor plasmids carry antibiotic resistance, start selection 48 hours post-transfection."),
                ui.tags.li("Use appropriate antibiotic (e.g., puromycin, hygromycin)."),
                ui.tags.li("Continue selection for 5–7 days to enrich for recombined cells."),
                ui.tags.li("Alternatively, if fluorescent or other reporters are present, sort positive cells by FACS."),
            ),
            ui.br(),
            ui.tags.span("Step 4: Validation of Integration"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Extract genomic or plasmid DNA from selected cells."),
                ui.tags.li("Perform PCR across attL/attR junctions formed after recombination."),
                ui.tags.li("Sequence PCR products to verify correct integration."),
            ),
            ui.br(),
            ui.tags.span("Optional: Removal of attL/attR Sites"),
            ui.br(),
            ui.tags.span(
                "Optionally, you can remove the attL and attR recombination sites generated after integration using the ",
                ui.tags.b("EditABLE"), " tool for further sequence refinement and editing."
            ),
            ui.br(),
            ui.br(),
        )
    ]

    return ui_card("Detailed Bxb1 Integrase-Mediated Insertion Protocol", 'integrase_protocol_section', *integrase_section)

def generate_prime_del_protocols_section(guides_df):
    prime_del_section = [
        ui.help_text(
            ui.tags.b("PRIME-del Workflow:", style="text-decoration: bold"),
            ui.br(),

            # Components to Order
            ui.tags.span("Components to Order:"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("PE2 Prime Editor plasmid (Addgene #199267)"),
                ui.tags.li("pegRNA backbones (Addgene #172657, #172658)"),
                ui.tags.li("Custom pegRNA oligos or gene blocks (from EditABLE output)"),
            ),
            ui.br(),

            # Step 1
            ui.tags.span("Step 1: Clone pegRNAs into Backbone"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Gene block method (recommended for pegRNAs with PBS + RTT):"),
                ui.tags.ul(
                    ui.tags.li(
                        "Order pegRNA sequences as synthetic gene blocks (IDT, Twist, etc.) with flanking BsaI sites."),
                    ui.tags.li(
                        "Digest both gene block and pegRNA backbone plasmid with BsaI-HFv2 (NEB) at 37 °C for 1 hr."),
                    ui.tags.li("Heat inactivate at 80 °C for 20 min."),
                    ui.tags.li("Gel purify digested products (e.g., Zymo D4007 kit)."),
                    ui.tags.li("Ligate insert into backbone with T4 DNA ligase."),
                ),
                ui.br(),
                ui.tags.li("Oligo annealing method (for short pegRNA spacers ≤24 nt):"),
                ui.tags.ul(
                    ui.tags.li("Order paired oligos with EditABLE-provided CACC/AAAC overhangs."),
                    ui.tags.li("Anneal oligos: 95 °C for 5 min, then cool to 25 °C."),
                    ui.tags.li("Ligate into BsaI-digested backbone plasmid."),
                ),
                ui.br(),
                ui.tags.li("Transformation & Screening:"),
                ui.tags.ul(
                    ui.tags.li("Transform ligation into competent E. coli (DH5α)."),
                    ui.tags.li("Plate on LB agar + antibiotic (per plasmid resistance)."),
                    ui.tags.li("Screen colonies by colony PCR or direct sequencing."),
                ),
            ),
            ui.br(),

            # Step 2
            ui.tags.span("Step 2: Cell Culture and Transfection"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li(
                    "Plate cells (e.g., HEK293T, U2OS, or target line) 24 hrs before transfection to reach ~70% confluency."),
                ui.tags.li("Transfect with:"),
                ui.tags.ul(
                    ui.tags.li("PE2 plasmid (Addgene #199267): 1–2 µg per well in a 6-well plate"),
                    ui.tags.li("Paired pegRNA plasmids: 0.5–1 µg each"),
                ),
                ui.tags.li("Use Lipofectamine 3000 (ThermoFisher) or electroporation for harder-to-transfect cells."),
                ui.tags.li("Include a control well (PE2 + no pegRNAs)."),
            ),
            ui.br(),

            # Step 3
            ui.tags.span("Step 3: Expression and Editing"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Incubate cells for 72–96 hours to allow editing."),
                ui.tags.li("For difficult loci, extend incubation up to 7 days."),
                ui.tags.li("For stable expression, consider lentiviral delivery of pegRNAs."),
            ),
            ui.br(),

            # Step 4
            ui.tags.span("Step 4: Harvest and DNA Analysis"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Extract genomic DNA from cells."),
                ui.tags.li("PCR across the target locus using primers flanking the expected deletion."),
                ui.tags.li("Run PCR products on agarose gel: deletion allele runs smaller than wild-type."),
                ui.tags.li("Sequence PCR products (Sanger or NGS) to confirm precise deletion junction."),
            ),
            ui.br(),
        )
    ]

    return ui_card("PRIME-del Protocol", 'prime_del_protocol_section', prime_del_section)


def generate_donor_protocols_section(guides_df):
    integrase_section = [
        ui.help_text(
            ui.tags.b("Donor Plasmid Construction Protocol (Insert + attP Site):", style="text-decoration: bold"),
            ui.br(),

            ui.tags.span("Part 1: Prepare Your Insert"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Design your insert: include your gene or sequence of interest."),
                ui.tags.li(
                    "Include regulatory elements if needed: Promoter (e.g., CMV, EF1a), Kozak sequence, and PolyA tail."),
                ui.tags.li("Add flanking restriction sites (e.g., EcoRI/BamHI) or Gibson overlaps."),
            ),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Amplify your insert via high-fidelity PCR, or order it as a synthetic gBlock."),
                ui.tags.li("Run the PCR product or gBlock on an agarose gel."),
                ui.tags.li("Excise and purify using a gel extraction kit."),
            ),
            ui.br(),

            ui.tags.span("Part 2: Prepare the Vector"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Use pLenti-attP-DEST (Addgene #60561) or another attP-containing vector."),
                ui.tags.li(
                    "Digest 1–2 µg of the vector using two restriction enzymes compatible with your insert (e.g., EcoRI and BamHI)."),
                ui.tags.li("(Optional) Treat the digested vector with alkaline phosphatase to prevent re-ligation."),
                ui.tags.li("Run the digested vector on an agarose gel and purify the linearized backbone."),
            ),
            ui.br(),

            ui.tags.span("Part 3: Clone the Insert into the Vector"),
            ui.br(),
            ui.tags.span("Option A: Restriction-Ligation Cloning"),
            ui.tags.ul(
                ui.tags.li("Set up ligation reaction using T4 DNA ligase with a 3:1 molar ratio of insert to vector."),
                ui.tags.li("Incubate at 16°C overnight or 1–2 hours at room temperature."),
                ui.tags.li("Transform 5–10 µL of the ligation mix into competent E. coli (e.g., DH5α)."),
                ui.tags.li("Plate on LB agar with appropriate antibiotic (e.g., Ampicillin)."),
                ui.tags.li("Incubate overnight at 37°C."),
            ),
            ui.br(),
            ui.tags.span("Option B: Gibson Assembly"),
            ui.tags.ul(
                ui.tags.li("Design insert and vector fragments with 20–40 bp overlaps."),
                ui.tags.li("Mix insert and linearized vector with Gibson Assembly mix."),
                ui.tags.li("Incubate at 50°C for 1 hour."),
                ui.tags.li("Transform into competent E. coli and plate as above."),
            ),
            ui.br(),

            ui.tags.span("Part 4: Screen and Confirm"),
            ui.br(),
            ui.tags.ul(
                ui.tags.li("Pick 4–8 colonies and perform colony PCR to screen for insert presence."),
                ui.tags.li("Alternatively, miniprep and perform a diagnostic restriction digest."),
                ui.tags.li("Miniprep positive clones and submit for Sanger sequencing."),
                ui.tags.li("Verify correct insert sequence, orientation, and presence of attP site."),
            ),
            ui.br(),
        )
    ]

    return ui_card("Donor Plasmid Construction (Insert + attP)", 'integrase_protocol_section',
                   *integrase_section)
def generate_experimental_validation_section(guides_df, pam_type):
    # Extract the guide RNAs from the DataFrame
    guide_rnas = guides_df["Base Editing Guide"].tolist()

    # Process the guide RNAs
    processed_guide_rnas = process_guide_rnas(guide_rnas)

    # Create the validation section UI elements
    validation_section = [
        ui.help_text(
            ui.tags.span("1. Order the following paired oligos from "),
            ui.tags.a("IDT", href="https://www.idtdna.com", target="_blank"),
            ui.tags.span(" (or other preferred vendor):")
        ),
        ui.br(),
        ui.br(),
    ]

    # Add each processed guide RNA to the validation section
    for guide_rna in processed_guide_rnas:
        guide_rna_parts = guide_rna.split('\n')
        for part in guide_rna_parts:
            validation_section.append(ui.help_text(part))
            validation_section.append(ui.br())
        validation_section.append(ui.br())  # Add a double space between guide RNAs

    validation_section.extend([
        ui.help_text(
            ui.tags.span("2. Follow Cloning Protocol "),
            ui.tags.a("here.",
                      href="https://drive.google.com/file/d/1RBkiXhl7CCKxHuDNkaogvjkEykxkH4ws/view?usp=sharing",
                      target="_blank"),
        ),
    ])

    return ui_card("Experimental Validation of Base Editing Guide RNAs", 'validation_section', *validation_section)

#-------------------------------------------------------------------------------------
#VISUALIZATION

def generate_prime_del_visualization(guides_df, ref_sequence_input, substitution_position, PAM):

    prime_del_visualization_elements = []
    # Add the introductory text with colored formatting
    intro_text = ui.help_text(
        "This section visualizes PRIME-Del (paired prime editing). ",
        "The ", ui.tags.b("red", style="color: red;"), " and ",
        ui.tags.b("blue", style="color: blue;"),
        " segments simply distinguish the 3′ DNA flaps synthesized from each pegRNA’s RT template—the two ends that will be joined. ",
        "The dotted region marks the DNA segment programmed for deletion. ",
        "Steps shown: nick by PE2 (Cas9 H840A), reverse transcription to form 3′ flaps, and DNA repair to create a precise deletion (optionally with a short, encoded insert)."
    )



    # Center the image

    prime_del_visualization_elements.append(intro_text)
    prime_del_visualization_elements.append(
        ui.div(
            ui.output_image("prime_del_visualization_image"),
            {
                "style": "margin:0; padding:0; line-height:0;"  # kill card body line-height
            }
        )

    )
    prime_del_visualization_elements.append(ui.br())

    return prime_del_visualization_elements

def generate_prime_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM):
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()

    processed_pegRNA_oligos = visualization_peg_rnas(pegRNA_oligo_top,
                                                     pegRNA_oligo_extension_bottom) if pegRNA_oligo_top and any(
        pegRNA_oligo_top) else None
    processed_ngRNA_oligos = visualization_ng_rnas(ngRNA_oligos) if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[
        0] != 'n/a' else None

    prime_visualization_elements = []

    # Add the introductory text with colored formatting
    intro_text = ui.help_text(
        "This section provides a visualization of the recommended prime editing guide RNA design for installing a Bxb1 integration site. ",
        "The ", ui.tags.b("red", style="color: red;"), " characters represent the pegRNA spacer, ",
        "the ", ui.tags.b("blue", style="color: blue;"), " characters represent the SpCas9 scaffold, ",
        "the ", ui.tags.b("orange", style="color: orange;"), " characters represent the pegRNA extension sequence, ",
        "and the ", ui.tags.b("violet", style="color: violet;"),
        " characters represent the nicking guide RNA (ngRNA, if applicable). Note that the ngRNA is only applicable to certain edits. ",
    )



    prime_visualization_elements.append(intro_text)
    prime_visualization_elements.append(ui.br())
    prime_visualization_elements.append(ui.br())

    # Center the image
    prime_visualization_elements.append(
        ui.div(
            ui.output_image("prime_visualization_image"),
            {"style": "text-align: left;"}
        )
    )

    prime_visualization_elements.append(ui.br())

    if processed_pegRNA_oligos:
        for guide_index, guide_rna in enumerate(processed_pegRNA_oligos):
            guide_rna_parts = guide_rna.split()
            formatted_guide_rna = []

            constant_sequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
            pegRNA_top = str(
                [oligo[4:24] for oligo in pegRNA_oligo_top if oligo and oligo != 'n/a'][0]).upper().replace('-', '')
            pegRNA_extension_bottom = str(
                [oligo[4:] for oligo in pegRNA_oligo_extension_bottom if oligo and oligo != 'n/a'][0]).upper().replace(
                '-', '')

            for part in guide_rna_parts:
                if part == constant_sequence:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: blue; font-family: Courier,courier"))
                elif part == pegRNA_top:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: red; font-family: Courier,courier"))
                elif part == pegRNA_extension_bottom:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: orange; font-family: Courier,courier"))
                else:
                    formatted_guide_rna.append(ui.tags.b(part, style="font-family: Courier,courier"))

            complete_guide_text = ui.help_text(*formatted_guide_rna)

            prime_visualization_elements.append(ui.help_text(f"pegRNA:"))
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(complete_guide_text)
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(ui.br())

    if processed_ngRNA_oligos:
        for ngRNA_index, ng_rna in enumerate(processed_ngRNA_oligos):
            ng_rna_parts = ng_rna.split('\n')
            formatted_ng_rna = []

            for part in ng_rna_parts:
                sub_parts = part.split()
                for sub_part in sub_parts:
                    if sub_part in ["5'", "3'", "CACC", "AAAC", "-"]:
                        formatted_ng_rna.append(ui.tags.b(sub_part, style="font-family: Courier,courier"))
                    else:
                        formatted_ng_rna.append(
                            ui.tags.b(sub_part, style="color: violet; font-family: Courier,courier"))
                formatted_ng_rna.append(ui.br())

            complete_ng_rna_text = []
            for part in formatted_ng_rna:
                complete_ng_rna_text.append(part)

            prime_visualization_elements.append(ui.help_text(f"ngRNA:"))
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(ui.help_text(*complete_ng_rna_text))
            prime_visualization_elements.append(ui.br())
            # prime_visualization_elements.append(ui.br())

    return prime_visualization_elements

def generate_twin_prime_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM):
    pegRNA_oligo_top = guides_df["pegRNA Spacer Oligo Top"].tolist()
    pegRNA_oligo_extension_bottom = guides_df['pegRNA Extension Oligo Bottom'].tolist()
    ngRNA_oligos = guides_df["ngRNA Oligo Top"].tolist()
    print("pegRNA_oligo_top", pegRNA_oligo_top)
    print("pegRNA_oligo_extension_bottom", pegRNA_oligo_extension_bottom)


    processed_pegRNA_oligos = visualization_peg_rnas(pegRNA_oligo_top,
                                                     pegRNA_oligo_extension_bottom) if pegRNA_oligo_top and any(
        pegRNA_oligo_top) else None
    processed_ngRNA_oligos = visualization_ng_rnas(ngRNA_oligos) if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[
        0] != 'n/a' else None
    print("processed_pegRNA_oligos", processed_pegRNA_oligos)
    prime_visualization_elements = []
    # Add the introductory text with colored formatting
    intro_text = ui.help_text(
        "This section provides a visualization of the recommended prime editing guide RNA design. ",
        "The ", ui.tags.b("red", style="color: red;"), " characters represent the pegRNA spacer, ",
        "the ", ui.tags.b("blue", style="color: blue;"), " characters represent the SpCas9 scaffold, ",
        "the ", ui.tags.b("orange", style="color: orange;"), " characters represent the pegRNA extension sequence, ",
        "and the ", ui.tags.b("violet", style="color: violet;"),
        " characters represent the nicking guide RNA (ngRNA, if applicable). Note that the ngRNA is only applicable to certain edits. ",
    )

    prime_visualization_elements.append(intro_text)
    prime_visualization_elements.append(ui.br())
    prime_visualization_elements.append(ui.br())

    # Center the image
    prime_visualization_elements.append(
        ui.div(
            ui.output_image("prime_visualization_image"),
            {"style": "text-align: left;"}
        )
    )

    prime_visualization_elements.append(ui.br())

    if processed_pegRNA_oligos:
        for guide_index, guide_rna in enumerate(processed_pegRNA_oligos):
            guide_rna_parts = guide_rna.split()
            formatted_guide_rna = []

            constant_sequence = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
            pegRNA_top = str(
                [oligo[4:24] for oligo in pegRNA_oligo_top if oligo and oligo != 'n/a'][0]).upper().replace('-', '')
            pegRNA_extension_bottom = str(
                [oligo[4:] for oligo in pegRNA_oligo_extension_bottom if oligo and oligo != 'n/a'][0]).upper().replace(
                '-', '')

            for part in guide_rna_parts:
                if part == constant_sequence:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: blue; font-family: Courier,courier"))
                elif part == pegRNA_top:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: red; font-family: Courier,courier"))
                elif part == pegRNA_extension_bottom:
                    formatted_guide_rna.append(ui.tags.b(part, style="color: orange; font-family: Courier,courier"))
                else:
                    formatted_guide_rna.append(ui.tags.b(part, style="font-family: Courier,courier"))

            complete_guide_text = ui.help_text(*formatted_guide_rna)

            prime_visualization_elements.append(ui.help_text(f"pegRNA:"))
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(complete_guide_text)
            print("complete_guide_text", complete_guide_text)
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(ui.br())

    if processed_ngRNA_oligos:
        for ngRNA_index, ng_rna in enumerate(processed_ngRNA_oligos):
            ng_rna_parts = ng_rna.split('\n')
            formatted_ng_rna = []

            for part in ng_rna_parts:
                sub_parts = part.split()
                for sub_part in sub_parts:
                    if sub_part in ["5'", "3'", "CACC", "AAAC", "-"]:
                        formatted_ng_rna.append(ui.tags.b(sub_part, style="font-family: Courier,courier"))
                    else:
                        formatted_ng_rna.append(
                            ui.tags.b(sub_part, style="color: violet; font-family: Courier,courier"))
                formatted_ng_rna.append(ui.br())

            complete_ng_rna_text = []
            for part in formatted_ng_rna:
                complete_ng_rna_text.append(part)

            prime_visualization_elements.append(ui.help_text(f"ngRNA:"))
            prime_visualization_elements.append(ui.br())
            prime_visualization_elements.append(ui.help_text(*complete_ng_rna_text))
            prime_visualization_elements.append(ui.br())
            # prime_visualization_elements.append(ui.br())

    return prime_visualization_elements

def generate_integrase_visualization(guides_df, ref_sequence_input, substitution_position, PAM):

    integrase_visualization_elements = []
    # Add the introductory text with colored formatting
    intro_text = ui.help_text(
        "This section visualizes the BxB1 integrase-mediated DNA integration process. ",
        "The ", ui.tags.b("orange", style="color: orange;"), " characters represent the BxB1 attB site in the genome, ",
        "the ", ui.tags.b("blue", style="color: blue;"),
        " characters represent the BxB1 attP site on the donor plasmid, ",
        "and the ", ui.tags.b("red", style="color: red;"),
        " characters represent the inserted DNA sequence being integrated. ",
        "After integration, the hybrid sites ", ui.tags.b("attL/attR", style="color: black;"),
        " flank the inserted sequence. ",
        " This process enables precise, site-specific insertion of large DNA payloads into the genome."
    )


    # Center the image

    integrase_visualization_elements.append(intro_text)
    integrase_visualization_elements.append(ui.br())
    integrase_visualization_elements.append(ui.br())
    integrase_visualization_elements.append(
        ui.div(
            ui.output_image("integrase_visualization_image"),
            {"style": "text-align: left;"}
        )
    )
    integrase_visualization_elements.append(ui.br())

    return integrase_visualization_elements






def generate_prime_del_section(prime_editing_guides_df):
    """
    Generate a styled UI section for displaying Prime-Del pegRNAs and related info.
    """
    section_elements = []

    section_elements.append(ui.header("Prime-Del pegRNAs Identified"))
    section_elements.append(ui.text("These pegRNAs are optimized for deletions >80bp using Prime-Del technology."))

    if prime_editing_guides_df.empty:
        section_elements.append(ui.text("No pegRNAs available."))
        return ui.container(section_elements)

    # Iterate rows and display pegRNA info in consistent style
    for _, row in prime_editing_guides_df.iterrows():
        section_elements.append(
            ui.card([
                ui.text(f"**pegRNA ID:** {row.get('pegRNA ID', 'N/A')}"),
                ui.text(f"**Sequence:** {row.get('pegRNA Sequence', 'N/A')}"),
                ui.text(f"**PAM:** {row.get('PAM', 'N/A')}"),
                ui.text(f"**Predicted Efficiency:** {row.get('Efficiency', 'N/A')}"),
                ui.text(f"**Notes:** {row.get('Notes', '')}")
            ])
        )

    return ui.container(section_elements)
def parse_base_editing_window(window_str):
    start, end = map(int, window_str.split('-'))
    return start - 1, end  # Adjust for 0-based index


def generate_base_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM, window):
    start, end = parse_base_editing_window(window)
    list_of_guides_to_display = []
    for index, row in guides_df.iterrows():
        if row['Editing Technology'] == 'Base Editing':
            guide = row["Base Editing Guide"]
            orientation = row["Base Editing Guide Orientation"]

            ref_sequence_almost_rc = almost_reverse_complement(ref_sequence_input)
            base_to_edit = ref_sequence_input[substitution_position] if orientation != 'reverse' else \
            ref_sequence_almost_rc[substitution_position]
            if orientation == 'reverse':
                guide = guide[::-1]
                all_guide_occurance_starts = [m.start() for m in re.finditer(guide, ref_sequence_almost_rc)]
            else:
                all_guide_occurance_starts = [m.start() for m in re.finditer(guide, ref_sequence_input)]

            true_starting_positions = []
            for start_pos in all_guide_occurance_starts:
                end_pos = start_pos + len(guide) - 1
                if substitution_position >= start_pos and substitution_position <= end_pos:
                    true_starting_positions.append(start_pos)
            assert len(true_starting_positions) == 1, (
                "Error! Guide cannot be aligned properly to input original sequence", guide, ref_sequence_almost_rc,
                orientation, all_guide_occurance_starts, substitution_position)
            guide_start = true_starting_positions[0]

            list_of_guides_to_display.append(ui.help_text(f"Guide {index + 1}"))
            list_of_guides_to_display.append(ui.br())
            if orientation == 'reverse':
                reverse_guide_rna_with_formatting = []
                for i in range(len(guide)):
                    base_position = guide_start + len(guide) - 1 - i
                    base = ref_sequence_almost_rc[base_position]
                    if i >= start and i < end:  # Bases within the selected window
                        if base_position == substitution_position:  # Base to edit
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base,
                                                                               style="text-decoration: underline; color: green; font-family: Courier,courier"))
                        elif base == base_to_edit:
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base,
                                                                               style="text-decoration: underline; color: red; font-family: Courier,courier"))
                        else:
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base,
                                                                               style="text-decoration: underline; color: violet; font-family: Courier,courier"))
                    else:
                        reverse_guide_rna_with_formatting.append(
                            ui.tags.b(base, style="color: violet; font-family: Courier,courier"))

                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Forward Strand: 5'-" + ref_sequence_input[:substitution_position],
                                  style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_input[substitution_position],
                                  style="color: green; font-family: Courier,courier"),  # base to edit
                        ui.tags.b(ref_sequence_input[substitution_position + 1:] + "-3'",
                                  style="font-family: Courier,courier"),
                    )
                )
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Reverse Strand: 3'-" + ref_sequence_almost_rc[:guide_start - len(PAM)],
                                  style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_almost_rc[guide_start - len(PAM):guide_start],
                                  style="color: blue; font-family: Courier,courier"),  # PAM
                        *reverse_guide_rna_with_formatting[::-1],  # reverse guide RNA
                        ui.tags.b(ref_sequence_almost_rc[guide_start + len(guide):] + "-5'",
                                  style="font-family: Courier,courier"),
                    )
                )
            else:
                guide_rna_with_formatting = []
                for i in range(len(guide)):
                    base = ref_sequence_input[guide_start + i]
                    if i >= start and i < end:  # Bases within the selected window
                        if guide_start + i == substitution_position:  # Base to edit
                            guide_rna_with_formatting.append(ui.tags.b(base,
                                                                       style="text-decoration: underline; color: green; font-family: Courier,courier"))
                        elif base == base_to_edit:
                            guide_rna_with_formatting.append(ui.tags.b(base,
                                                                       style="text-decoration: underline; color: red; font-family: Courier,courier"))
                        else:
                            guide_rna_with_formatting.append(ui.tags.b(base,
                                                                       style="text-decoration: underline; color: violet; font-family: Courier,courier"))
                    else:
                        guide_rna_with_formatting.append(
                            ui.tags.b(base, style="color: violet; font-family: Courier,courier"))

                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Forward Strand: 5'-" + ref_sequence_input[:guide_start],
                                  style="font-family: Courier,courier"),
                        *guide_rna_with_formatting,
                        ui.tags.b(ref_sequence_input[len(guide) + guide_start:len(guide) + guide_start + len(PAM)],
                                  style="color: blue; font-family: Courier,courier"),  # PAM
                        ui.tags.b(ref_sequence_input[guide_start + len(guide) + len(PAM):] + "-3'",
                                  style="font-family: Courier,courier"),
                    )
                )
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Reverse Strand: 3'-" + ref_sequence_almost_rc[:substitution_position],
                                  style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_almost_rc[substitution_position],
                                  style="color: green; font-family: Courier,courier"),  # base to edit
                        ui.tags.b(ref_sequence_almost_rc[substitution_position + 1:] + "-5'",
                                  style="font-family: Courier,courier"),
                    )
                )
            if index != guides_df.shape[0] - 1:
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(ui.br())

    return list_of_guides_to_display


def server(input, output, session):
    # Define default values
    DEFAULT_BASE_EDITING_WINDOW_START = 4
    DEFAULT_BASE_EDITING_WINDOW_END = 9

    # DEFAULT_RS3_TOGGLE = False

    # Function to get base editing window values
    def get_base_editing_window(mutation_type=None):
        if mutation_type == "CGB":
            return DEFAULT_BASE_EDITING_WINDOW_START, DEFAULT_BASE_EDITING_WINDOW_END  # Override for CGB
        if input.advanced_settings_toggle() % 2 == 1:
            return input.base_editing_window_start(), input.base_editing_window_end()
        else:
            return DEFAULT_BASE_EDITING_WINDOW_START, DEFAULT_BASE_EDITING_WINDOW_END

    # def get_calculate_rs3():
    #     if input.advanced_settings_toggle() % 2 == 1:
    #         return input.rs3_toggle()
    #     else:
    #         return DEFAULT_RS3_TOGGLE

    @output
    @render.image
    def prime_visualization_image():
        img: ImgData = {"src": str(Path(__file__).parent / "prime_editing_diagram.png"), "width": "500px"}
        return img

    @output
    @render.image
    def prime_del_visualization_image():
        return {
            "src": str(Path(__file__).parent / "Prime_del_diagram.png"),
            "style": (
                "display:block;"  # no inline baseline gap
                "margin:0;"
                "padding:0;"
                "border:0;"
                "max-width:100%;"
                "height:auto;"
                "vertical-align:top;"  # <-- critical to kill the white gap under the image
            )
        }
        return img

    @output
    @render.image
    def integrase_visualization_image():
        img: ImgData = {"src": str(Path(__file__).parent / "integrase_diagram.png"),"width": "1000px"}
        return img

    def input_check(ref_sequence_input, edited_sequence_input, df=None, edit_start=4, edit_end=9):
        nonlocal input_file
        if input_file and not (ref_sequence_input or edited_sequence_input):
            try:
                df = pd.read_csv(input_file)
            except Exception as e:
                return False, "Input file is not a properly formed CSV file. Please input a proper CSV file."

            if len(df.columns) != 2 or df.columns.tolist() != ['Original Sequence', 'Desired Sequence']:
                return False, 'Uploaded csv does not have the proper columns. Your csv must have two columns with names "Original Sequence" and "Desired Sequence"'

            counter = 1
            for index, row in df.iterrows():
                ref_sequence = "".join(row['Original Sequence'].split()).upper()
                edited_sequence = "".join(row['Desired Sequence'].split()).upper()
                check, message = check_ref_edited_pair(ref_sequence, edited_sequence, edit_start, edit_end)
                if not check:
                    return check, f"Error row {counter}: {message}"
                counter += 1
            return True, "Input CSV verified. Proceed to get guides."
        elif ref_sequence_input and edited_sequence_input and not input_file:
            check, message = check_ref_edited_pair("".join(ref_sequence_input.split()).upper(),
                                                   "".join(edited_sequence_input.split()).upper(), edit_start, edit_end)
            return check, message
        elif ref_sequence_input and edited_sequence_input and input_file:
            return False, "Error: Fill in both text input fields or upload a CSV file but do not do both."
        else:
            return False, "Error: Fill in both text input fields or upload a CSV file."

    input_file = None
    user_selected_pam = reactive.Value("")

    @output
    @render.ui
    @reactive.event(input.upload)
    def upload_status():
        nonlocal input_file
        file_infos = input.file1()
        if file_infos:
            input_file = file_infos[0]['datapath']
            return ui.div(ui.br(),
                          ui.tags.b("File Successfully Uploaded", style="color: grey;", id='upload_status_message'))
        else:
            return ui.div(ui.br(),
                          ui.tags.b("Error: No file selected", style="color: red;", id='upload_status_message'))

    @reactive.Effect()
    def clear():
        value = input.clear()
        if value > 0:
            ui.update_text_area("ref_sequence_input", value="")
            ui.update_text_area("edited_sequence_input", value="")
            ui.update_select("pam_type", selected='NGN')
            ui.remove_ui(selector="div:has(> #results)")
            ui.remove_ui(selector="div:has(> #upload_status_message)")
            ui.remove_ui(selector="#error_message")

            # Clear the file input UI
            ui.remove_ui(selector="div:has(> #file1)")

            # Re-render the file input UI
            @output
            @render.ui
            def ui_input_file():
                return ui.input_file("file1", 'Choose a CSV File of Sequences to Upload:', accept='.csv',
                                     multiple=False, width="100%")

            nonlocal input_file
            input_file = None

    @output
    @render.ui
    def ui_input_file():
        return ui.input_file(f"file1", 'Choose a CSV File of Sequences to Upload:', accept='.csv', multiple=False,
                             width="100%")

    @reactive.Effect
    @reactive.event(input.file1)
    def process_file_upload():
        nonlocal input_file
        file_infos = input.file1()
        if file_infos:
            input_file = file_infos[0]['datapath']
            # Call the function to process the file
            process_uploaded_file(input_file)

    def process_uploaded_file(file_path):
        try:
            # read the file
            df = pd.read_csv(file_path)
        except Exception as e:
            print(f"Error reading CSV: {e}")

    @output
    @render.ui
    def advanced_settings():
        if input.advanced_settings_toggle() % 2 == 1:
            return ui.TagList(
                ui.input_select("pam_type", "Select Desired Base Editing PAM",
                                {"NGG": "NGG (Most Efficient)", "NGN": "NGN (Recommended)", "NGA": "NGA",
                                 "NNGRRT": "NNGRRT (SaCas9)", "NNNRRT": "NNNRRT (SaCas9-KKH)", "NRN": "NRN (SpRY)"}),
                ui.div(
                    {"style": "display: flex; gap: 10px;"},
                    ui.input_numeric("base_editing_window_start", "Base Editing Window Start", value=4, min=1, max=20),
                    ui.input_numeric("base_editing_window_end", "Base Editing Window End", value=9, min=1, max=20)
                ),
                #     ui.input_checkbox("rs3_toggle", "Calculate RuleSet3 Scores?", value=True)
            )
        else:
            return ui.div()

    @reactive.Effect()
    def handle_pam_selection():
        if input.advanced_settings_toggle() % 2 == 1:
            user_selected_pam.set(input.pam_type())
        else:
            user_selected_pam.set("")

    def display_error_message(message):
        return ui.div(ui.tags.b(message, style="color: red;", id='error_message'))

    # function to run if the user inputs a CSV
    @output
    @render.ui
    @reactive.event(input.find_guides_csv)
    def run_with_csv_input():
        nonlocal input_file
        if not input_file:
            return display_error_message("No CSV file uploaded")

        df = pd.read_csv(input_file)

        # Get base editing window values
        base_editing_window_start, base_editing_window_end = get_base_editing_window()
        base_editing_window = f"{base_editing_window_start}-{base_editing_window_end}"

        valid_inputs, message = input_check(None, None, df, base_editing_window_start, base_editing_window_end)
        if not valid_inputs:
            return display_error_message(message)

        selected_PAM = user_selected_pam.get() or "NGG"  # Default to NGG if no PAM is selected
        PAM = selected_PAM

        @output
        @render.data_frame
        def render_Results():
            nonlocal to_display_guides_df
            return render.DataGrid(
                to_display_guides_df,
                width="100%",
                filters=False,
                summary=False,
            )

        @session.download(
            filename=lambda: f"guides-{date.today().isoformat()}-{datetime.now().strftime('%H-%M-%S')}.csv")
        async def download_Results():
            nonlocal guides_df
            # Ensure the "Input CSV Row Number" is the first column
            cols = ['Input CSV Row Number'] + [col for col in guides_df.columns if col != 'Input CSV Row Number']
            guides_df = guides_df[cols]
            yield guides_df.to_csv(index=False)

        dfs_to_merge_download = []
        dfs_to_merge_display = []
        counter = 1

        with ui.Progress(min=1, max=df.shape[0] + 1) as p:
            p.set(message="Finding guides", detail="This may take a while...")
            for index, row in df.iterrows():
                p.set(counter, message="Finding guides")
                ref_sequence_input = "".join(row['Original Sequence'].split()).upper()
                edited_sequence_input = "".join(row['Desired Sequence'].split()).upper()
                to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input, selected_PAM,
                                                             base_editing_window_start, base_editing_window_end)
                index_column = [str(counter)] * to_display_guides_df.shape[0]
                to_display_guides_df.insert(loc=0, column='Input CSV Row Number', value=index_column)
                guides_df.insert(loc=0, column='Input CSV Row Number', value=index_column)
                dfs_to_merge_download.append(guides_df)
                dfs_to_merge_display.append(to_display_guides_df)
                counter += 1

        to_display_guides_df = pd.concat(dfs_to_merge_display)

        cols_to_drop = [c for c in ['Original Sequence', 'Desired Sequence'] if c in to_display_guides_df.columns]
        if cols_to_drop:
            to_display_guides_df = to_display_guides_df.drop(columns=cols_to_drop)

        guides_df = pd.concat(dfs_to_merge_download)

        ui_elements = [
            ui.help_text(
                '''Note: for base editing, we will show multiple guide RNAs, but for prime editing, we will only show the recommended PrimeDesign guide RNA. The off-target score is calculated using the CFD score algorithm (Doench et al. 2014) where a higher score indicates a lower likelihood of off-target activity. The on-target score is calculated using the RuleSet1 algorithm where a higher score indicates greater efficiency of guide RNA binding to the genomic target sequence. Both algorithms are from ''',
                ui.tags.a('Doench et al. 2014 Nat Biotechnol.',
                          {'href': 'https://pubmed.ncbi.nlm.nih.gov/25184501/', 'target': '_blank'}),
                ui.tags.b(" A higher score is better for both algorithms.", style="text-decoration: bold")
            ),
            ui.br(),
            ui.br(),
            ui.output_data_frame("render_Results"),
            ui.br(),
            ui.br(),
        ]

        ui_elements.append(ui.download_button("download_Results", "Download Results as CSV File"))

        return ui.TagList(
            ui_card(
                "Batch Sequence Results",
                'results',
                *ui_elements
            )
        )

    # function to run if the user inputs a text sequence
    @output
    @render.ui
    @reactive.event(input.find_guides_text)
    def run_with_text_input():
        ref_sequence_input = "".join(input.ref_sequence_input().split()).upper()
        #print(ref_sequence_input)

        edited_sequence_input = "".join(input.edited_sequence_input().split()).upper()
        #print(edited_sequence_input)

        # Get base editing window values
        base_editing_window_start, base_editing_window_end = get_base_editing_window()
        input_window = f"{base_editing_window_start}-{base_editing_window_end}"
        # print(f'input_window: {input_window}')

        selected_PAM = user_selected_pam.get() or "NGG"  # Default to NGG if no PAM is selected
        PAM = selected_PAM

        # calculate_rs3 = get_calculate_rs3()

        valid_inputs, message = check_ref_edited_pair(ref_sequence_input, edited_sequence_input,
                                                      base_editing_window_start, base_editing_window_end)
        if not valid_inputs:
            return display_error_message(message)

        to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input, selected_PAM,
                                                     base_editing_window_start, base_editing_window_end)
        to_display_guides_df = to_display_guides_df.drop(columns=['Original Sequence', 'Desired Sequence'])
        to_display_guides_df.insert(loc=0, column='Guide',
                                    value=[f"Guide {i + 1}" for i in range(to_display_guides_df.shape[0])])

        base_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Base Editing']
        prime_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Prime Editing']
        twin_prime_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Prime Editing (Creating a Bxb1 Site)']
        filtered_guides_df = to_display_guides_df[
            to_display_guides_df['Editing Technology'].isin(['Base Editing', 'Prime Editing'])]

        # Tells users to try other PAMs if default isn't available
        if not base_editing_guides_df.empty:
            pams = guides_df["PAM"].tolist()
            PAM = pams[0]
        else:
            PAM = selected_PAM

        # Define variables to find enzyme protein on Addgene in suggested Addgene plasmid section
        editor_information, mutation_type = get_editor_info(ref_sequence_input, edited_sequence_input, PAM)
        editor_name, editor_id, editor_url = editor_information
        editor_info = f"{editor_name} (Addgene: {editor_id})"

        base_editing_window_start, base_editing_window_end = get_base_editing_window(
            mutation_type)  # in case its a cgb edit
        base_editing_window = f"{base_editing_window_start}-{base_editing_window_end}"
        # print(f'actual window: {base_editing_window}')

        # Define variables to find guide RNA cloning plasmid on Addgene in suggested Addgene plasmid section
        cloning_name, cloning_id, cloning_url = get_cloning_url(PAM)
        plasmid_info = f"{cloning_name} (Addgene: {cloning_id})"

        # If no prime or base editing guides found, output a message
        message = check_for_BE_or_PE_guides(guides_df)
        if message:
            return ui.div(ui.tags.b(f"{message}", style="color: red;"),
                          ui.tags.a("https://www.nature.com/articles/nprot.2013.143",
                                    href="https://www.nature.com/articles/nprot.2013.143", target="_blank",
                                    style="color: red;"))

        @output
        @render.data_frame
        def render_results():
            nonlocal to_display_guides_df
            return render.DataGrid(
                to_display_guides_df,
                width="100%",
                filters=False,
                summary=False,
            )

        @render.download(
            filename=lambda: f"guides-{date.today().isoformat()}-{datetime.now().strftime('%H-%M-%S')}.csv"
        )
        def download_results():
            nonlocal guides_df
            yield guides_df.to_csv(index=False)

        substitution_position = None
        for i in range(len(ref_sequence_input)):
            ref_base = ref_sequence_input[i]
            edited_base = edited_sequence_input[i]
            if ref_base != edited_base:
                substitution_position = i
                break

        if len(ref_sequence_input) > 51:
            ref_sequence_input = ref_sequence_input[substitution_position - 25:substitution_position + 25 + 1]
            substitution_position = 25

        list_of_guides_to_display = generate_base_editing_visualization(guides_df, ref_sequence_input,
                                                                        substitution_position, PAM, base_editing_window)

        ui_elements = []
        if get_base_editing_window():
            message = check_CGB_window(base_editing_window, input_window)
            if message:
                ui_elements.append(ui.help_text(ui.tags.b(f"{message}", style="color: red;")))
                ui_elements.append(ui.br(), )
        if user_selected_pam.get():
            message = check_availible_pam(PAM, selected_PAM)
            if message:
                ui_elements.append(ui.help_text(ui.tags.b(f"{message}", style="color: red;")))
                ui_elements.append(ui.br(), )

        # Determine the guide title and help text based on the editing technology
        if "Base Editing" in filtered_guides_df['Editing Technology'].values:
            guide_title = "Recommended Base Editing Guide RNAs"
            help_text = ui.help_text(
                '''Note: Multiple guide RNAs are shown with the ability to toggle based on optimal on-target and off-target scoring. The off-target score is calculated using the CFD score algorithm where a higher score indicates a lower likelihood of off-target activity. The on-target score is calculated using the RuleSet1 algorithm where a higher score indicates greater efficiency of guide RNA binding to the genomic target sequence. For additional information on these algorithms, please see the original publication ''',
                ui.tags.a('Doench et al. 2014 Nat Biotechnol.',
                          {'href': 'https://pubmed.ncbi.nlm.nih.gov/25184501/', 'target': '_blank'}),
                ui.tags.b(" A higher score is better for both algorithms.", style="text-decoration: bold")
            )
        elif any(col.startswith("Prime-Del ") for col in guides_df.columns):
            guide_title = "Recommended PRIME-Del pegRNA sequences"
            help_text = ui.help_text("Note: PRIME-Del is a genome editing method based on prime "
                                     "editing that facilitates the precise deletion of DNA sequences. "       
                                    "We use the Shendure Labs's algorithm with default parameters, "
                                    "to identify the single most optimal prime editing guide RNAs. For more advanced usage "
                                    "please visit the ",
                                    ui.tags.a(
                                        "Shendure Labs's PRIME-Del portal",
                                        {"href": "https://shendurelab.github.io/Prime-del/", "target": "_blank"}
                                    ),
                                    " or see the original publication ",
                                    ui.tags.a(
                                        "(Choi, J., et al. 2021 Nat Biotech)",
                                        {"href": "https://www.nature.com/articles/s41587-021-01025-z", "target": "_blank"}
                                    ),
                                    ui.br(),
                                     "."
                                     )
            ui_elements.append(
                ui_card(
                    guide_title,
                    'guides_df',
                    help_text,
                    ui.br(),
                    ui.output_data_frame("render_results")
                )
            )


        else:
            guide_title = "Recommended Prime Editing Guide RNAs"
            help_text = ui.help_text(
                "Note: We use the PrimeDesign algorithm with default parameters (including NGG PAM), "
                "to identify the single most optimal prime editing guide RNA. For more advanced usage "
                "please visit the ",
                ui.tags.a(
                    "PrimeDesign portal",
                    {"href": "https://prime-design-766275861038.us-east1.run.app/", "target": "_blank"}
                ),
                " or see the original publication ",
                ui.tags.a(
                    "(Hsu et al. 2021 Nat Commun)",
                    {"href": "https://pubmed.ncbi.nlm.nih.gov/33589617/", "target": "_blank"}
                ),
                "."
            )

            explore_button = ui.tags.button(
                "Explore on PrimeDesign",
                id="toggle_prime",
                style="margin:10px; padding:6px 12px; border-radius:4px;"
            )

            prime_iframe = ui.div(
                ui.tags.iframe(
                    src="https://prime-design-766275861038.us-east1.run.app/",
                    style="width:100%; height:80vh; border:1px solid #ccc; border-radius:8px; display:none;",
                    id="prime_iframe",
                    title="PrimeDesign Portal",
                    loading="lazy"
                )
            )

            # One event listener handles open/close + table hide/show
            toggle_script = ui.tags.script("""
            (function () {
              function $(id){ return document.getElementById(id); }
              document.addEventListener("click", function(e) {
                if (!e.target) return;
                if (e.target.id !== "toggle_prime" && e.target.id !== "toggle_twin_prime") return;

                var isTwin = e.target.id === "toggle_twin_prime";
                var iframe = $(isTwin ? "twin_prime_iframe" : "prime_iframe");
                var tableCard = $("guides_table_card");

                var opening = (iframe.style.display === "" || iframe.style.display === "none");
                iframe.style.display = opening ? "block" : "none";
                e.target.innerText = opening ? "Close PrimeDesign" : "Explore on PrimeDesign";
                if (tableCard) tableCard.style.display = opening ? "none" : "";
              }, true);
            })();
            """)

            ui_elements.append(help_text)
            ui_elements.append(explore_button)
            ui_elements.append(prime_iframe)
            ui_elements.append(toggle_script)

            ui_elements.append(
                ui.div(
                    ui_card(
                        guide_title,
                        "guides_df",
                        ui.output_data_frame("render_results")
                    ),
                    {"id": "guides_table_card"}
                )
            )

        if 'Base Editing' in filtered_guides_df['Editing Technology'].values:
            ui_elements.append(
                ui_card(
                    "Suggested Addgene Plasmids",
                    "recommended_editor",
                    ui.help_text(f"Recommended Guide RNA Plasmid: {plasmid_info}"),
                    ui.br(),
                    ui.help_text(f"Recommended Base Editor Plasmid: {editor_info}"),
                    ui.br(),
                    ui.br(),
                    ui.div(
                        ui.tags.a("View Guide RNA Plasmid", href=cloning_url, target="_blank",
                                  class_="btn btn-primary me-2"),
                        ui.tags.a("View Base Editor Plasmid", href=editor_url, target="_blank",
                                  class_="btn btn-primary"),
                        style="display: flex; gap: 5px; justify-content: left; margin-bottom: 30px;"
                    ),
                )
            )
            ui_elements.append(
                ui_card(
                    "Visualization of Base Editing Guides",
                    "base_editing_visualization",
                    ui.help_text(
                        "For each base editing guide, your input will be displayed with the guide sequence highlighted on the appropriate strand.",
                        ui.tags.b(" Green", style="color: green"),
                        " characters represent your edited base, ",
                        ui.tags.b("blue", style="color: blue"),
                        " characters represent the PAM nucleotides, ",
                        ui.tags.b("violet", style="color: violet"),
                        " characters represent other nucleotides in the guide, ",
                        ui.tags.b("underlined", style="text-decoration: underline"),
                        " characters represent characters the base editing window, and ",
                        ui.tags.b("red", style="color: red"),
                        " characters represent any bystander edits within the editing window. Grey characters represent nucleotides not spanned by the guide.  NOTE: Only 25 bp of sequence upstream and downstream of the desired edit is shown. Also, toggling the order of the guides, the guide number above corresponds to the guide number within the visualization."
                    ),
                    ui.br(),
                    ui.br(),
                    *list_of_guides_to_display
                )
            )
            ui_elements.append(generate_experimental_validation_section(base_editing_guides_df, PAM))

        if 'Prime Editing' in filtered_guides_df['Editing Technology'].values:
            prime_editing_plasmid_card = create_prime_editing_plasmid_card(guides_df, editor_info, editor_url)

            ui_elements.append(prime_editing_plasmid_card)

            ui_elements.append(ui_card(
                "Visualization of Prime Editing Guides",
                "prime_editing_visualization",
                *generate_prime_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM)
            ))

            ui_elements.append(generate_prime_protocols_section(prime_editing_guides_df))

        if 'Prime Editing (Creating a Bxb1 Site)' in guides_df['Editing Technology'].values:
            twin_prime_editing_plasmid_card = create_twin_prime_editing_plasmid_card(twin_prime_editing_guides_df,
                                                                                     editor_info, editor_url)
            ui_elements.append(ui_card(
                "Visualization of Prime Editing Guides - Creating Bxb1 Integration site",
                "prime_editing_visualization",
                *generate_twin_prime_editing_visualization(twin_prime_editing_guides_df, ref_sequence_input,
                                                           substitution_position, PAM)
            ))
            ui_elements.append(generate_twin_prime_protocols_section(twin_prime_editing_guides_df))


            integrase_plasmid_card = create_integrase_plasmid_card(guides_df, editor_info, editor_url)
            ui_elements.append(integrase_plasmid_card)

            ui_elements.append(generate_donor_protocols_section(twin_prime_editing_guides_df))

            ui_elements.append(ui_card(
                "Visualization of Bxb1 Integrase Step",
                "prime_editing_visualization",
                *generate_integrase_visualization(guides_df, ref_sequence_input, substitution_position, PAM)
            ))
    
            ui_elements.append(generate_integrase_protocols_section(twin_prime_editing_guides_df))

        if 'Deletion >80bp: Use Prime-Del' in guides_df['Editing Technology'].values or any(col.startswith('Prime-Del ') for col in guides_df.columns):
            # This is where I want to output the pegrnas from prime-del and any other good info from the output, it should be formatted well and in the same style as the rest o the application
            #This is just a placeholder
            prime_del_plasmid_card = create_prime_del_plasmid_card(guides_df, editor_info, editor_url)
            ui_elements.append(prime_del_plasmid_card)
            ui_elements.append(ui_card(
                "Visualization of PRIME-Del",
                "prime_editing_visualization",
                *generate_prime_del_visualization(guides_df, ref_sequence_input, substitution_position, PAM)
            ))


            ui_elements.append(generate_prime_del_protocols_section(twin_prime_editing_guides_df))


        ui_elements.append(ui.download_button("download_results", "Download Results as CSV File"))
        return ui.TagList(
            ui_card(
                "Single Sequence Results",
                'results',
                *ui_elements
            )
        )

    @output
    @render.image
    def display_logo():
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "EditABLE-logos_transparent.png"), "width": "300px"}
        return img

    @output
    @render.image
    def stanford_logo():
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "SOM_Web_vert_LG.png"), "width": "300px"}
        return img


app = App(app_ui, server)
