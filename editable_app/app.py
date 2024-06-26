import asyncio
import os
import pandas as pd
import re
from datetime import date, datetime
from pathlib import Path
from utils import get_guides, almost_reverse_complement, get_editor_info, process_guide_rnas, get_cloning_url, process_ng_rnas,  process_peg_rnas
from shiny import App, render, ui, reactive
from shiny.types import ImgData

bases = {"A", "C", "G", "T"}
accepted_bases = {"A", "C", "G", "T", "-"}


# A card component wrapper.
def ui_card(title, id, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body", "id" : id}, *args),
        ),
    )

app_ui = ui.page_fluid(
    {"id": "main-content"},
    ui.output_image("display_logo", inline=True),
    ui.output_image("stanford_logo", inline=True),
    ui.br(),
    ui.help_text(
        '''Welcome to EditABLE! We have designed this tool as a starting point for determining the optimal gene editing tool for a single gene edit. 
        Our algorithm is pending publication, but at a high level, it prioritizes base editing due to greater efficiency followed by prime editing. 
        Please refer to the following papers for more information on base and prime editing:'''
    ),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a('Komor AC, Kim YB, Packer MS, Zuris JA, Liu DR. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. Nature. 2016;533(7603):420-424.', {'href' : 'https://pubmed.ncbi.nlm.nih.gov/27096365/', 'target' : '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a('Gaudelli NM, Komor AC, Rees HA, et al. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage. Nature. 2017;551(7681):464-471.', {'href' : 'https://pubmed.ncbi.nlm.nih.gov/29160308/', 'target' : '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a('Anzalone AV, Randolph PB, Davis JR, et al. Search-and-replace genome editing without double-strand breaks or donor DNA. Nature. 2019;576(7785):149-157.', {'href' : 'https://pubmed.ncbi.nlm.nih.gov/31634902/', 'target' : '_blank'})),
    ui.br(),
    ui.br(),
    ui.help_text(ui.tags.a('Doman JL, Pandey S, Neugebauer ME, et al. Phage-assisted evolution and protein engineering yield compact, efficient prime editors. Cell. 2023;186(18):3983-4002.e26.', {'href' : 'https://pubmed.ncbi.nlm.nih.gov/37657419/', 'target' : '_blank'})),
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
            ui.tags.b("T", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For changes that result in deletions, use a string of "-" characters in the desired sequence to denote the 
            area of the deletion. For example, this would be a valid set of inputs for a deletion:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Original Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("ATT", style="color: red; font-family: Courier,courier"), 
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
        ui.tags.span("Please limit the length of the insertion or deletion to less than 40 nucleotides as prime editors have been primarily tested in this context. Refer to  "),
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
        ui.input_text_area("ref_sequence_input", "Original Sequence", placeholder="Enter sequence", height="50%", width="100%"),
        ui.input_text_area("edited_sequence_input", "Desired Sequence", placeholder="Enter sequence", height="50%", width="100%"),
        ui.output_ui("ui_input_file"),
        ui.input_action_button("upload", "Upload File", class_="btn-primary"),
        ui.output_ui("upload_input"),
        ui.br(),
        ui.input_select("pam_type", "Select Desired Base Editing PAM", {"NGN": "NGN (Recommended)", "NGG": "NGG (Most Efficient)", "NGA" : "NGA", "NNGRRT" : "NNGRRT (SaCas9)", "NNNRRT" : "NNNRRT (SaCas9-KKH)", "NRN" : "NRN (SpRY)"}),
        ui.input_action_button("get_guides", "Find Guides", class_="btn-primary"),
        ui.help_text(" "),
        ui.input_action_button("clear", "Clear Inputs", class_="btn-danger"),
    ),
    ui.output_ui("run"),
    ui.br(),
    ui.help_text(
        '''For base editing, EditABLE finds guide RNAs where the editable base is in positions 4-9 starting from the 5' 
           end of the guide RNA. Note that guide RNAs with base edits outside of this editing range may still work but 
           tend to be significantly less efficient.'''
    ),
    ui.br(),
    ui.br(),
    ui.help_text(
        '''Note that we use PrimeDesign (Hsu, J.Y., Grünewald, J., Szalay, R. et al. 
        PrimeDesign software for rapid and simplified design of prime editing guide RNAs. 
        Nat Commun 12, 1034 (2021)) for prime editing guide calcualtions. We run PrimeDesign 
        with default parameters and take only the suggested guides. For more advanced usage 
        please use the ''',
        ui.tags.a('PrimeDesign portal', {'href' : 'https://primedesign.pinellolab.partners.org/', 'target' : '_blank'}),
        ''' (https://primedesign.pinellolab.partners.org/) to design your Prime Editing guides.'''
    ),
    ui.br(),
    ui.br(),
    ui.help_text(
        '''For troubleshooting and suggested revisions, please contact the ''',
        ui.tags.a("Bhalla Lab", {'href' : 'https://med.stanford.edu/bhallalab.html', 'target' : '_blank'}),
        " at vbhalla@stanford.edu"
    ),
    ui.br(),
    ui.br(),
    ui.br(),
    ui.br(),
)

#checks if no guides found
def check_for_BE_or_PE_guides(guides_df):
    no_BE_or_PE = guides_df[guides_df['Editing Technology'].isin(['No Base or Prime Editing Guides Found'])]
    if not no_BE_or_PE.empty:
        return "No Prime Editing or Base Editing Guide RNAs available for this site. In this case, we recommend the use of CRISPR with homologous recombination. See this paper here: "


#checks if pam defaults to nrn or nyn bc no guide rnas for selected pam
def check_availible_pam(PAM, selected_PAM):
    if (PAM == "NRN") or (PAM == "NYN"): 
        if (selected_PAM != "NRN"):
            return f"No Guide RNAs available with {selected_PAM} PAM. See available NRN PAMs below:"

def check_ref_edited_pair(ref_sequence, edited_sequence):
    if len(ref_sequence) == 0 or len(edited_sequence) == 0:
        return False, "Both the original sequence and the edited sequence must be of nonzero length."
    if len(ref_sequence) != len(edited_sequence):
        return False, "The length of the original sequence and the edited sequence must be the same."
    if ref_sequence == edited_sequence:
        return False, "The original sequence and the edited sequence are the same."
    if len(set(ref_sequence) - accepted_bases) != 0 or len(set(edited_sequence) - accepted_bases) != 0:
        return False, "You may only have the following characters in your sequences {A, C, G, T, a, c, g, t, -}."
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

def create_prime_editing_plasmid_card(guides_df, editor_info, editor_url):
    pegRNA_oligos = guides_df["PrimeDesign pegRNA Spacer Oligo Top"].tolist()
    ngRNA_oligos = guides_df["PrimeDesign ngRNA Oligo Top"].tolist()
    
    processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligos) if pegRNA_oligos and any(pegRNA_oligos) else None
    processed_ngRNA_oligos = process_ng_rnas(ngRNA_oligos) if ngRNA_oligos and any(ngRNA_oligos) and ngRNA_oligos[0] != 'n/a' else None
    
    prime_card_elements = [
        ui.help_text(f"Recommended Prime Editing Plasmid: {editor_info}"),
        ui.br(),

    ]
    
    # Add the recommended pegRNA plasmid information if pegRNA oligos are not empty
    if processed_pegRNA_oligos:
        prime_card_elements.append(ui.help_text("Recommended pegRNA Plasmid: pU6-tevopreq1-GG-acceptor (Addgene: 174038)"))
        prime_card_elements.append(ui.br())
    
    # Add the recommended ngRNA plasmid information if ngRNA oligos are not empty
    if processed_ngRNA_oligos:
        prime_card_elements.append(ui.help_text("Recommended ngRNA Plasmid: pmCherry-U6-empty (Addgene: 140580)"))
        prime_card_elements.append(ui.br())
    
    prime_card_elements.append(ui.br())
    
    # Add buttons at the bottom next to each other
    prime_card_elements.append(
        ui.div(
            {"class": "d-flex"},
            ui.tags.a("View Prime Editing Plasmid", href=editor_url, target="_blank", class_="btn btn-primary me-2"),
            ui.tags.a("View pegRNA Plasmid", href="https://www.addgene.org/174038/", target="_blank", class_="btn btn-primary me-2") if processed_pegRNA_oligos else None,
            ui.tags.a("View ngRNA Plasmid", href="https://www.addgene.org/65777/", target="_blank", class_="btn btn-primary") if processed_ngRNA_oligos else None,
        )
    )
    
    return ui_card(
        "Suggested Addgene Plasmids",  
        "recommended_prime_editor",
        *[element for element in prime_card_elements if element is not None]  # Filter out None elements
    )


def generate_prime_protocals_section(guides_df):
    # Extract the oligos from the DataFrame
    pegRNA_oligos = guides_df["PrimeDesign pegRNA Spacer Oligo Top"].tolist()
    ngRNA_oligos = guides_df["PrimeDesign ngRNA Oligo Top"].tolist()

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
    if pegRNA_oligos and any(pegRNA_oligos):
        processed_pegRNA_oligos = process_peg_rnas(pegRNA_oligos)
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
        ui.tags.a("BsaI enzyme.", href="https://www.neb.com/en-us/products/r3733-bsai-hf-v2", target="_blank"),
    ),
])
    prime_section.extend([
    ui.help_text(
        ui.br(),
        ui.br(),
        ui.tags.span("3. Follow Cloning Protocol here for "),
        ui.tags.a("pegRNA plasmid.", href="https://drive.google.com/file/d/1kKD7EVwS7nZMbiZ6UH4LftK7GepoS8si/view", target="_blank"),
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
        ) )
            #displaying the ng rna oligo top and bottom
            prime_section.append(ui.br())
            prime_section.append(ui.br())
            for guide_rna in processed_ngRNA_oligos:
                guide_rna_parts = guide_rna.split('\n')
                for part in guide_rna_parts:
                    prime_section.append(ui.help_text(part))
                    prime_section.append(ui.br())
                prime_section.append(ui.br())
                
                #adding cloning protocal for ng rna
                prime_section.extend([
                    ui.help_text(
                        ui.tags.span("2. Follow Cloning Protocol here for "),
                        ui.tags.a("ngRNA plasmid.", href="https://media.addgene.org/cms/filer_public/6d/d8/6dd83407-3b07-47db-8adb-4fada30bde8a/zhang-lab-general-cloning-protocol-target-sequencing_1.pdf", target="_blank"),
                    ),
                ])

    prime_section.append(ui.br())
    prime_section.append(ui.br())

    return ui_card("Experimental Validation of Prime Editing", 'prime_section', *prime_section)

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
            ui.tags.a("here.", href="https://media.addgene.org/cms/filer_public/6d/d8/6dd83407-3b07-47db-8adb-4fada30bde8a/zhang-lab-general-cloning-protocol-target-sequencing_1.pdf", target="_blank"),
        ),
    ])
    
    return ui_card("Experimental Validation of Base Editing Guide RNAs", 'validation_section', *validation_section)

def generate_base_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM):
    list_of_guides_to_display = []
    for index, row in guides_df.iterrows():
        if row['Editing Technology'] == 'Base Editing':
            guide = row["Base Editing Guide"]
            orientation = row["Base Editing Guide Orientation"]

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

            list_of_guides_to_display.append(ui.help_text(f"Guide {index + 1}"))
            list_of_guides_to_display.append(ui.br())
            if orientation == 'reverse':
                reverse_guide_rna_with_formatting = []
                for i in range(len(guide)):
                    base_position = guide_start + len(guide) - 1 - i
                    base = ref_sequence_almost_rc[base_position]
                    if i >= 3 and i < 9:  # Bases 4-9 (0-based index)
                        if base_position == substitution_position:  # Base to edit
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: green; font-family: Courier,courier"))
                        elif base == base_to_edit:
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: red; font-family: Courier,courier"))
                        else:
                            reverse_guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: violet; font-family: Courier,courier"))
                    else:
                        reverse_guide_rna_with_formatting.append(ui.tags.b(base, style="color: violet; font-family: Courier,courier"))

                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Forward Strand: 5'-" + ref_sequence_input[:substitution_position], style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_input[substitution_position], style="color: green; font-family: Courier,courier"),  # base to edit
                        ui.tags.b(ref_sequence_input[substitution_position + 1:] + "-3'", style="font-family: Courier,courier"),
                    )
                )
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Reverse Strand: 3'-" + ref_sequence_almost_rc[:guide_start - len(PAM)], style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_almost_rc[guide_start - len(PAM):guide_start], style="color: blue; font-family: Courier,courier"),  # PAM
                        *reverse_guide_rna_with_formatting[::-1],  # reverse guide RNA
                        ui.tags.b(ref_sequence_almost_rc[guide_start + len(guide):] + "-5'", style="font-family: Courier,courier"),
                    )
                )
            else:
                guide_rna_with_formatting = []
                for i in range(len(guide)):
                    base = ref_sequence_input[guide_start + i]
                    if i >= 3 and i < 9:  # Bases 4-9 (0-based index)
                        if guide_start + i == substitution_position:  # Base to edit
                            guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: green; font-family: Courier,courier"))
                        elif base == base_to_edit:
                            guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: red; font-family: Courier,courier"))
                        else:
                            guide_rna_with_formatting.append(ui.tags.b(base, style="text-decoration: underline; color: violet; font-family: Courier,courier"))
                    else:
                        guide_rna_with_formatting.append(ui.tags.b(base, style="color: violet; font-family: Courier,courier"))

                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Forward Strand: 5'-" + ref_sequence_input[:guide_start], style="font-family: Courier,courier"),
                        *guide_rna_with_formatting,
                        ui.tags.b(ref_sequence_input[len(guide) + guide_start:len(guide) + guide_start + len(PAM)], style="color: blue; font-family: Courier,courier"),  # PAM
                        ui.tags.b(ref_sequence_input[guide_start + len(guide) + len(PAM):] + "-3'", style="font-family: Courier,courier"),
                    )
                )
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(
                    ui.help_text(
                        ui.tags.b("Reverse Strand: 3'-" + ref_sequence_almost_rc[:substitution_position], style="font-family: Courier,courier"),
                        ui.tags.b(ref_sequence_almost_rc[substitution_position], style="color: green; font-family: Courier,courier"),  # base to edit
                        ui.tags.b(ref_sequence_almost_rc[substitution_position + 1:] + "-5'", style="font-family: Courier,courier"),
                    )
                )
            if index != guides_df.shape[0] - 1:
                list_of_guides_to_display.append(ui.br())
                list_of_guides_to_display.append(ui.br())
    return list_of_guides_to_display

def server(input, output, session):
    
    def input_check(ref_sequence_input, edited_sequence_input):
        nonlocal input_file
        if input_file and not (ref_sequence_input or edited_sequence_input):
            try:
                df = pd.read_csv(input_file)
            except:
                return False, "Input file is not a properly formed CSV file. Please input a proper CSV file."

            if len(df.columns) != 2 or df.columns.tolist() != ['Original Sequence', 'Desired Sequence']:
                return False, 'Uploaded csv does not have the proper columns. Your csv must have two columns with names "Original Sequence" and "Desired Sequence"'
            
            counter = 1
            for index, row in df.iterrows():
                ref_sequence = "".join(row['Original Sequence'].split()).upper()
                edited_sequence = "".join(row['Desired Sequence'].split()).upper()
                check, message = check_ref_edited_pair(ref_sequence, edited_sequence)
                if not check:
                    return check, f"Error row {counter}: {message}"
                counter += 1
            return True, "Input CSV verified. Proceed to get guides."
        elif ref_sequence_input and edited_sequence_input and not input_file:
            check, message = check_ref_edited_pair("".join(ref_sequence_input.split()).upper(), "".join(edited_sequence_input.split()).upper())
            return check, message
        elif ref_sequence_input and edited_sequence_input and input_file:
            return False, "Error: Fill in both text input fields or upload a CSV file but do not do both."
        else:
            return False, "Error: Fill in both text input fields or upload a CSV file."
    
    input_file = None

    @output
    @render.ui
    @reactive.event(input.upload)
    def upload():
        nonlocal input_file
        file_infos = input.file1()
        if file_infos:
            input_file = file_infos[0]['datapath']
            return ui.div(ui.br(), ui.tags.b("File Successfully Uploaded", style="color: grey;", id='upload_status'))
        else:
            return ui.div(ui.br(), ui.tags.b("Error: No file selected", style="color: red;", id='upload_status'))

    @reactive.Effect()
    def clear():
        value = input.clear()    
        if value > 0:
            ui.update_text_area("ref_sequence_input", value = "")
            ui.update_text_area("edited_sequence_input", value = "")
            ui.update_select("pam_type", selected='NGN')
            ui.remove_ui(selector="div:has(> #results)")
            ui.remove_ui(selector="div:has(> #upload_status)")

            nonlocal input_file
            input_file = None

    @output
    @render.ui
    def ui_input_file():
        input.clear()  
        return ui.input_file(f"file1", 'Choose a CSV File of Sequences to Upload (note that you must click the blue "Upload File" button even if the progress bar under the file browser says "Upload complete". Also, clicking the button will cause the screen to scroll up to the top which is annoying and we are trying to fix that):', accept='.csv', multiple=False, width="100%"),
        
    @output
    @render.ui
    @reactive.event(input.get_guides)
    def run():
        ref_sequence_input = "".join(input.ref_sequence_input().split()).upper()
        edited_sequence_input = "".join(input.edited_sequence_input().split()).upper()
        
        # Call check_ref_edited_pair directly
        valid, message = check_ref_edited_pair(ref_sequence_input, edited_sequence_input)
        if not valid:
            return ui.div(ui.tags.b(f"Error: {message}", style="color: red;"))
        
        selected_PAM = input.pam_type()
        
        to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input, selected_PAM)
        base_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Base Editing']
        prime_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Prime Editing']
        
        message = check_for_BE_or_PE_guides(guides_df)
        if message: 
            return ui.div(ui.tags.b(f"{message}", style="color: red;"), ui.tags.a("https://www.nature.com/articles/nprot.2013.143", href=" https://www.nature.com/articles/nprot.2013.143", target="_blank", style="color: red;"),)
        
        if not base_editing_guides_df.empty:
            pams = guides_df["Base Editing Guide Pam"].tolist()
            PAM = pams[0]
        else:
            PAM = selected_PAM

         #defining variables to find enzyme protein on addgene in suggested addgene plasmid section
        editor_name, editor_id, editor_url = get_editor_info(ref_sequence_input, edited_sequence_input, PAM)
        editor_info = f"{editor_name} (Addgene: {editor_id})"
        
        #defining variables to find guide rna cloning plasmid on addgene in suggested addgene plasmid section
        cloning_name, cloning_id, cloning_url = get_cloning_url(PAM)
        plasmid_info = f"{cloning_name} (Addgene: {cloning_id})"

        @output
        @render.data_frame
        def render_results():
            nonlocal to_display_guides_df
            return render.DataGrid(
                to_display_guides_df,
                width="100%",
                filters=False,
                summary=True,
            )

        @session.download(filename=lambda: f"guides-{date.today().isoformat()}-{datetime.now().strftime('%H-%M-%S')}.csv")
        async def download_results():
            nonlocal guides_df
            yield guides_df.to_csv()

        nonlocal input_file
        valid_inputs, message = input_check(ref_sequence_input, edited_sequence_input)

        if valid_inputs:
            if input_file and not (ref_sequence_input or edited_sequence_input):
                df = pd.read_csv(input_file)
                dfs_to_merge_download = list()
                dfs_to_merge_display = list()
                counter = 1

                with ui.Progress(min=1, max=df.shape[0] + 1) as p:
                    p.set(message="Finding guides", detail="This may take a while...")
                    for index, row in df.iterrows():
                        p.set(counter, message="Finding guides")
                        ref_sequence_input = "".join(row['Original Sequence'].split()).upper()
                        edited_sequence_input = "".join(row['Desired Sequence'].split()).upper()
                        to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input, PAM)
                        index_column = [str(counter)] * to_display_guides_df.shape[0]
                        to_display_guides_df.insert(loc=0, column='Input CSV Row Number', value=index_column)
                        dfs_to_merge_download.append(guides_df)
                        dfs_to_merge_display.append(to_display_guides_df)
                        counter += 1

                to_display_guides_df = pd.concat(dfs_to_merge_display)
                to_display_guides_df = to_display_guides_df.drop(columns=['Original Sequence', 'Desired Sequence'])
                guides_df = pd.concat(dfs_to_merge_download)
                base_editing_guides_df = guides_df[guides_df['Editing Technology'] == 'Base Editing']
                
                ui_elements = [
                    ui.help_text(
                        "Note: for base editing, there can be more than one guide RNA for a single desired edit, but for prime editing, we will only show the recommended PrimeDesign guide RNA"
                    ),
                    ui.br(),
                    ui.br(),
                    ui.output_data_frame("render_results"),
                    ui.br(),
                    ui.br(),
                ]

                ui_elements.append(ui.download_button("download_results", "Download Results as CSV File"))

                return ui.TagList(
                    ui_card(
                        "Recommended Guide RNAs",
                        'results',
                        *ui_elements
                    )
                )

            elif ref_sequence_input and edited_sequence_input and not input_file:
                ref_sequence_input = "".join(ref_sequence_input.split()).upper()
                edited_sequence_input = "".join(edited_sequence_input.split()).upper()
                to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input, PAM)
                to_display_guides_df = to_display_guides_df.drop(columns=['Original Sequence', 'Desired Sequence'])
                to_display_guides_df.insert(loc=0, column='Guide', value=[f"Guide {i + 1}" for i in range(to_display_guides_df.shape[0])])

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

                list_of_guides_to_display = generate_base_editing_visualization(guides_df, ref_sequence_input, substitution_position, PAM)

                # Filter for Base Editing or Prime Editing guides only
                filtered_guides_df = to_display_guides_df[to_display_guides_df['Editing Technology'].isin(['Base Editing', 'Prime Editing'])]

                ui_elements = [
                    ui.help_text(
                        '''Note: for base editing, there can be more than one guide RNA for a single desired edit, but for prime editing, we will only show the recommended PrimeDesign guide RNA. If your selected pam is unavailible, see NRN or NYN pams for base editing. The off-target score is calculated using the CFD score algorithm (Doench et al. 2014) where a higher score indicates a lower likelihood of off-target activity (higher is better). The on-target score is calculated using the RuleSet1 algorithm where a higher score indicates greater efficiency of guide RNA binding to the genomic target sequence (higher is better). Both algorithms are from ''',
                        ui.tags.a('Doench et al. 2014 Nat Biotechnol.', {'href': 'https://pubmed.ncbi.nlm.nih.gov/25184501/', 'target': '_blank'})
                    ),
                    ui.br(),
                    ui.br(),
                ]
                
                message = check_availible_pam(PAM, selected_PAM)
                if message:
                    ui_elements.append(ui.help_text(ui.tags.b(f"{message}", style="color: red;")))
                    ui_elements.append(ui.br(),)
                    
                ui_elements.append(
                    ui.output_data_frame("render_results"),
                )
                ui_elements.append(ui.br(),)
                
                if 'Base Editing' in filtered_guides_df['Editing Technology'].values:
                    ui_elements.append(
                        ui_card(
                            "Suggested Addgene Plasmids",  
                            "recommended_editor",
                            ui.help_text(f"Recommended Base Editing Plasmid: {editor_info}"),
                            ui.br(),
                            ui.help_text(f"Recommended Guide RNA Plasmid: {plasmid_info}"),
                            ui.br(),
                            ui.br(),
                            ui.tags.a("View Base Editing Plasmid", href=editor_url, target="_blank", class_="btn btn-primary"),
                            ui.help_text(" "),
                            ui.tags.a("View Guide RNA Plasmid", href=cloning_url, target="_blank", class_="btn btn-primary"),
                        )
                    )
                    
                if 'Prime Editing' in filtered_guides_df['Editing Technology'].values:
                    prime_editing_plasmid_card = create_prime_editing_plasmid_card(guides_df, editor_info, editor_url)
                    ui_elements.append(prime_editing_plasmid_card)


                if 'Base Editing' in filtered_guides_df['Editing Technology'].values:
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
                                " characters represent a span for potential bystander edits, and ",
                                ui.tags.b("red", style= "color: red"),
                                " characters represent the potential bystander edits, themselves. Grey characters represent nucleotides not spanned by the guide. Only 25 bp of sequence upstream and downstream of the desired edit is shown. NOTE: when toggling the order of the guides, the guide number above corresponds to the guide number within the visualization."
                            ),
                            ui.br(),
                            ui.br(),
                            *list_of_guides_to_display
                        )
                    )
                    ui_elements.append(generate_experimental_validation_section(base_editing_guides_df, PAM))
                    
                if 'Prime Editing' in filtered_guides_df['Editing Technology'].values:
                    ui_elements.append(generate_prime_protocals_section(prime_editing_guides_df))
                    
                ui_elements.append(ui.download_button("download_results", "Download Results as CSV File"))

                return ui.TagList(
                    ui_card(
                        "Recommended Guide RNAs",
                        'results',
                        *ui_elements
                    )
                )
        else:
            return ui.div(ui.tags.b(message, style="color: red;"))

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
