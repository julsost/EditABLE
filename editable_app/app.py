import asyncio
import pandas as pd
from datetime import date
from datetime import datetime
from utils import get_guides
from shiny import App, render, ui, reactive

import numpy as np

bases = {"A", "C", "G", "T"}
accepted_bases = {"A", "C", "G", "T", "-"}

# A card component wrapper.
def ui_card(title, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )

app_ui = ui.page_fluid(
    {"id": "main-content"},
    ui.panel_title('editABLE'),
    ui.br(),
    ui.help_text(
        '''Welcome to editABLE! This Shiny app was built primarily to find CRISPR base editing guides, 
         but can also do a first pass analysis to find prime editing guides.'''),
    ui.br(),
    ui.br(),
    ui_card(
        "How to use this app",
        ui.help_text(
            '''In CRISPR editing experiments, one is trying to induce some change in a DNA sequence.
             Therefore, you have a reference sequence (the original sequence you are trying to change)
             and a edited sequence (what you want your sequence to look like after the CRISPR edit).'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("There are two ways to use this app:"),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''1. If you want to find guides for a single CRISPR edit. For this use case, please enter 
            in your reference sequence and edited sequence in their respective input boxes.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''2. If you have more than one CRISPR edit you want to make, you can upload a CSV file with 
            two columns, named "Reference Sequence" and "Edited Sequence" that contain your reference and 
            edited sequences, with each row representing one edit you would like to make.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''After specifing your input, click the "Find Guides" button and editABLE will try to find CRISPR 
            guides that will induce your desired edits. Once editABLE finishes running, a table will appear 
            displaying either the guides that editABLE has found for each of your desired edits or a suggestion to 
            use an alternative CRISPR technology if base or prime editing guides can't be found. Lastly, you can 
            download a CSV of the guides found by editABLE by clicking on the "Download Results as CSV File" 
            button. Note that if you uploaded a large CSV file with many edits, you may have to wait awhile.'''
        )
    ),
    ui_card(
        "Input requirements",
        ui.help_text(
            '''Your reference sequence(s) and the edited sequence(s) must be the same length. Only single edits 
            (SNV, insertion, deletion) are supported at this time.'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For single nucleotide variant (SNV) edits, the input sequences are the most straightforward. You can 
            input the reference and edited sequences in without modification. For example, this would be a valid set of inputs:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Reference Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("A", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Edited Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("T", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For changes that result in deletions, use a string of "-" characters in the edited sequence to denote the 
            area of the deletion. For example, this would be a valid set of inputs for a deletion:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Reference Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("ATT", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Edited Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("---", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''For changes that result in insertions/duplications, use a string of "-" characters in the reference sequence 
            to denote the area of the insertion/duplication. For example, this would be a valid set of inputs for a insertions/duplications:'''
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Reference Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("---", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text("Edited Sequence:"),
        ui.br(),
        ui.help_text(
            ui.tags.b("GATAGCTCAGCTAGCCTAGTCAAACCTATC", style="font-family: Courier,courier"), 
            ui.tags.b("GCG", style="color: red; font-family: Courier,courier"), 
            ui.tags.b("ACGTCGATCGATCGATCACACCGCCTAATC", style="font-family: Courier,courier"),
        ),
        ui.br(),
        ui.br(),
        ui.help_text(
            '''Lastly, we require at least 25 base pairs of sequence to the left and right of your desired edit. 
            So in each of the examples above, there must be 25 or more base pairs to the right and left 
            of the red highlighted regions.'''
        )
    ),
    ui_card(
        "Input",
        ui.input_text_area("ref_sequence_input", "Reference Sequence", placeholder="Enter sequence", height="50%", width="100%"),
        ui.input_text_area("edited_sequence_input", "Edited Sequence", placeholder="Enter sequence", height="50%", width="100%"),
        #ui.tags.a(ui.input_file("file1", "Choose a csv file of sequences to upload:", accept='.csv', multiple=False, width="100%"), {"position" : "sticky"}),
        ui.input_file("file1", "Choose a CSV File of Sequences to Upload (note that clicking the button will cause the screen to scroll up to the top which is annoying and we are trying to fix that):", accept='.csv', multiple=False, width="100%"),
        ui.input_action_button("get_guides", "Find Guides", class_="btn-success"),
    ),
    ui_card(
        "Results",
        ui.help_text(
            "Note: there can be more than one base editing guide for a single desired edit, but we will only show the recommended PrimeDesign prime editing guide"
        ),
        ui.br(),
        ui.br(),
        ui.output_ui("run")
    ),
    ui.help_text(
        '''Note that we use PrimeDesign (Hsu, J.Y., Grünewald, J., Szalay, R. et al. 
        PrimeDesign software for rapid and simplified design of prime editing guide RNAs. 
        Nat Commun 12, 1034 (2021)) for prime editing guide calcualtions. We run PrimeDesign 
        with default parameters and take only the suggested guides. For more advanced usage 
        please use the ''',
        ui.tags.a('PrimeDesign portal', {'href' : 'https://primedesign.pinellolab.partners.org/', 'target' : '_blank'}),
        ''' (https://primedesign.pinellolab.partners.org/) to design your Prime Editing guides.'''
    ),
)

def check_ref_edited_pair(ref_sequence, edited_sequence):
    if len(ref_sequence) == 0 or len(edited_sequence) == 0:
        return False, "Both the reference sequence and the edited sequence must be of nonzero length."
    if len(ref_sequence) != len(edited_sequence):
        return False, "The length of the reference sequence and the edited sequence must be the same."
    if ref_sequence == edited_sequence:
        return False, "The reference sequence and the edited sequence are the same."
    if len(set(ref_sequence) - accepted_bases) != 0 or len(set(edited_sequence) - accepted_bases) != 0:
        return False, "You may only have the following characters in your sequences {A, C, G, T, -}."
    if len(set(ref_sequence) - bases) == 0 and len(set(edited_sequence) - bases) == 0:
        substitution_position = None
        for i in range(len(ref_sequence)):
            ref_base = ref_sequence[i]
            edited_base = edited_sequence[i]
            if ref_base != edited_base:
                if substitution_position is not None:
                    return False, "The reference sequence and the edited sequence contain more than one SNV. You can only have one SNV edit."
                else:
                    substitution_position = i
        if substitution_position < 25:
            return False, f"There must be at least 25 base pairs of sequence before the desired edit. {substitution_position} base pairs were found before your edit."
        if len(ref_sequence) - 1 - substitution_position < 25:
            return False, f"There must be at least 25 base pairs of sequence after the desired edit. {len(ref_sequence) - 1 - substitution_position} base pairs were found after your edit."
    else:
        if '-' not in ref_sequence and '-' not in edited_sequence:
            return False, 'The lengths of the reference sequence and the edited sequence are not the same but neither has a "-" in it.'
        if '-' in ref_sequence and '-' in edited_sequence:
            return False, 'You cannot have a "-" in both the reference and edited sequences.'
        elif '-' in ref_sequence:
            start_dash_position = None
            current_dash_position = None
            for i in range(len(ref_sequence)):
                if ref_sequence[i] == '-':
                    if start_dash_position is None:
                        start_dash_position = i
                    if current_dash_position is not None and i - current_dash_position != 1:
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple insertions. You can only have one insertion.'
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
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple insertions. You can only have one insertion.'
                    else:
                        current_dash_position = i
            if start_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence before the desired edit. {start_dash_position} base pairs were found before your edit."
            if len(ref_sequence) - 1 - current_dash_position < 25:
                return False, f"There must be at least 25 base pairs of sequence after the desired edit. {len(ref_sequence) - 1 - current_dash_position} base pairs were found after your edit."
    return True, "Inputs verified. Proceed to get guides."
    
def server(input, output, session):    
    MAX_SIZE = 1000
    
    def input_check(ref_sequence_input, edited_sequence_input, file_infos):
        if file_infos and not (ref_sequence_input or edited_sequence_input):
            uploaded_fule = file_infos[0]
            try:
                df = pd.read_csv(uploaded_fule['datapath'])
            except:
                return False, "Input file is not a properly formed CSV file. Please input a proper CSV file."
            #else:
            if len(df.columns) != 2 and df.columns.tolist() != ['Reference Sequence', 'Edited Sequence']:
                return False, 'Uploaded csv does not have the proper columns. Your csv must have two columns with names "Reference Sequence" and "Edited Sequence"'
            counter = 1
            for index, row in df.iterrows():
                ref_sequence = row['Reference Sequence'].strip()
                edited_sequence = row['Edited Sequence'].strip()
                check, message = check_ref_edited_pair(ref_sequence, edited_sequence)
                if not check:
                    return check, f"Error row {counter}: {message}"
                counter += 1
            return True, "Input CSV verified. Proceed to get guides."
        elif ref_sequence_input and edited_sequence_input and not file_infos:
            check, message = check_ref_edited_pair(ref_sequence_input.strip(), edited_sequence_input.strip())
            return check, message
        else:
            return False, "Error: Fill in both text input fields or upload a CSV file but do not do both."
    
    @output
    @render.ui
    @reactive.event(input.get_guides)
    def run():
        
        @output
        @render.data_frame
        def render_results():
            nonlocal to_display_guides_df
            return render.DataGrid(
                        to_display_guides_df,
                        row_selection_mode='none',
                        width="100%",
                        filters=False,
                        summary = True,
            )

        @session.download(filename=lambda: f"guides-{date.today().isoformat()}-{datetime.now().strftime('%H-%M-%S')}.csv")
        async def download_results():
            nonlocal guides_df
            yield guides_df.to_csv()
            
        ref_sequence_input = input.ref_sequence_input()
        edited_sequence_input = input.edited_sequence_input()
        file_infos = input.file1()

        valid_inputs, message = input_check(ref_sequence_input, edited_sequence_input, file_infos)

        if valid_inputs:
            if file_infos and not (ref_sequence_input or edited_sequence_input):
                uploaded_file = file_infos[0]
                df = pd.read_csv(uploaded_file['datapath'])
                dfs_to_merge_download = list()
                dfs_to_merge_display = list()
                counter = 1
                
                with ui.Progress(min=1, max=df.shape[0]) as p:
                    p.set(message="Finding guides", detail="This may take a while...")
                    for index, row in df.iterrows():
                        p.set(counter, message="Finding guides")
                        ref_sequence_input = row['Reference Sequence'].strip()
                        edited_sequence_input = row['Edited Sequence'].strip()
                        guides_df = get_guides(ref_sequence_input, edited_sequence_input)
                        to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input)
                        index_column = [str(counter)] * to_display_guides_df.shape[0]
                        to_display_guides_df.insert(loc=0, column='Input CSV Row Number', value=index_column)
                        dfs_to_merge_download.append(guides_df)
                        dfs_to_merge_display.append(to_display_guides_df)
                        counter += 1

                to_display_guides_df = pd.concat(dfs_to_merge_display)
                to_display_guides_df = to_display_guides_df.drop(columns=['Original Sequence', 'Edited Sequence'])
                guides_df = pd.concat(dfs_to_merge_download)
            elif ref_sequence_input and edited_sequence_input and not file_infos:
                to_display_guides_df, guides_df = get_guides(ref_sequence_input.strip(), edited_sequence_input.strip())
                to_display_guides_df = to_display_guides_df.drop(columns=['Original Sequence', 'Edited Sequence'])
            return ui.TagList(
                ui.output_data_frame("render_results"),
                ui.br(),
                ui.br(),
                ui.download_button("download_results", "Download Results as CSV File")
            )
        else:
            return ui.div(ui.tags.b(message, style="color: red;"))

app = App(app_ui, server)
