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
    ui.h2("Guide Finder"),
    ui.help_text('Please enter in the original reference sequence that you will start with and the edited sequence that contains the edit that you wish to make. The reference sequence and the edited sequences must be the same length. Only single edits (SNV, insertion, deletion) are supported. For deletions, use a string of "-" characters in the edited sequence to denote the deletion. For insertions/duplications, use a string of "-" characters in the reference sequence to denote the insertion/duplication. You can also upload a CSV file with two columns, named "Reference Sequence" and "Edited Sequence" that contain your reference and edited sequences.'),
    ui_card(
        "Input",
        ui.input_text_area("ref_sequence_input", "Reference Sequence", placeholder="Enter text", height="50%", width="100%"),
        ui.input_text_area("edited_sequence_input", "Edited Sequence", placeholder="Enter text", height="50%", width="100%"),
        ui.input_file("file1", "Choose a csv file of sequences to upload:", multiple=False, width="100%"),
        ui.input_action_button("get_guides", "Find Guides", class_="btn-success"),
    ),
    ui_card(
        "Results",
        ui.help_text('Please note that guide finding may take some time'),
        ui.output_ui("run")
    ),
    ui.help_text('Note that we use PrimeDesign (Hsu, J.Y., Grünewald, J., Szalay, R. et al. PrimeDesign software for rapid and simplified design of prime editing guide RNAs. Nat Commun 12, 1034 (2021)) for prime editing guide calcualtions. We run PrimeDesign with default parameters and take only the suggested guides. For more advanced usage please use the PrimeDesign portal to design your Prime Editing guides.'),
    
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
    else:
        if '-' not in ref_sequence and '-' not in edited_sequence:
            return False, 'The lengths of the reference sequence and the edited sequence are not the same but neither has a "-" in it.'
        if '-' in ref_sequence and '-' in edited_sequence:
            return False, 'You cannot have a "-" in both the reference and edited sequences.'
        elif '-' in ref_sequence:
            current_dash_position = None
            for i in range(len(ref_sequence)):
                if ref_sequence[i] == '-':
                    if current_dash_position is not None and i - current_dash_position != 1:
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple insertions. You can only have one insertion.'
                    else:
                        current_dash_position = i
        else:
            current_dash_position = None
            for i in range(len(edited_sequence)):
                if ref_sequence[i] == '-':
                    if current_dash_position is not None and i - current_dash_position != 1:
                        return False, 'The "-" characters are not contiguous, indicating that there are multiple insertions. You can only have one insertion.'
                    else:
                        current_dash_position = i
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
                ref_sequence = row['Reference Sequence']
                edited_sequence = row['Edited Sequence']
                check, message = check_ref_edited_pair(ref_sequence, edited_sequence)
                if not check:
                    return check, f"Error row {counter}: {message}"
                counter += 1
            return True, "Input CSV verified. Proceed to get guides."
        elif ref_sequence_input and edited_sequence_input and not file_infos:
            check, message = check_ref_edited_pair(ref_sequence_input, edited_sequence_input)
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
                uploaded_fule = file_infos[0]
                df = pd.read_csv(uploaded_fule['datapath'])
                #dfs_to_merge = list()
                dfs_to_merge_download = list()
                dfs_to_merge_display = list()
                for index, row in df.iterrows():
                    ref_sequence_input = row['Reference Sequence']
                    edited_sequence_input = row['Edited Sequence']
                    guides_df = get_guides(ref_sequence_input, edited_sequence_input)
                    #dfs_to_merge.append(guides_df)
                    to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input)
                    dfs_to_merge_download.append(guides_df)
                    dfs_to_merge_display.append(to_display_guides_df)

                to_display_guides_df = pd.concat(dfs_to_merge_display)
                guides_df = pd.concat(dfs_to_merge_download)
            elif ref_sequence_input and edited_sequence_input and not file_infos:
                to_display_guides_df, guides_df = get_guides(ref_sequence_input, edited_sequence_input)
            return ui.TagList(
                ui.output_data_frame("render_results"),
                ui.download_button("download_results", "Download Results as CSV File")
            )
        else:
            return ui.div(ui.tags.b(message, style="color: red;"))

app = App(app_ui, server)
