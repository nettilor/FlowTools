import PySimpleGUI as sg
import os
import sys
import warnings

warnings.filterwarnings("ignore")
sys.stderr = open(os.devnull, 'w')
# Set the theme to match macOS better
sg.theme('SystemDefault')

# Default values
default_workspace = ""
default_fcs_path = ""
default_groups = ""
default_keywords = ""
default_gates_fluoro = ""
default_output = ""

# Define the layout for the UI
layout = [
    [sg.Text("FlowJo Workspace (.wsp):", background_color='white'), 
     sg.Input(default_text=default_workspace, key="-WORKSPACE-", background_color='white'), 
     sg.FileBrowse(file_types=(("WSP Files", "*.wsp"),))],
    [sg.Text("FCS Files Folder:", background_color='white'),
     sg.Input(default_text=default_fcs_path, key="-FCS-FOLDER-", background_color='white'), 
     sg.FolderBrowse()],
    [sg.Text("Output CSV File:", background_color='white'), 
     sg.Input(default_text=default_output, key="-OUTPUT-", background_color='white'), 
     sg.FileSaveAs(file_types=(("CSV Files", "*.csv"),))],
    [sg.Text("Selected Groups (comma-separated):", background_color='white'), 
     sg.Input(default_text=default_groups, key="-GROUPS-", background_color='white')],
    [sg.Text("Keywords (comma-separated):", background_color='white'), 
     sg.Input(default_text=default_keywords, key="-KEYWORDS-", background_color='white')],
    [sg.Text("Gates and Fluorophores (Gate1: Fluo1, Fluo2...):", background_color='white')],
    [sg.Multiline(default_text=default_gates_fluoro, key="-GATES-FLUORO-", size=(50, 5), background_color='white')],
    [sg.Text("MFI Statistics:", background_color='white')],
    [sg.Checkbox('Mean', key='-MEAN-', background_color='white', default=False),
     sg.Checkbox('Median', key='-MEDIAN-', background_color='white', default=False),
     sg.Checkbox('Geometric Mean', key='-GEO_MEAN-', background_color='white', default=False)],
    [sg.Button("Run"), sg.Button("Exit")]
]

# Create the window with white background
window = sg.Window("FlowJo data extractor 🔧", layout, background_color='white')

# Event loop
while True:
    event, values = window.read()
    if event == sg.WINDOW_CLOSED or event == "Exit":
        window.close()
        sys.exit()   # Terminate the entire script
    if event == "Run":
        # Check if any required fields are blank
        if not values["-WORKSPACE-"] or not values["-FCS-FOLDER-"] or not values["-OUTPUT-"] or not values["-GROUPS-"]:
            sg.popup_ok("All required fields must be filled out!", title="Error", background_color='white')
        else:
            fj_workspace_path = values["-WORKSPACE-"]
            fcs_path = values["-FCS-FOLDER-"]
            output_file = values["-OUTPUT-"]
            selected_groups = [group.strip() for group in values["-GROUPS-"].split(",")]
            fields = [keyword.strip() for keyword in values["-KEYWORDS-"].split(",") if keyword.strip()]
            
            # Process the Gates and Fluorophores input
            gates_fluoro_input = values["-GATES-FLUORO-"].strip()
            list_mfi = []
            if gates_fluoro_input:
                for line in gates_fluoro_input.split("\n"):
                    if ':' in line:
                        gate, fluorophores = line.split(":", 1)
                        gate = gate.strip()
                        fluorophores = [f.strip() for f in fluorophores.split(",") if f.strip()]
                        list_mfi.append((gate, fluorophores))
                
                # Only check for statistics if Gates and Fluorophores is not empty
                statistics = []
                if values['-MEAN-']:
                    statistics.append('mean')
                if values['-MEDIAN-']:
                    statistics.append('median')
                if values['-GEO_MEAN-']:
                    statistics.append('geo_mean')
                
                if not statistics:
                    sg.popup_error("Error", "At least one MFI statistic must be selected when Gates and Fluorophores are provided.", background_color='white')
                    continue  # Go back to the main loop without closing the window
            else:
                statistics = []  # Empty list if no Gates and Fluorophores are provided
            
            # Close the window
            window.close()
            
            # Print all variables to the terminal
            print("\n--- Variables ---")
            print(f"FlowJo Workspace Path: {fj_workspace_path}")
            print(f"FCS Files Path: {fcs_path}")
            print(f"Output File: {output_file}")
            print(f"Selected Groups: {selected_groups}")
            print(f"Keywords: {fields}")
            print(f"Gates and Fluorophores: {list_mfi}")
            print(f"MFI Statistics: {statistics}")
            print("--- End of Variables ---\n")
            
            # Exit the loop
            break

window.close()

### DO NOT CHANGE ANYTHING BELOW THIS LINE, Thank you :) ###
############################################################

# Define accessory functions

def get_fcs_file_list(fcs_dir):
    """
    Get a list of all the fcs files in the directory
    """
    import os
    fcs_files=[]
    for root, dirs, files in os.walk(fcs_dir):
        for file in files:
            if file.endswith('.fcs') and "._" not in file:
                fcs_files.append(os.path.join(root, file))
    return fcs_files

#-----------------------------------------
# Import necessary packages

import numpy as np
import flowkit as fk
import pandas as pd
import scipy as sp

#-----------------------------------------
# Read in FlowJo workspace and FCS files
print('Reading FlowJo workspace...', end='\r')
## Build file list
fcs_file_list = get_fcs_file_list(fcs_path)
## Read FlowJo workspace and FCS files using FlowKit
wsp = fk.Workspace(fj_workspace_path ,fcs_samples=fcs_file_list)
print("                                          ", end='\r')
print("Reading FlowJo workspace ✅")

print('Compiling results from FlowJo gating...', end='\r')
## Retrieve FlowJo analysis from workspace file
all_groups = wsp.get_sample_groups()
#selected_groups = ['Blood samples']
selected_groups = [group for group in all_groups if group in selected_groups]
# Yields a dataframe with all gates for all samples in all groups selected
group_results_report = pd.DataFrame()
for group in selected_groups:
    wsp.analyze_samples(group_name=group, verbose=False, use_mp=True)
    results_report = wsp.get_analysis_report(group)
    group_results_report = pd.concat([group_results_report, results_report])

## Set the index of the dataframe to the sample names
group_results_report = group_results_report.set_index('sample')
print("                                          ", end='\r')
print("Compiling results from FlowJo gating ✅")

#-----------------------------------------
# Handle metadata (keywords extraction)
if (fields != ['']):
    print('Extracting keywords from metadata...', end='\r')
    ## Extract the sample data
    sample_metadata = fk.extract_wsp_sample_data(fj_workspace_path)

    ## Make list of all samples in selected groups
    all_samples_list = []
    for group in selected_groups:
        all_samples_list += wsp.get_sample_ids(group_name=group)

    ## Extract keywords from the defined list at beginning of script
    metadata_df = pd.DataFrame()
    for sample in all_samples_list:
        keywords = sample_metadata[sample]["keywords"]
        values = [keywords.get(field, np.NaN) for field in fields ]
        temp = pd.DataFrame([values], columns=fields, index=[sample])
        metadata_df = pd.concat([metadata_df, temp])

    ## Merge the metadata with the group_results_report and add metadata columns to the group_results_report where the indicies match
    group_results_report = group_results_report.merge(metadata_df, left_index=True, right_index=True)
    print("                                          ", end='\r')
    print("Extracting keywords from metadata ✅")

#-----------------------------------------
# Compute MFIs
print('Computing MFIs...', end='\r')
## First of all this needs to be done per group because we assume all groups have same gates
#selected_groups = 'Blood samples'
# Now we can specify gates and channels of interest in a format that is Gate - [channels]
#NEED TO ADD ONE MORE FOR LOOP HERE IF MULTIPLE SAMPLES ARE SELECTED
for selected_group in selected_groups:
    sample_ids = wsp.get_sample_ids(group_name=selected_group)
    # add columns for MFI to df
    for statistic in statistics:
        i=0
        for elements in list_mfi:
            gate = elements[0]
            fluorophores = elements[1]
            if i==0:
                for fluo in fluorophores:
                    i=i+1
                    group_results_report[fluo + '_' + statistic] = np.nan
            # need to add gates and channels of interest in the loop
            for sample_id in sample_ids:
                selected_gates = group_results_report[group_results_report['gate_name'] == gate].loc[sample_id]
                if type(selected_gates) == pd.core.series.Series:
                    # this means there is only one gate that matches the selection in this sample
                    selected_gate = selected_gates.gate_name
                    selected_gate_path = selected_gates.gate_path
                    gate_events = wsp.get_gate_membership(sample_id, gate_name=selected_gate, gate_path=selected_gate_path)
                    # Grab a sample & its compensation and apply it
                    sample = wsp.get_sample(sample_id)
                    sample_comp = wsp.get_comp_matrix(sample_id)
                    sample.apply_compensation(sample_comp)
                    # consider only selected fluorophores and calculate mfi for each
                    mfis = []
                    for fluo in fluorophores:
                        # Apply the comp matrix & extract the compensated CD4 events
                        # using marker name
                        channel_index = pd.Index(sample.channels["pns"]).get_loc(fluo)
                        #Using fluo name
                        #channel_index = sample.get_channel_index(fluo)
                        comp_events = sample.get_channel_events(channel_index, source='comp')
                        # Apply the Boolean gate membership from above to filter for the gated events 
                        pos_comp_events = comp_events[gate_events]
                        if statistic == 'mean':
                            gate_mfi = pos_comp_events.mean()
                        elif statistic == 'median':
                            gate_mfi = np.median(pos_comp_events)
                        elif statistic == 'geo_mean':
                            gate_mfi = sp.stats.gmean(abs(pos_comp_events))
                        
                        # put mfi in place
                        group_results_report.loc[(group_results_report.index==sample_id) & (group_results_report['gate_name'] == gate),
                                                fluo + '_' + statistic] = gate_mfi
                    
                else:
                    # this means there are multiple gates that match the selection in this sample so do mfi for each of them
                    for index, row in selected_gates.iterrows():
                        # Take each gate and identify it by its unique gating_path (aka gating strategy)
                        selected_gate = row.gate_name
                        selected_gate_path = row.gate_path
                        gate_events = wsp.get_gate_membership(sample_id, gate_name=selected_gate, gate_path=selected_gate_path)
                        # Grab a sample & its compensation and apply it
                        sample = wsp.get_sample(sample_id)
                        sample_comp = wsp.get_comp_matrix(sample_id)
                        sample.apply_compensation(sample_comp)
                        # consider only selected fluorophores and calculate mfi for each
                        mfis = []
                        for fluo in fluorophores:
                            # Apply the comp matrix & extract the compensated events
                            # using marker name
                            channel_index = pd.Index(sample.channels["pns"]).get_loc(fluo)
                            # Using fluo name (i dont use this but it would get the index of entry for fluo)
                            #channel_index = sample.get_channel_index(fluo)
                            comp_events = sample.get_channel_events(channel_index, source='comp')
                            # Apply the Boolean gate membership from above to filter for the gated events 
                            pos_comp_events = comp_events[gate_events]
                            if statistic == 'mean':
                                gate_mfi = pos_comp_events.mean()
                            elif statistic == 'median':
                                gate_mfi = np.median(pos_comp_events)
                            elif statistic == 'geo_mean':
                                gate_mfi = sp.stats.gmean(abs(pos_comp_events))
                            # put mfi in place
                            group_results_report.loc[(group_results_report.index==sample_id) & (group_results_report['gate_name'] == gate) & (group_results_report['gate_path'] == selected_gate_path),
                                                    fluo + '_' + statistic] = gate_mfi
                        
print("                                               ", end='\r')
print("Computing MFIs ✅")
#-----------------------------------------
# Export processed data
print('Exporting results to csv...', end='\r')
# rename the first column to "Sample"
group_results_report.index.name = 'Sample'
 
group_results_report.to_csv(output_file)
print("                                          ", end='\r')
print("Exporting results to csv ✅")