from __future__ import print_function
import mdtraj as md
import numpy as np
import pandas as pd
from bokeh.io import output_file, save
from bokeh.plotting import figure, show, output_file, from_networkx
from bokeh.models import ColumnDataSource, PrintfTickFormatter,LabelSet
from bokeh.layouts import column, gridplot
from vmd import molecule
from bokeh.models import (BoxSelectTool, Circle, EdgesAndLinkedNodes, HoverTool, Label,
                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool, ResetTool)
import networkx as nx
from bokeh.palettes import Spectral4
import sys
import os

"Usage: python MD_Analysis_workflow.py <trajectoryFile(clickme.dtr)> <in.cms File> <output.hmtl> <number_of_residues>"
"Take care of the relative position in directory"

# =====================FUNCTIONS DEFINITION===========================


def create_bokeh_plot_percentage_of_ss(outputFile):
    p = figure(plot_width=950, plot_height=400)
    p.line(x='Frame',y='C%',line_width=2, color="gray",source=source,legend_label="C")
    p.line(x='Frame',y='H%',line_width=2,color="red",source=source,legend_label="H")
    p.line(x='Frame',y='310%',line_width=2,color="blue",source=source,legend_label="310")
    p.line(x='Frame',y='π%',line_width=2,color="black",source=source,legend_label="π")
    p.line(x='Frame',y='S%',line_width=2,color="pink",source=source,legend_label="S")
    p.line(x='Frame',y='B%',line_width=2,color="olive",source=source,legend_label="B")
    p.line(x='Frame',y='E%',line_width=2,color="yellow",source=source,legend_label="E")
    p.line(x='Frame',y='T%',line_width=2,color="purple",source=source,legend_label="T")
    p.xaxis.axis_label = "Frame number"
    p.yaxis.axis_label = "Percentage of SS"
    p.yaxis.ticker = [25, 50, 75, 100]
    p.legend.click_policy="hide"
    return p

def create_bokeh_plot_mean_percentage(outputFile):
    mean_percentage = figure(x_range=aminoacids,plot_width=950, plot_height=500)
    mean_percentage.vbar_stack(column_names[1::2], x='aminoacid', width=0.9, color=colors, source=source_2 ,legend_label=column_names[:-1:2])
    mean_percentage.legend.click_policy="hide"
    mean_percentage.yaxis[0].formatter = PrintfTickFormatter(format="%0.1f")
    label_opts = dict(
        x=0, y=0,
        x_units='screen', y_units='screen'
    )
    msg1 = 'Overall_Helicity: ' + str(overall_helicity)
    caption1 = Label(text=msg1,text_font_size = '8pt', **label_opts)
    mean_percentage.add_layout(caption1, 'below')
    msg2 = 'Helical_presence: ' + str(Helical_presence)
    caption2 = Label(text=msg2,text_font_size = '8pt', **label_opts)
    mean_percentage.add_layout(caption2, 'below')
    msg3 = 'Mean percentage of each type of SS accross all frames'
    caption3 = Label(text=msg3, **label_opts)
    mean_percentage.add_layout(caption3, 'above')
    msg5 = "H: α-helix. Min 3 residues, B: residue in isolated β-bridge (single pair beta-sheet hydrogen bond formation),"
    caption5 = Label(text=msg5,text_font_size = '8pt', **label_opts)
    mean_percentage.add_layout(caption5, 'below')
    msg7 =  "E: beta sheet in parallel and/or anti-parallel sheet conformation (extended strand). Min 2 residues"
    caption7 = Label(text=msg7,text_font_size = '8pt', **label_opts)
    mean_percentage.add_layout(caption7, 'below')
    msg6 =  "310: 3-helix (310 helix). Min 3 residues, π: 5 helix (π-helix). Min 5 residues, T: Hydrogen bonded turn, S: Bend, C: Random Coil "
    caption6 = Label(text=msg6,text_font_size = '8pt', **label_opts)
    mean_percentage.add_layout(caption6, 'below')
    return mean_percentage

# You can change the size of the plot by changinge the parameters plot_width and plot_height
def create_bokeh_network_interactions_plot(outputFile,graph,title):
    plot = Plot(plot_width=900, plot_height=900,
                x_range=Range1d(-1.1,1.1), y_range=Range1d(-1.1,1.1))
    plot.add_tools(HoverTool(), TapTool(), BoxSelectTool(), ResetTool())
    graph_renderer = from_networkx(graph, nx.circular_layout, scale=0.9, center=(0,0))
    graph_renderer.node_renderer.glyph = Circle(size=20, fill_color='color')
    graph_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="color", line_alpha=0.8, line_width=1)
    graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color="color", line_width=5)
    graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color="color", line_width=5)
    graph_renderer.selection_policy = NodesAndLinkedEdges()
    graph_renderer.inspection_policy = EdgesAndLinkedNodes()
    plot.renderers.append(graph_renderer)
    pos = graph_renderer.layout_provider.graph_layout
    x,y=zip(*pos.values())
    modifier_x = [i+0.02 if i>0  else i-0.1 for i in x]
    modifier_y = [i+0.04 if i>0  else i-0.12 for i in y]
    field = residues
    label_source = ColumnDataSource({'x': modifier_x,'y':modifier_y, 'field': residues})
    labels = LabelSet(x='x', y='y', text='field', source=label_source, text_font_size='1.1em')
    plot.renderers.append(labels)    
    label_opts = dict(
        x=0, y=0,
        x_units='screen', y_units='screen'
    )
    msg1 = 'sb : red, pc : blue, ps : pink, ts : yellow, vdw : orange, hbbb : green, hbsb : purple, hbss : brown'
    caption1 = Label(text=msg1,text_font_size = '8pt', **label_opts)
    plot.add_layout(caption1, 'below')
    msg2 = ''
    caption2 = Label(text=msg2,text_font_size = '8pt', **label_opts)
    plot.add_layout(caption2, 'below')
    msg4 = ''
    caption4 = Label(text=msg4,text_font_size = '8pt', **label_opts)
    plot.add_layout(caption4, 'above')
    msg3 = title
    caption3 = Label(text=msg3, **label_opts)
    plot.add_layout(caption3, 'above')
    return plot

# =====================DICTIONARIES DEFINITION===========================

DSSP_codes = {
    'H': 'α-helix. Min length 3 residues',
    'B': 'residue in isolated β-bridge. (single pair beta-sheet hydrogen bond formation)',
    'E': 'beta sheet in parallel and/or anti-parallel sheet conformation (extended strand). Min length 2 residues.',
    '310': '3-helix (310 helix). Min length 3 residues.',
    'π': '5 helix (π-helix).Min length 5 residues',
    'T': 'Hydrogen bonded turn',
    'S': 'Bend',
    'C': 'Random Coil'
}

bond_colors = {
    "sb" : "red",
    "pc" : "blue",
    "ps" : "pink",
    "ts" : "yellow",
    "vdw" : "orange",
    "hbbb" : "green",
    "hbsb" : "purple",
    "hbss" : "brown",
    "lwb" : "grey",
    "lwb2" : "black"
}
# Modified cinema colouring scheme https://www.bioinformatics.nl/~berndb/aacolour.html
aminoacids_dict = {
    'CYS': 'yellow',
    'ASP': 'red',
    'SER': 'green',
    'GLN': 'green',
    'LYS': 'blue',
    'TRP': 'magenta',
    'THR': 'green',
    'ASN': 'green',
    'PRO': 'brown',
    'PHE': 'magenta',
    'ALA': 'gray',
    'GLY': 'brown',
    'ILE': 'grey',
    'LEU': 'grey',
    'HIS': 'cyan',
    'HIE': 'cyan',
    'ARG': 'blue',
    'MET': 'grey',
    'VAL': 'grey',
    'GLU': 'red',
    'TYR': 'magenta'}

replacements = {
    'HIE': 'HIS'
}

colors=["red","olive","yellow","blue","black","purple","pink","gray"]

# =====================EXTRACT DESMOND OUTPUT FILE DATA TO PANDA DATAFRAMES===========================


assert len(sys.argv) > 2, "Usage: python MD_Analysis_workflow.py <trajectoryFile(clickme.dtr)> <in.cms File> <output.hmtl> <number_of_residues>"

dtrFile = sys.argv[1]
cmsFile = sys.argv[2]

# Create a topology(pdb) file from Desmonds output files(-in.cms) based on vmd-python module
molid = molecule.load("mae", cmsFile)
molecule.write(molid, "pdb", "topology.pdb", first=0)

# Convert Desmond trajectory files into a suitable format
trajectory = md.load_dtr(dtrFile, top='topology.pdb')
trajectory[:].save('tmp1.h5')
os.system('mdconvert tmp1.h5 -o tmp1.dcd')

# Load a trajectory
traj = md.load('tmp1.h5')

# Pick peptide atom indices
# Depending on the ions and water models you use, the values below have to be changed. This step cleans the trajectory from all het groups and solvent molecules
peptide_atoms_indices = traj.topology.select('not resname T4P and not resname NA and not resname CL')

# Clean traj from het and solvent
peptide_traj = traj.atom_slice(peptide_atoms_indices)
peptide_residues=peptide_traj.topology.residues

# Clean traj from het and solvent
peptide_traj = traj.atom_slice(peptide_atoms_indices[:-1])

# Calculate DSSP
dssp = md.compute_dssp(peptide_traj)
dssp_data = pd.DataFrame(dssp, columns=[str(i) for i in peptide_traj.topology.residues])
dssp_data = dssp_data.loc[:,~dssp_data.columns.duplicated()]
dssp_data_transposed = dssp_data.iloc[:,0:int(sys.argv[4])].T

# occurrence --->  The number of occurrences of a single SS type related to a particular amino acid across all frames of the MD Simulation
# dssp_data --->  The SS type of a particular amino acid across in a particular frame of the MD
occurrence = pd.DataFrame()
dssp_data.replace(['G','Ι'],['310','π'])
numb_aa = dssp_data.shape[1]
frames = dssp_data.shape[0]
for key in DSSP_codes:
    dssp_data[str(key)] = dssp_data.isin([key]).sum(1)
    dssp_data[str(key)+'%'] = (dssp_data[str(key)]/numb_aa)*100
    occurrence[str(key)] = dssp_data.iloc[:,0:numb_aa].isin([key]).sum(0)
    occurrence[str(key)+ ' %'] = (occurrence[str(key)]/frames)*100
occurrence = occurrence.iloc[:int(sys.argv[4]),:]
dssp_data['Frame'] = range(1,dssp_data.shape[0]+1)
source = ColumnDataSource(dssp_data.iloc[:,-17:])
occurrence['aminoacid']= occurrence.index.tolist()
source_2 = ColumnDataSource(occurrence)

# Get the column names = SS type (add [::2] or [1::2] for SS total %)
column_names = occurrence.columns.values.tolist()

# Get the row names = aminoacid names
aminoacids = occurrence['aminoacid']
mean_percentage = figure(x_range=aminoacids,plot_width=950, plot_height=500, title="Mean percentage of each type of SS accross all frames")
mean_percentage.vbar_stack(column_names[1::2], x='aminoacid', width=0.9, color=colors, source=source_2 ,legend_label=column_names[:-1:2])
mean_percentage.legend.click_policy="hide"
mean_percentage.yaxis[0].formatter = PrintfTickFormatter(format="%0.1f")

# Calculate helix related measures
number_of_frames = dssp_data.Frame.max() + 1
overall_helicity = (occurrence['H'].sum()/(int(sys.argv[4])*number_of_frames))*100
count = 0
for i in range(frames):
    if "H" in dssp_data_transposed[i].tolist():
        count+=1
Helical_presence = (count/number_of_frames)*100

# ==================== EXTRACT ALL IMPORTANT DYNAMIC CONTACTS FROM THE GETCONTACTS LIBRARY ============================

# Adjust the desired bond interactions to calculate - see GetContacts documentation
os.system('python ./getcontacts/get_dynamic_contacts.py --topology topology.pdb --trajectory tmp1.dcd --itypes sb hb pc --output peptide_result.tsv')

# =====================FILE EDITING===========================

# Remove the GetContacts output files first 2 lines
with open('peptide_result.tsv', 'r') as f:
    reader = f.readlines()[2:]
    reader.insert(0, 'Frame\tInteraction_type\tAtom_1\tAtom_2\n')
with open('peptide_result.tsv', 'w') as fout:
    fout.writelines(reader)

# Create filtered file (without hbbb and hbsb)
slined = []
with open("peptide_result.tsv", "r") as file:
    for line in file:
        sline = line.split("\t")
#   ADJUST WHAT KIND OF BONDS YOU WANT TO FILTER HERE!!! (based on the GetContacts nomenclature) 
        if sline[1] != 'hbbb' and sline[1] != 'hbsb' :
            "\t".join(sline)
            slined.append(line)
with open("filtered_peptide_result.tsv", "w") as filed:
    filed.writelines(slined)

# Import GetContacts results into a pandas Dataframe
contact_results = pd.read_csv('peptide_result.tsv', delimiter = '\t')
filtered_contact_results = pd.read_csv('filtered_peptide_result.tsv', delimiter = '\t')

#=========================== DATA EDITING ===============================

# Correct 3-letter residue name e.g. HIE --> HIS
contact_results.replace(replacements, regex=True, inplace=True)

# Make data suitable for plotting
contact_results['Residue_1'] =contact_results['Atom_1'].str[2:]
contact_results['Residue_2'] =contact_results['Atom_2'].str[2:]
contact_results['Atom_1'] = contact_results['Atom_1'].str[2:]
contact_results['Atom_2'] = contact_results['Atom_2'].str[2:]
contact_results['Residue_1'] =contact_results['Residue_1'].str.extract('([^:]*:[^:]*)')
contact_results['Residue_1'] =contact_results['Residue_1'].str.replace('[:]', '') 
contact_results['Residues'] =contact_results['Residue_1'].str[:3]
contact_results['Residue_2'] =contact_results['Residue_2'].str.extract('([^:]*:[^:]*)')
contact_results['Residue_2'] =contact_results['Residue_2'].str.replace('[:]', '') 
contact_results['residues_pairs'] = [i for i in zip(contact_results['Residue_1'],contact_results['Residue_2'])]
interacting_residues = contact_results.Residue_1.unique()
contact_results['Color'] =  contact_results['Interaction_type'].map(bond_colors)


# Correct 3-letter residue name e.g. HIE --> HIS
filtered_contact_results.replace(replacements, regex=True, inplace=True)

# Make data suitable for plotting
filtered_contact_results['Residue_1'] =filtered_contact_results['Atom_1'].str[2:]
filtered_contact_results['Residue_2'] =filtered_contact_results['Atom_2'].str[2:]
filtered_contact_results['Atom_1'] = filtered_contact_results['Atom_1'].str[2:]
filtered_contact_results['Atom_2'] = filtered_contact_results['Atom_2'].str[2:]
filtered_contact_results['Residue_1'] =filtered_contact_results['Residue_1'].str.extract('([^:]*:[^:]*)')
filtered_contact_results['Residue_1'] =filtered_contact_results['Residue_1'].str.replace('[:]', '') 
filtered_contact_results['Residues'] =filtered_contact_results['Residue_1'].str[:3]
filtered_contact_results['Residue_2'] =filtered_contact_results['Residue_2'].str.extract('([^:]*:[^:]*)')
filtered_contact_results['Residue_2'] =filtered_contact_results['Residue_2'].str.replace('[:]', '') 
filtered_contact_results['residues_pairs'] = [i for i in zip(filtered_contact_results['Residue_1'],filtered_contact_results['Residue_2'])]
interacting_residues = filtered_contact_results.Residue_1.unique()
filtered_contact_results['Color'] =  filtered_contact_results['Interaction_type'].map(bond_colors)

# Filter filtered_contact_results database based on desired conditions: 1) Show only highly occurring interactions(>29%)(based on atom interactions)
number_of_frames = contact_results.Frame.max() + 1
#   ADJUST PERCENTAGE OF PRESENCE OVER THE WHOLE SIMULATION HERE!!! 
alpha = pd.DataFrame((filtered_contact_results['residues_pairs'].value_counts()/number_of_frames > 0.20))
alpha.rename(columns = {'residues_pairs':'Counts'}, inplace = True)
alpha['residues_pairs'] = alpha.index
filtered_contacts = alpha[alpha['Counts'] == True]
beta = filtered_contacts.residues_pairs.tolist()

super_filtered_contact_results = contact_results[contact_results['residues_pairs'].isin(beta)]

# Extract Residues (in order) from MD traj objects
peptide_residues=peptide_traj.topology.residues
residues = [str(i) for i in peptide_residues]
for value in list(replacements.keys()):
    if contact_results["Residues"].str.contains(value).any():
        residues = [replacements[value]+x[3:] if x[:3]==value else x for x in residues]

# Utilize Networkx lib to create a network of interaction between residues
G=nx.MultiGraph()
for res in residues:
    G.add_node(res,color= aminoacids_dict[res[:3]])
for pair, c, a1, a2, f in zip(contact_results['residues_pairs'],contact_results['Color'],contact_results['Atom_1'],contact_results['Atom_2'],contact_results['Frame']):
    G.add_edge(pair[0],pair[1],color= c, atom_1=a1, atom_2 = a2, frames = f)

# Utilize Networkx lib to create a network of interaction between residues for filtered data
filtered_G=nx.MultiGraph()
for res in residues:
    filtered_G.add_node(res,color= aminoacids_dict[res[:3]])
for pair, c, a1, a2, f in zip(super_filtered_contact_results['residues_pairs'],super_filtered_contact_results['Color'],super_filtered_contact_results['Atom_1'],super_filtered_contact_results['Atom_2'],super_filtered_contact_results['Frame']):
    filtered_G.add_edge(pair[0],pair[1],color = c, atom_1 = a1, atom_2 = a2, frames = f)

output_file(sys.argv[3])

a = create_bokeh_plot_percentage_of_ss(sys.argv[3])
b = create_bokeh_plot_mean_percentage(sys.argv[3])
c = create_bokeh_network_interactions_plot(sys.argv[3],G,'Peptide Intramolecular Interactions')
d = create_bokeh_network_interactions_plot(sys.argv[3],filtered_G,'Peptide Intramolecular Interactions without hbbb and >30%')
grid = gridplot([[a, b], [c, d]], merge_tools=False)
save(grid)

for i in ["topology.pdb","tmp1.h5","tmp1.dcd","peptide_result.tsv", "filtered_peptide_result.tsv"]:
    if os.path.exists(i):
        os.remove(i)