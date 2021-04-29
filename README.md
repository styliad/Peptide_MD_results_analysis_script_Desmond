## Peptide Molecular Dynamics Simulation Results Analysis Script (Desmond)

The present script can be used to assist in results analysis for an MD Simulation of a peptide run in Desmond. 

# Usage

*"python peptide_md_analysis_script.py <trajectoryFile(clickme.dtr)> <in.cms File> <output.hmtl> <number_of_residues>"*

Take care of the relative position in directory

# Parameters to change before running the script
+ Plot size by modifying the **plot_width** and **plot_height parameters**
+ Depending on the ions and water models you use, the values in the peptide_atoms_indices variable have to be changed (e.g. T4P,NA,CL). This step cleans the trajectory from all het groups and solvent molecules
+ Desired bond interactions to calculate from GetContacts library, --itypes argument
+ Percentage of presence over the whole simulation (for Presence calculation)

# Definitions
Occurrence: The number of Occurrences of a single SS type related to a particular amino acid across all frames of the MD Simulation
dssp_data: The secondary structure (SS) type of a particular amino acid across in a particular frame of the MD

# Dependencies

+ vmd-python
+ bokeh
+ numpy
+ pandas
+ mdtraj
+ getcontacts
+ networkx

# References
+ Cinema colouring scheme: https://www.bioinformatics.nl/~berndb/aacolour.html
