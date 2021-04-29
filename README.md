# Peptide Molecular Dynamics Simulation Results Analysis Script (Desmond)

The present script can be used to assist in results analysis for MD Simulation of a peptide, run in Desmond. 

![Percentage of secondary structure per aminoacid](https://raw.githubusercontent.com/styliad/Peptide_MD_results_analysis_script_Desmond/master/images/Plot1.png?token=AH274VXNBCWFPBOGRL3WO6DARLCEW)
![Intrapeptide interactions Network Plot](https://raw.githubusercontent.com/styliad/Peptide_MD_results_analysis_script_Desmond/master/images/Plot2.png?token=AH274VQ26RBS2RQALDIM6Z3ARLCFM)

## Usage

*"python peptide_md_analysis_script.py <trajectoryFile(clickme.dtr)> <in.cms File> <output.hmtl> <number_of_residues>"*

## Parameters to change before running the script
+ Plot size by modifying the **plot_width** and **plot_height parameters**
+ Depending on the ions and water models you use, the values in the peptide_atoms_indices variable have to be changed (e.g. T4P,NA,CL). This step cleans the trajectory from all het groups and solvent molecules
+ Desired bond interactions to calculate from GetContacts library, --itypes argument
+ Percentage of presence over the whole simulation (for Presence calculation)

## Definitions
Occurrence: The number of occurrences of a single SS type related to a particular amino acid across all frames of the MD Simulation
dssp_data: The secondary structure (SS) type of a particular amino acid across in a particular frame of the MD

## Dependencies

+ vmd-python
+ bokeh
+ numpy
+ pandas
+ mdtraj
+ getcontacts
+ networkx

## References
+ Cinema colouring scheme: https://www.bioinformatics.nl/~berndb/aacolour.html
