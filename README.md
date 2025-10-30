# Semi-automated workflow for docking small molecules
This simple workflow for docking small molecule ligands into receptors is ideal for hit-to-lead or lead optimization. Its core modules are:
1. Protein preparation with PDBFixer, rdkit, MDAnalysis, Biopython
2. Ligand preparation with Scrubber (molscrub)
4. Docking with Gnina (based on Autodock Vina and Smina)
5. Visual analysis with py3Dmol
This workflow can be used by beginner and advanced users alike. The only input needed is the PDB id. The ligand(s) will be automatically identified after entering the PDB id in the ligand preparation module.
Scope and limitations:
1. rigid and flexible docking (later: free energy perturbation);
2. small molecule ligand docking;
3. docking into known binding site of co-crystallized ligand
4. docking into unknown binding site(s)

# Requirements/packages needed to be installed
Biopython, Gnina, MDAnalysis, Numpy, OpenBabel, OpenMM, Os, PDBFixer, Requests, Rdkit utils, Scrubber


