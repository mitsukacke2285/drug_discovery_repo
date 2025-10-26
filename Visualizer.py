#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#from google.colab import files

# Upload {pdb_id}_A.pdbqt, {pdb_id}_A_fixed.pdb and {pdb_id}_A.pdb files from local PC to your Colab VM
#files.upload("molecular_docking/protein_files")


# In[ ]:


# Upload {ligand_id}_ideal.sdf, {ligand_id}_corrected_pose.sdf and ligands_to_dock.sdf from local machine to google colab VM
#files.upload("molecular_docking/ligand_structures")


# In[1]:


import os
import requests
from util import visualize_poses
import py3Dmol

pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at
ligand_id = input("Enter ligand code: ") # The ID of the co-crystallized ligand

# Start by making a directory for us to work in and stage our intermediate files
protein_directory = "molecular_docking/protein_files"
protein_filename = f"{pdb_id}.pdb"
ligand_directory = "molecular_docking/ligand_structures"
ideal_ligand_filename = f"{ligand_id}_ideal.sdf" # Name of target ligand downloaded from RCSB PDB
docking_results_directory = "molecular_docking/docking_results"


# In[22]:


selection = input(f"Select molecule to view: {ligand_id}_fromPDB: a\
                                            {ligand_id}_ideal: b\
                                            {ligand_id}_corrected_pose: c\
                                            {ligand_id}_docked_{pdb_id}: d\
                                            multiple_ligands_docked: e\
                ")

if selection == 'a':
    print(f"{ligand_id}_fromPDB") 
    v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{ligand_directory}/{ligand_id}_fromPDB.pdb"
    )
    v.show()
elif selection == 'b':
    print(f"{ligand_id}_ideal") 
    v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{ligand_directory}/{ligand_id}_ideal.sdf"
    )
    v.show()
elif selection == 'c':
    print(f"{ligand_id}_corrected_pose") 
    v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{ligand_directory}/{ligand_id}_corrected_pose.sdf"
    )
    v.show()
elif selection == 'd':
    print(f"{ligand_id}_docked_{pdb_id}")
    v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{docking_results_directory}/{ligand_id}_docked_{pdb_id}.sdf",
    cognate_file=f"{ligand_directory}/{ligand_id}_corrected_pose.sdf"
    )
    v.show()
else:
    print("multiple_ligands_docked")
    v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{docking_results_directory}/multiple_ligands_docked.sdf",
    cognate_file=f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    animate=True,
    )  # Change to True to see an animation of all of the poses
    v.show()


# In[5]:


from util import visualize_poses

v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{ligand_directory}/{ligand_id}_ideal.sdf"
)
v.show()


# In[6]:


from util import visualize_poses

v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{ligand_directory}/{ligand_id}_corrected_pose.sdf"
)
v.show()


# In[7]:


import py3Dmol

# Compare the two ligand models
v = py3Dmol.view()
v.addModel(open(f"{ligand_directory}/{ligand_id}_fromPDB.pdb").read())
v.addModel(open(f"{ligand_directory}/{ligand_id}_corrected_pose.sdf").read())
#v.addModel(open(f"{ligand_directory}/{ligand_id}_ideal.sdf").read())
#v.setStyle()
v.setStyle({"model":0}, {'stick': {'color': '#0e9674'}})
v.setStyle({"model":1}, {'stick': {'color': '#c46225'}})
v.zoomTo({'model':0})


# In[10]:


v = visualize_poses(
    f"{protein_directory}/{pdb_id}_A_fixed.pdb",
    f"{docking_results_directory}/multiple_ligands_docked.sdf",
    cognate_file=f"{ligand_directory}/{ligand_id}_corrected_pose.sdf",
    animate=True,
)  # Change to True to see an animation of all of the poses
v.show()

