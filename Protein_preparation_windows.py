#!/usr/bin/env python
# coding: utf-8

# ## Protein preparation

# # Step 1 Installation of dependecies

# Please install the following packages:
# 
# - PDBFixer
# - Biopython
# - MDAnalysis
# - RDKit
# - OpenMM (and OpenMMForceFields)
# - OpenBabel
# - Scrubber (package: "molscrub")
# - py3Dmol

# In[2]:


import subprocess
subprocess.run(['pip', 'install', 'rdkit', 'pdbfixer', 'openmm', 'mdanalysis', 'molscrub', 'py3dmol', 'biopython', 'openbabel'])


# In[16]:


#!pip install rdkit pdbfixer openmm mdanalysis molscrub py3dmol biopython


# In[3]:


# Install in Powershell
#!conda install -c conda-forge openbabel
get_ipython().system('pip install openbabel')


# # Step 2: Building Atomistic Ligand Model

# ## Step 2.1 Download PDB file

# In[17]:


import os
import requests

pdb_id = input("Enter PDB code: ") # The Protein ID we're looking at

# Start by making a directory for us to work in and stage our intermediate files
protein_directory = "molecular_docking/protein_files"
protein_filename = f"{pdb_id}.pdb"
protein_filepath = os.path.join(protein_directory, protein_filename)

# Actually make the directory, the exist_ok flag lets the command execute even if the folder already exists. It does NOT overwrite existing data.
os.makedirs(protein_directory, exist_ok=True)

print(protein_filepath)


# In[18]:


# Download the protein file
print(f"Downloading protein {pdb_id}...")
protein_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

# Send the request and save the returned JSON blob as a variable
protein_request = requests.get(protein_url)
protein_request.raise_for_status() # Check for errors


# In[19]:


# Save the actual text of the returned JSON blob as the PDB file we're used to
with open(protein_filepath, "w") as f:
    f.write(protein_request.text)
print(f"Saved protein to {protein_filepath}")


# In[20]:


# Display raw PDB file
from IPython.display import display, HTML

def render_text(text_blob):
  # Helper function for displaying text in Jupyter Notebooks in a scrollable object
  html = f"""
      <div style="height:400px; overflow:auto;">
          <pre>{text_blob}</pre>
      </div>
      """
  display(HTML(html))

render_text(protein_request.text)


# ## Step 2.2 Select one domain/chain of protein to work with

# In[21]:


# From now on we will work with only one domain/chain of the target protein
from Bio.PDB import PDBParser, Select, PDBIO

print("Selecting chain A if protein contains multiple chains...")
class ChainSelector(Select):
    def __init__(self, target_chain):
        self.target_chain = target_chain
    def accept_chain(self, chain):
        return chain.id == self.target_chain

# Load structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", f"{protein_directory}/{pdb_id}.pdb")

# Save each chain (monomer) as a separate PDB
#io = PDBIO()
#for model in structure:
#    for chain in model:
#        chain_id = chain.id
#        io.set_structure(structure)
#        io.save(f"monomer_{chain_id}.pdb", ChainSelector(chain_id))

# Save chain A as a separate PDB file
io = PDBIO()
io.set_structure(structure)
io.save(f"{protein_directory}/{pdb_id}_A.pdb", ChainSelector("A"))

print(f"Chain A of {pdb_id} was selected and saved as {pdb_id}_A.pdb!" )


# # Step 3 Preparation of receptor PDB structure

# ## Step 3.1 Fixing the PDB structure

# In[22]:


# Load the PDB into the PDBFixer class

from pdbfixer import PDBFixer

fixer = PDBFixer(filename=f"{protein_directory}/{pdb_id}_A.pdb")

print("Starting PDBFixer")
print("Fixing protein...")
# Fixing the structure at pH 7.4
fixer.findMissingResidues()
fixer.missingResidues
fixer.findNonstandardResidues()
print(fixer.nonstandardResidues)
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
print(fixer.missingAtoms)
print(fixer.missingTerminals)
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)

print("Fixing protein complete!")


# In[10]:


# Adding missing heavy atoms to PDB structure
from openmm.app import PDBFile

print("Adding missing heavy atoms to PDB structure...")
with open(f"{protein_directory}/{pdb_id}_A_fix_heavy.pdb", 'w') as f:
    
  # Toplology, Positions, file stream, and keep chain ID's
  PDBFile.writeFile(fixer.topology, fixer.positions, f, True)
print("Missing heavy atoms added!")
print(f"Structure was saved as {pdb_id}_A_fix_heavy.pdb")


# ## Step 3.2 Simple energy minimization

# In[11]:


from openmm.app import ForceField
from Bio.PDB import PDBParser, PDBIO

# Load and fix the structure
fixer = PDBFixer(filename=f"{protein_directory}/{pdb_id}_A_fix_heavy.pdb")
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.4)

# Load a force field (e.g., Amber)
forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')

# Create OpenMM system for minimization
#system = forcefield.createSystem(fixer.topology, ignoreExternalBonds=True)
system = forcefield.createSystem(fixer.topology, ignoreExternalBonds=False)

# Show our forces as OpenMM understands them
system.getForces()
print("Force field loaded!")


# In[12]:


# Loop through residues and print residue numbers

print(f"Loop through all residues of {pdb_id}_A")
# Load structure
print("Loading structure...")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", f"{protein_directory}/{pdb_id}_A_fix_heavy.pdb")
print("Structure loaded!")

for model in structure:
    for chain in model:
        for residue in chain:
            res_id = residue.get_id()
            res_num = res_id[1]  # residue number
            res_name = residue.get_resname()
            print(f"Chain {chain.id}, Residue {res_name} {res_num}")


# In[13]:


from openmm import VerletIntegrator
from openmm.app import Simulation
import openmm.unit as unit

# Use a generic VerletIntegrator
integrator = VerletIntegrator(0.001 * unit.picoseconds)

# Create simulation for minimization
# Optional, if you have access to a CUDA GPU, comment out the next line and uncomment the one after it
platform = None
# platform = Platform.getPlatformByName('CUDA')

# Define the OpenMM Simulation object, which serves as a convenince wrapper for the OpenMM Context and Reporter objects
simulation = Simulation(fixer.topology, system, integrator, platform)

# Set the position of our atoms
simulation.context.setPositions(fixer.positions)

# Minimize energy
print('Minimizing energy...')
simulation.minimizeEnergy()

# Get minimized positions. We have to copy the positions of the compiled object back into Python's memory
minimized_positions = simulation.context.getState(getPositions=True).getPositions()

# Write minimized structure to a PDB file
with open(f"{protein_directory}/{pdb_id}_A_fixed.pdb", 'w') as output:
    PDBFile.writeFile(fixer.topology, minimized_positions, output)

print(f'Minimization complete. Minimized structure saved to {protein_directory}/{pdb_id}_A_fixed_.pdb')


# ## Step 3.3 Adding Partial Charge Information to the Receptor for Docking

# In[14]:


print("Generating pdbqt file...")
# Invoke OpenBabel's CLI from Python. Can also use subprocess as its safer, but os.system works fine here.
receptor_pdbqt_path = f"{protein_directory}/{pdb_id}_A.pdbqt"
receptor_fixed_path = f"{protein_directory}/{pdb_id}_A_fixed.pdb"

# Generate the PDBQT file.
# We could have "--partialcharge <method>" as a flag if we wanted to compute the partial charges, but this will just assume they are all "0"
# The "-xh" flag preserves the hydrogens we worked so hard to get.
os.system(f"obabel -ipdb {receptor_fixed_path} -opdbqt -O {receptor_pdbqt_path}")
print(f"{pdb_id}_A.pdbqt has been generated and saved!")
print(f"{pdb_id}_A.pdbqt is ready for docking!")
# If you get a status code "2" here, rerun it. You want status code 0


# # Visualization with py3Dmol

# In[15]:


# view raw PDB
import py3Dmol

v = py3Dmol.view()
v.addModel(open(f"{protein_directory}/{pdb_id}_A_fixed.pdb").read())
v.setStyle({'chain':'A'}, {'cartoon': {'color': '#0e9674'}})
v.setStyle({'chain':'B'}, {'cartoon': {'color': '#c46225'}})
v.zoomTo({'model':0})
v.rotate(90, "z")
v.rotate(-25, "y")


# In[ ]:


# Voila! The protein is ready for docking!

