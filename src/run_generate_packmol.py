import os 
import subprocess 

import numpy as np
import ase.io
from ase import data

from generate_packmol import generate_packmol_input

def call_packmol_generator(solutes, solvent, max_n_atoms):
    # solutes: list of tuple: iupac name, number of atoms
    root = 'all'
    
    packmol_file = "packmol.inp"
    
    n_atoms = 0
    n_ions = 0
    out_name = ""
    for solute in solutes:
        name, n = solute
        atoms = ase.io.read(os.path.join(root, name + '.xyz'))
        n_atoms += n * len(atoms)
        n_ions += n
        out_name += f"{n}{name}_" 
    
    solvent_name = solvent[0]
    out_name += solvent_name
    out_fname = out_name + ".xyz"
    # get allowed number of solvent molecules
    solvent_atoms = ase.io.read(os.path.join(root, solvent_name + '.xyz'))
    n_solvent = (max_n_atoms - n_atoms) // len(solvent_atoms)
    masses = [data.atomic_masses[atomic_number] for atomic_number in
              solvent_atoms.get_atomic_numbers()]
    molar_mass = np.sum(masses)
    density = solvent[1] # g/cm^3

    # calculate box lengths
    mass = molar_mass * (n_solvent + n_ions) / 6.022E23 # correct for ion
    #deplacement
    vol = mass / density
    # assuming a cubic box
    L = np.cbrt(vol) * 1E8 # from cm to A
    box_lengths = np.array([L,L,L])  # Box size in Angstroms
    
    generate_packmol_input(box_lengths, 
                           solutes,
                           solvent_name,
                           n_solvent,
                           out_fname, 
                           packmol_file,
                           )

    with open(packmol_file, "r") as f:
        run = subprocess.run(["packmol"], stdin=f, capture_output=True, text=True)
        print(run.stdout)
    # save as lammps.data
    structures = ase.io.read(out_fname)
    #structures.set_cell(box_lengths+1) # Add a small buffer to avoid periodic boundary issues
    #structures.set_pbc(True)
    #specorder = list(symbol_to_id.keys())[:102]
    #write(out_fname, structures, format='extxyz')
    #write("lammps_indexing.data", structures, specorder=[ielement for ielement in range(102)], format='lammps-data')
    #write("lammps_indexing.xyz", structures, format='extxyz')
    print(f"Packmol run completed. Output written to {out_fname}.")
