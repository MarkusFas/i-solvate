import os 
import subprocess 

import numpy as np
import ase.io
from ase import data
from math import gcd
from generate_packmol import generate_packmol_input

def smallest_stoichiometry(c_charge, a_charge):
    # We want integers m, n > 0 so that m*c_charge + n*a_charge == 0
    # let m = abs(a_charge) / gcd(|c|,|a|) * k, n = abs(c_charge) / gcd(|c|,|a|) * k
    g = gcd(abs(c_charge), abs(a_charge))
    m = abs(a_charge) // g
    n = abs(c_charge) // g
    return m, n

def call_packmol_generator(solutes, solvent, max_n_atoms):
    # solutes: list of tuple: iupac name, number of atoms
    in_root = '../structures'
    out_root = 'solvated'
    if not os.path.exists(out_root):
        os.makedirs(out_root)
    packmol_file = "packmol.inp"
    
    n_atoms = 0
    n_ions = 0
    out_name = ""
    for solute in solutes:
        name, n = solute
        atoms = ase.io.read(os.path.join(in_root, name + '.xyz'))
        n_atoms += n * len(atoms)
        n_ions += n
        out_name += f"{n}{name}_" 
    
    solvent_name = solvent[0]
    out_name += solvent_name
    out_fname = os.path.join(out_root, out_name + ".xyz")
    if os.path.exists(out_fname):
        print(f"{out_fname} already exists. Skipping packmol generation.")
        return
    
    # get allowed number of solvent molecules
    solvent_atoms = ase.io.read(os.path.join(in_root, solvent_name + '.xyz'))
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
        # write stdout and stderr to files
        with open(out_fname[:-3] + 'out', 'w') as fout:
            fout.write(run.stdout)
    # save as lammps.data
    structures = ase.io.read(out_fname)
    structures.info.clear()

    structures.set_cell(box_lengths) # Add a small buffer to avoid periodic boundary issues
    structures.set_pbc(True)
    ase.io.write(out_fname, structures, format='extxyz')