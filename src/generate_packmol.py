import ase
from ase.io import read, write
import numpy as np
import sys
import subprocess
import os 


def generate_packmol_input(
    box_lengths,
    solutes,
    solvent,
    n_solvent,
    output_xyz="packed.xyz",
    packmol_inp="packmol.inp",
    tolerance=2.0,
):
    root = '../structures/'
    Lx, Ly, Lz = box_lengths

    lines = []
    lines.append(f"tolerance {tolerance}")
    lines.append("filetype xyz")
    lines.append(f"output {output_xyz}")
    lines.append(f"pbc {Lx} {Ly} {Lz}")

    # Random placement for all solutes
    for sol in solutes:
        name, n = sol
        fname = os.path.join(root, name + '.xyz')
        lines.append(f"\nstructure {fname}")
        lines.append(f"  number {n}")
        lines.append(f"  inside box 0. 0. 0. {Lx} {Ly} {Lz}")
        lines.append("end structure")

    # Solvent
    fname = os.path.join(root, solvent + '.xyz')
    lines.append(f"\nstructure {fname}")
    lines.append(f"  number {n_solvent}")
    lines.append(f"  inside box 0. 0. 0. {Lx} {Ly} {Lz}")
    lines.append("end structure\n")

    with open(packmol_inp, "w") as f:
        f.write("\n".join(lines))

    #print(f"✓ packmol input written to {packmol_inp}")


def center_molecule(mol_file, box_lengths):
    mol = read(mol_file)
    # Center the molecule in the middle of the box
    mol.center(vacuum=0.0)
    # Shift molecule to the exact center of the box
    mol.translate(box_lengths / 2 - mol.get_center_of_mass())
    write(mol_file, mol, format='xyz')

def get_num(vol):
    mass = rho*vol
    num = (mass/molmass) * 6.022E23
    num -= 52 # correcting for the volume taken up by the dota molecules
    return round(num)


