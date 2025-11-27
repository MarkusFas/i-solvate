#!/usr/bin/env python3
import argparse
import os
from glob import glob

from ase.io import read, write
from ase.md.npt import NPT
from ase import units
from ase.md.bussi import Bussi
from ase.md.verlet import VelocityVerlet

from metatomic.torch.ase_calculator import MetatomicCalculator
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from tqdm import tqdm
import time


def run_dyn(atoms, model_path, out_file, ensemble, total_time):
    """
    Run NVT or NPT for 10 ps at T = 330 K and P = 1 bar.
    """

    atoms.calc = MetatomicCalculator(model_path)
    MaxwellBoltzmannDistribution(atoms, temperature_K=330)

    # 10 ps simulation: choose timestep = 0.5 fs
    dt = 0.5 * units.fs
    total_time *= 1000 * units.fs # convert ps to fs
    nsteps = int(total_time / dt)

    if ensemble == 'nvt':
        dyn = Bussi(
            atoms,
            timestep=dt,
            temperature_K = 330,
            taut = 50 * units.fs
            )
    elif ensemble == 'npt':
        dyn = NPT(
            atoms,
            timestep = dt,
            temperature_K = 330 ,
            externalstress = 1/1602180,   # 1 bar
            ttime= 50 * units.fs,
            pfactor= 50 * units.fs,
            )
    else:
        raise ValueError(f"Ensemble {ensemble} not recognized. Use 'nvt' or 'npt'.")

    trajectory = []
    atoms.info['energy'] = atoms.get_potential_energy()
    atoms.info['temperature'] = atoms.get_temperature() 
    atoms.info['stress'] = atoms.get_stress()
    atoms.arrays['forces'] = atoms.get_forces()
    trajectory.append(atoms.copy())

    for step in tqdm(range(nsteps), desc=f"Running {ensemble} dynamics"):
        # run a single simulation step
        dyn.run(1)
        if step%100==0:
            atoms.info['energy'] = atoms.get_potential_energy()
            atoms.info['temperature'] = atoms.get_temperature() 
            atoms.info['stress'] = atoms.get_stress()
            atoms.arrays['forces'] = atoms.get_forces()
            trajectory.append(atoms.copy())
    write(out_file,trajectory)

def main():
    parser = argparse.ArgumentParser(
            description="Run 50 ps NPT simulations on XYZ files."
            )
    parser.add_argument(
        "--root",
        required=True,
        help="Root directory containing .xyz files (searched recursively)."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--model",
        required=True,
        help="Path to the Metatomic .pt model file."
    )
    parser.add_argument(
        "--ensemble",
        required=True,
        default='nvt',
        help="Ensemble to use: nvt or npt"
    )
    parser.add_argument(
        "--time",
        type=float,
        default=10.0,
        help="Total simulation time in ps (default: 50 ps)."
    )

    args = parser.parse_args()

    root = args.root
    outdir = args.output
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    model_path = args.model
    ensemble = args.ensemble
    total_time = args.time

    xyz_files = sorted(glob(os.path.join(root, "*.xyz")))[:10]
    os.makedirs(outdir, exist_ok=True)

    if not xyz_files:
        raise ValueError(f"No .xyz files found in {args.root}")

    for f in xyz_files:
        atoms = read(f)

        base = os.path.basename(f).replace(".xyz", "")
        out_file = os.path.join(outdir, f"{base}-{ensemble}.xyz")

        run_dyn(atoms, model_path, out_file, ensemble, total_time)

if __name__ == "__main__":
    main()

