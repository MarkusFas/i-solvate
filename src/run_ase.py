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


def run_dyn(atoms, model_path, out_file, ensemble, total_time, temperature, pressure):
    """
    Run NVT or NPT simulation.
    """

    atoms.calc = MetatomicCalculator(model_path)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

    # choose timestep = 0.5 fs
    dt = 0.5 * units.fs
    total_time *= 1000 * units.fs # convert ps to fs
    nsteps = int(total_time / dt)

    if ensemble == 'nvt':
        dyn = Bussi(
            atoms,
            timestep=dt,
            temperature_K = temperature,
            taut = 200 * units.fs
            )
    elif ensemble == 'npt':
        dyn = NPT(
            atoms,
            timestep = dt,
            temperature_K = temperature,
            externalstress = pressure / 1602180,   # 1 bar
            ttime= 200 * units.fs,
            pfactor= 200 * units.fs,
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
        if step%10==0:
            atoms.info['energy'] = atoms.get_potential_energy()
            atoms.info['temperature'] = atoms.get_temperature() 
            atoms.info['stress'] = atoms.get_stress()
            atoms.arrays['forces'] = atoms.get_forces()
            trajectory.append(atoms.copy())
    write(out_file,trajectory)
