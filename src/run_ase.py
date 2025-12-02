#!/usr/bin/env python3
import argparse
import os
from glob import glob
import numpy as np

from ase.io import read, write
from ase.md.melchionna import MelchionnaNPT
from ase.md.langevinbaoab import LangevinBAOAB
from ase.md.nose_hoover_chain import IsotropicMTKNPT
from ase import units
from ase.md.bussi import Bussi
from ase.md.verlet import VelocityVerlet
from ase.md import MDLogger

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
    elif ensemble == 'npt-melchionna':
        dyn = MelchionnaNPT(
            atoms,
            timestep = dt,
            temperature_K = temperature,
            externalstress = pressure / 1602180,   # 1 bar
            mask = np.diag(np.ones(3,dtype=bool)),
            ttime= 100 * dt,
            pfactor= 100 * dt,
            )
    elif ensemble == 'npt-langevin':
        dyn = LangevinBAOAB(
            atoms,
            timestep = dt,
            temperature_K = temperature,
            externalstress = pressure * units.bar,   # 1 bar
            hydrostatic = 'True',
            #T_tau = 20 * units.fs,
            #P_tau = 20 * units.fs,
            rng = np.random.default_rng(int(time.time()))
        )
    elif ensemble == 'npt-mtk':
        dyn = IsotropicMTKNPT(
            atoms,
            timestep = dt,
            temperature_K = temperature,
            pressure_au = pressure / 1602180,   # convert bar to eV/Å³
            tdamp= 100 * dt,
            pdamp= 1000 * dt,
            )
        dyn.attach(MDLogger(dyn, atoms, out_file.replace('.xyz', '.log'), header=True, stress=True, peratom=False), interval=10)
    else:
        raise ValueError(f"Ensemble {ensemble} not recognized. Use 'nvt' or 'npt'.")

    trajectory = []
    atoms.info['energy'] = atoms.get_potential_energy()
    atoms.info['temperature'] = atoms.get_temperature() 
    atoms.info['stress'] = atoms.get_stress()
    atoms.arrays['forces'] = atoms.get_forces()
    trajectory.append(atoms.copy())

    for step in range(nsteps):
        # run a single simulation step
        dyn.run(1)
        if step%100==0:
            atoms.info['energy'] = atoms.get_potential_energy()
            atoms.info['temperature'] = atoms.get_temperature() 
            atoms.info['stress'] = atoms.get_stress()
            atoms.arrays['forces'] = atoms.get_forces()
            trajectory.append(atoms.copy())
    write(out_file,trajectory)
