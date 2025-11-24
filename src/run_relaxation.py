import os
from ase.io import read, write
from ase.optimize import BFGS
from metatomic.torch.ase_calculator import MetatomicCalculator

root = "solvated"   # change this

for dirpath, dirnames, filenames in os.walk(root):
    for fname in filenames:
        if fname.endswith(".xyz"):
            fpath = os.path.join(dirpath, fname)
            print(f"Relaxing {fpath}")

            # Read the structure (one structure per file)
            atoms = read(fpath)

            # Attach calculator
            atoms.calc = MetatomicCalculator("pet-mad-1.5.pt",device='cpu')

            # Relax
            dyn = BFGS(atoms, logfile=os.path.join(dirpath, fname.replace(".xyz", "-relax.out")))
            dyn.run(fmax=0.01)

            # Output filename
            base, ext = os.path.splitext(fname)
            outname = f"{base}-relax{ext}"
            outpath = os.path.join(dirpath, outname)

            # Write relaxed structure
            write(outpath, atoms)

            print(f" → Saved: {outpath}")

