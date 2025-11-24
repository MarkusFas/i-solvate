import warnings
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in scalar add",
    category=RuntimeWarning,
)

import os
import argparse
from tqdm import tqdm
from ase.io import read, write
from ase.optimize import BFGS
from metatomic.torch.ase_calculator import MetatomicCalculator
from ase.filters  import UnitCellFilter


def main():
    parser = argparse.ArgumentParser(
        description="Relax all .xyz structures under a root directory using a Metatomic model."
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
        "--fmax",
        type=float,
        default=0.05,
        help="Force tolerance for relaxation. Default: 0.01 eV/Å."
    )

    args = parser.parse_args()

    root = args.root
    outdir = args.output
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    model_path = args.model
    fmax = args.fmax

    # Collect all xyz files first
    xyz_files = []
    for dirpath, dirnames, filenames in os.walk(root):
        for fname in filenames:
            if fname.endswith(".xyz"):
                xyz_files.append(os.path.join(dirpath, fname))

    print(f"Found {len(xyz_files)} .xyz files.\n")
    if len(xyz_files) == 0:
        return

    # Progress bar over all files
    for fpath in tqdm(xyz_files, desc="Relaxing structures", unit="file"):

        dirpath = os.path.dirname(fpath)
        fname = os.path.basename(fpath)

        try:
            atoms = read(fpath)
        except Exception as e:
            print(f"\nError reading {fpath}: {e}")
            continue

        # Attach calculator
        atoms.calc = MetatomicCalculator(model_path)
        
        # Relax
        logfile = os.path.join(outdir, fname.replace(".xyz", "-relax.out"))
        trajfile = os.path.join(outdir, fname.replace(".xyz", "-relax.traj"))

        uc = UnitCellFilter(atoms)
        dyn = BFGS(uc, logfile=logfile, trajectory=trajfile)

        #dyn = BFGS(atoms, logfile=logfile, trajectory=trajfile)
        
        try:
            dyn.run(fmax=fmax)
        except Exception as e:
            print(f"\nError relaxing {fpath}: {e}")
            continue

        # Write output
        base, ext = os.path.splitext(fname)
        outpath = os.path.join(outdir, f"{base}-relax{ext}")
        try:
            write(outpath, atoms)
        except Exception as e:
            print(f"\nError writing {outpath}: {e}")
            continue

    print("\nAll relaxations completed.")


if __name__ == "__main__":
    main()
