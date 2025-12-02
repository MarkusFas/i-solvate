from glob import glob
import os
import argparse
from run_ase import run_dyn
import random
from ase.io import read


def main():
    parser = argparse.ArgumentParser(
            description="Run NPT or NVT simulations on XYZ files."
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
        default='npt',
        help="Ensemble to use: nvt or npt"
    )
    parser.add_argument(
        "--time",
        type=float,
        default=10.0,
        help="Total simulation time in ps (default: 10 ps)."
    )


    args = parser.parse_args()

    root = args.root
    outdir = args.output
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    model_path = args.model
    ensemble = args.ensemble
    total_time = args.time

    xyz_files = sorted(glob(os.path.join(root, "*.xyz")))[:]
    random.shuffle(xyz_files)
    os.makedirs(outdir, exist_ok=True)

    if not xyz_files:
        raise ValueError(f"No .xyz files found in {args.root}")

    temperature = [random.randint(300, 800) for _ in range(len(xyz_files))] # K
    pressure = [random.randint(0, 10) * 1e4 for _ in range(len(xyz_files))] # 0 to 1 GPa
    for i, f in enumerate(xyz_files[:10]):
        print(f"File: {f}, Temperature: {temperature[i]} K, Pressure: {pressure[i]/1e5} GPa")
        atoms = read(f)
        base = os.path.basename(f).replace(".xyz", "")
        out_file = os.path.join(outdir, f"{base}-{ensemble}.xyz")
        run_dyn(atoms, model_path, out_file, ensemble, total_time, temperature[i], pressure[i])
    

if __name__ == "__main__":
    main()
