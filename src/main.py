
from run_generate_packmol import call_packmol_generator, smallest_stoichiometry
import pandas as pd
from tqdm import tqdm

MAX_N_ATOMS = 200

# Reading database ...
# write all combinations of cations, anions and solvents
cations_list = ['../structures/atomic_cations.csv',
                '../structures/molecular_cations.csv'
        ]
anions_list = ['../structures/atomic_anions.csv',
               '../structures/molecular_anions.csv'
        ]
cations = pd.concat([pd.read_csv(f, delimiter=';') for f in cations_list])
anions = pd.concat([pd.read_csv(f, delimiter=';') for f in anions_list])

solvents = pd.read_csv('../structures/solvents.csv', delimiter=';')

# Total number of combinations
print("Total combinations to generate:")
total_combinations = len(cations) * len(anions) * len(solvents)
print(total_combinations)

# Progress bar
pbar = tqdm(total=total_combinations, desc="Generating systems")

# Generating input 
for ication, cation in cations.iterrows():
    for ianion, anion in anions.iterrows():
        m, n = smallest_stoichiometry(cation["Charge"], anion["Charge"])
        solutes = [
            (cation["IUPAC_Name"], m),
            (anion["IUPAC_Name"], n),
        ]
        for isolvent, solvent in solvents.iterrows():
            solvent = (solvent["IUPAC_Name"], solvent["Density(g/cm3)"])
            try:
                call_packmol_generator(solutes, solvent, MAX_N_ATOMS)
            except Exception as e:
                print(f"Error generating system {solutes} in {solvent[0]}: {e}")
            finally:
                pbar.update(1)

pbar.close()
