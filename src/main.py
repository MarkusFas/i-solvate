
from run_generate_packmol import call_packmol_generator 
# iupac name, number
solutes = [
        ("Li", 1),
        ("perchlorate", 1),
]

# iupac name, density
solvent = ('water', 1.0)

call_packmol_generator(solutes, solvent, 200)
