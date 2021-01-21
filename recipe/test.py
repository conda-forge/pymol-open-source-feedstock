import tempfile
from os.path import isfile, join
import requests
from pymol2 import SingletonPyMOL


PDB_ID = "1l2y"
CIF_PATH = join(tempfile.gettempdir(), PDB_ID+".cif")
PNG_PATH = join(tempfile.gettempdir(), "pymol_test.png")


r = requests.get(f"https://files.rcsb.org/download/{PDB_ID}.cif")
with open(CIF_PATH, "w") as file:
    file.write(r.text)


pymol = SingletonPyMOL()
pymol.start()
cmd = pymol.cmd

cmd.load(CIF_PATH)

# Test rendering
cmd.png(PNG_PATH, ray=1)
assert isfile(PNG_PATH)

# Test correct integration of NumPy
assert cmd.get_coordset(PDB_ID).shape == (304, 3)
