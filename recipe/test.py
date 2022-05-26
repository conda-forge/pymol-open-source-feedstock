import tempfile
from os.path import isfile, join

import numpy as np
import requests
from pymol2 import SingletonPyMOL

PDB_ID = "1l2y"

with tempfile.TemporaryDirectory() as td:
    CIF_PATH = join(td, PDB_ID + ".cif")
    PNG_PATH = join(td, "pymol_test.png")
    GRO_PATH = join(td, "ala_tripeptide.gro")
    G96_PATH = join(td, "ala_tripeptide.g96")

    r = requests.get(f"https://files.rcsb.org/download/{PDB_ID}.cif")
    with open(CIF_PATH, "w") as file:
        file.write(r.text)

    g96data = """TITLE
Unterminated alanine tripeptide
END
POSITION
    1 ALA   N          1    2.158437490    2.297964334    1.189144969
    1 ALA   CA         2    2.226037264    2.427364349    1.189144969
    1 ALA   C          3    2.376037359    2.409964323    1.189144969
    1 ALA   O          4    2.429137468    2.298264265    1.188044906
    1 ALA   CB         5    2.175237417    2.506564379    1.310844898
    1 ALA   H          6    2.213037491    2.204764366    1.189144969
    1 ALA   HA         7    2.199237347    2.481264353    1.096344948
    1 ALA   1HB        8    2.065637350    2.521564245    1.307244897
    1 ALA   2HB        9    2.197637320    2.455164194    1.406344891
    1 ALA   3HB       10    2.220837355    2.607064247    1.316545010
    2 ALA   N         11    2.451123714    2.514853001    1.190328956
    2 ALA   CA        12    2.596143007    2.498014450    1.190326452
    2 ALA   C         13    2.666078091    2.631841421    1.191806197
    2 ALA   O         14    2.603549242    2.738533974    1.193944693
    2 ALA   CB        15    2.633230448    2.412751913    1.067774415
    2 ALA   H         16    2.410219669    2.607194662    1.191243529
    2 ALA   HA        17    2.625542641    2.444583654    1.282609224
    2 ALA   1HB       18    2.584017992    2.313645840    1.070281863
    2 ALA   2HB       19    2.603428125    2.461166620    0.972738743
    2 ALA   3HB       20    2.741987944    2.394025326    1.062008619
    3 ALA   N         21    2.795035839    2.634985447    1.190808892
    3 ALA   CA        22    2.862635136    2.764377832    1.192239761
    3 ALA   C         23    3.012629986    2.746996880    1.190758109
    3 ALA   O         24    3.065718174    2.635332823    1.187519908
    3 ALA   CB        25    2.812823772    2.841701984    1.315544367
    3 ALA   H         26    2.848410606    2.549258471    1.189065218
    3 ALA   HA        27    2.835084438    2.819687605    1.100495100
    3 ALA   1HB       28    2.703198433    2.856742859    1.313062429
    3 ALA   2HB       29    2.835996628    2.788850307    1.410062551
    3 ALA   3HB       30    2.858469486    2.942108154    1.322410703
END
BOX
    3.413626671    3.413626671    2.413798571    0.000000000    0.000000000    0.000000000    0.000000000    1.706813335    1.706813335
END
"""
    with open(G96_PATH, "wt") as f:
        f.write(g96data)

    grodata = """Unterminated alanine tripeptide
   30
    1ALA      N    1   2.158   2.298   1.189
    1ALA     CA    2   2.226   2.427   1.189
    1ALA      C    3   2.376   2.410   1.189
    1ALA      O    4   2.429   2.298   1.188
    1ALA     CB    5   2.175   2.507   1.311
    1ALA      H    6   2.213   2.205   1.189
    1ALA     HA    7   2.199   2.481   1.096
    1ALA    1HB    8   2.066   2.522   1.307
    1ALA    2HB    9   2.198   2.455   1.406
    1ALA    3HB   10   2.221   2.607   1.317
    2ALA      N   11   2.451   2.515   1.190
    2ALA     CA   12   2.596   2.498   1.190
    2ALA      C   13   2.666   2.632   1.192
    2ALA      O   14   2.604   2.739   1.194
    2ALA     CB   15   2.633   2.413   1.068
    2ALA      H   16   2.410   2.607   1.191
    2ALA     HA   17   2.626   2.445   1.283
    2ALA    1HB   18   2.584   2.314   1.070
    2ALA    2HB   19   2.603   2.461   0.973
    2ALA    3HB   20   2.742   2.394   1.062
    3ALA      N   21   2.795   2.635   1.191
    3ALA     CA   22   2.863   2.764   1.192
    3ALA      C   23   3.013   2.747   1.191
    3ALA      O   24   3.066   2.635   1.188
    3ALA     CB   25   2.813   2.842   1.316
    3ALA      H   26   2.848   2.549   1.189
    3ALA     HA   27   2.835   2.820   1.100
    3ALA    1HB   28   2.703   2.857   1.313
    3ALA    2HB   29   2.836   2.789   1.410
    3ALA    3HB   30   2.858   2.942   1.322
   3.41363   3.41363   2.41380   0.00000   0.00000   0.00000   0.00000   1.70681   1.70681
"""

    with open(GRO_PATH, "wt") as f:
        f.write(grodata)

    pymol = SingletonPyMOL()
    pymol.start()
    cmd = pymol.cmd

    cmd.load(CIF_PATH)

    # Test rendering
    cmd.png(PNG_PATH, ray=1)
    assert isfile(PNG_PATH)

    # Test correct integration of NumPy
    assert cmd.get_coordset(PDB_ID).shape == (304, 3)

    # Test if we can load a .gro file
    cmd.load(GRO_PATH, 'grodata')
    assert cmd.get_coordset('grodata').shape == (30, 3)

    # Test if we can load a .g96 file
    cmd.load(G96_PATH, 'g96data')
    assert cmd.get_coordset('g96data').shape == (30, 3)

    # Check if the data from .gro and .g96 are the same
    diff = np.abs(cmd.get_coordset('g96data') - cmd.get_coordset('grodata')).max()
    assert diff < 0.008, diff  # difference should only come from rounding error
