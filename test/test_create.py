import os
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = script_dir + "/../"
sys.path.append(module_dir)

from simpy import create
from simpy import write
import numpy

def test_type(ltype):
    lat = create.Lattice(lattice_system=ltype)
    mol = create.Molecule(["Ag"], numpy.array([[0,0,0]]))

    cry = create.Crystal(mol, lat)

    k = cry.get_crystal()

    wri = write.xyz("out.xyz")
    wri.write(k[0], k[1])


test_type('fcc')
