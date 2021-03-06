import os
import sys
import math

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

def test_rotate(axis, theta):
    lat = create.Lattice()
    lat.set_lattice_miller(5, [1,1,0], math.radians(45))
    mol = create.Molecule(["Ag", "Ag"], numpy.array([[0,0,0],[0,0,1]],dtype=float))
    mol.rotate(axis, theta)

    cry = create.Crystal(mol, lat)

    k = cry.get_crystal()

    wri = write.xyz("out.xyz")
    wri.write(k[0], k[1])


test_rotate([0, 1, 0], math.radians(15))
