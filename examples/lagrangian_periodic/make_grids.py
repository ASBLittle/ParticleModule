""" Make example data for Lagrangian advection with periodic boundaries."""
import sys

from numpy import pi, linspace
import particle_model as pm

def vel(pos):
    """Define velocity."""
    del pos
    return (2.0*pi, 0, 0)

def pres(pos):
    """Define pressure."""
    del pos
    return 0

if len(sys.argv) < 2:
    S = '40'
else:
    S = sys.argv[1]

MESH = pm.IO.GmshMesh()
MESH.read('Structured%sx%s.msh'%(S, S))

for level, time in enumerate(linspace(0, 0.1, 11)):
    pm.IO.make_unstructured_grid(MESH, vel, pres, time, 'gyre%sx%s_%d.vtu'%(S, S, level))
