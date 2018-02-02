import particle_model as pm

from numpy import pi, cos, sin, linspace
import sys


def vel(x):
    return ( -2*pi*cos(pi*x[0]/2.0)*cos(2.0*pi*x[1]),
              pi*cos(pi*x[0])*sin(2.0*pi*x[1]), 0)

def pres(x):
    return 0

if len(sys.argv)<2:
    S='40'
else:
    S = sys.argv[1]

MESH = pm.IO.GmshMesh()
MESH.read('Structured%sx%s.msh'%(S,S))

for k, v in enumerate(linspace(0,0.1,11)):
    pm.IO.make_unstructured_grid(MESH,vel,pres,v,'gyre%sx%s_%d.vtu'%(S, S, k))
