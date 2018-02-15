### Makefile --- 

SIZE=40
NP=1

Structured${SIZE}x${SIZE}.msh: Structured.geo	
	gmsh -2 -setnumber NS ${SIZE} Structured.geo -o Structured${SIZE}x${SIZE}.msh

gyre${SIZE}x${SIZE}_0.vtu: Structured${SIZE}x${SIZE}.msh
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ python2 make_grids.py ${SIZE}

clean: 
	rm -rf *.msh gyre*x*_[0-9]* gyre*x*_trajectories gyre*x*_trajectories.*

run: 	gyre${SIZE}x${SIZE}_0.vtu
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np ${NP} python2 pm.py ${SIZE}


### Makefile ends here
