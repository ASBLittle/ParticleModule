### Makefile --- 

SIZE=40
NP=8

Structured${SIZE}x${SIZE}.msh: Structured.geo	
	gmsh -2 -setnumber NS ${SIZE} Structured.geo -o Structured${SIZE}x${SIZE}.msh

gyre${SIZE}x${SIZE}_0.vtu: Structured${SIZE}x${SIZE}.msh
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ python2 make_grids.py ${SIZE}

clean: 
	rm -rf *.msh gyre*x*_[0-9]* gyre*x*_trajectories gyre*x*_trajectories.*

run: 	gyre${SIZE}x${SIZE}_0.vtu
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np ${NP} python2 pm.py 0.0 ${SIZE}
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np ${NP} python2 pm.py 0.01 ${SIZE}
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np ${NP} python2 pm.py 0.001 ${SIZE}

test:
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 1 python2 pm_test.py 0.01 ${SIZE} 1
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 2 python2 pm_test.py 0.01 ${SIZE} 2
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 3 python2 pm_test.py 0.01 ${SIZE} 3
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 4 python2 pm_test.py 0.01 ${SIZE} 4
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 5 python2 pm_test.py 0.01 ${SIZE} 5
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 6 python2 pm_test.py 0.01 ${SIZE} 6
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 7 python2 pm_test.py 0.01 ${SIZE} 7
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 8 python2 pm_test.py 0.01 ${SIZE} 8
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 9 python2 pm_test.py 0.01 ${SIZE} 9
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 10 python2 pm_test.py 0.01 ${SIZE} 10
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 11 python2 pm_test.py 0.01 ${SIZE} 11
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 12 python2 pm_test.py 0.01 ${SIZE} 12
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 13 python2 pm_test.py 0.01 ${SIZE} 13
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 14 python2 pm_test.py 0.01 ${SIZE} 14
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 15 python2 pm_test.py 0.01 ${SIZE} 15
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ mpirun -np 16 python2 pm_test.py 0.01 ${SIZE} 16
	


### Makefile ends here
