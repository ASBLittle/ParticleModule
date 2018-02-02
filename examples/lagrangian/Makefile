### Makefile --- 

SIZE=40

Structured${SIZE}x${SIZE}.msh: Structured.geo	
	gmsh -2 -setnumber NS ${SIZE} Structured.geo -o Structured${SIZE}x${SIZE}.msh

input: Structured${SIZE}x${SIZE}.msh
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ python2 make_grids.py ${SIZE}

clean: 
	rm -rf gyre_*x*.vtu *.vtp

run: 	input
	PYTHONPATH=${PYTHONPATH}:${PWD}/../../ python2 pm.py ${SIZE}


### Makefile ends here
