N=2

clean:
	rm -rf gyre_*.vtu *.vtp *.pvtp gyre_*.pvtu *.halo *.msh


input: clean

	gmsh -2 -bin Structured.geo
	${FLUIDITY_PATH}fldecomp -n ${N} Structured
	python make_grids.py ${N}	
