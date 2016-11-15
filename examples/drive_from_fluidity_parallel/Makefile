N = 2

clean:
	rm -rf *.vtu *.pvtu *.vtp *.pvtp *.stat drive_from_fluidity*[0-9] *.msh *.halo


input:
	gmsh -2 -bin Structured.geo
	${FLUIDITY_PATH}fldecomp -n ${N} Structured
