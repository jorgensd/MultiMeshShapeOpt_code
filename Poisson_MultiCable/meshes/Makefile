cable.msh:
	# gmsh -2 cable_halo.geo
	# dolfin-convert cable_halo.msh cable_halo.xml
	gmsh -2 cable.geo
	dolfin-convert cable.msh cable.xml
	# gmsh -2 inner_cable.geo
	# dolfin-convert inner_cable.msh inner_cable.xml
	gmsh -2 inner_cable_halo.geo
	dolfin-convert inner_cable_halo.msh inner_cable_halo.xml
