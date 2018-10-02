all:
	mkdir -p meshes
	mkdir -p figures
	mkdir -p output
	python create_meshes.py 0.01	
	python Stokes_Optimization.py > output/Stokes_opt.txt

wedge:
	mkdir -p meshes
	mkdir -p figures
	mkdir -p output
	python wedge_meshes.py 0.01
	python Stokes_Optimization.py True > output/Stokes_wedge_opt.txt
