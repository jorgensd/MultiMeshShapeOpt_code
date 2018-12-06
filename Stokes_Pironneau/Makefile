all:
	mkdir -p meshes
	mkdir -p figures
	mkdir -p output
	python create_meshes.py 0.01
	python StokesSolver.py > output/Stokes_opt.txt
