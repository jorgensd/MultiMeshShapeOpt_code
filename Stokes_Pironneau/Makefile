all:
	mkdir -p meshes
	mkdir -p figures
	mkdir -p output
	python3 create_meshes.py 0.01
	python3 StokesSolver.py > output/Stokes_opt.txt
