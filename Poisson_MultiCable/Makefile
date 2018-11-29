
create_meshes:
	mkdir -p meshes
	python refine_mesh.py 4

taylortest:
	mkdir -p meshes
	python3 refine_mesh.py 2
	python3	taylor_test.py > taylor_data.txt
	python3 refine_mesh.py 4
	python3	taylor_test.py >> taylor_data.txt
	python3 refine_mesh.py 6
	python3	taylor_test.py >> taylor_data.txt
	python3 refine_mesh.py 8
	python3	taylor_test.py >> taylor_data.txt


equilateral:
	mkdir -p meshes
	python3 refine_mesh.py 4
	python3 equilateral.py

timings:
	mkdir -p meshes
	python3 refine_mesh.py 4
