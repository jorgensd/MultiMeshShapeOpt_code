## Overview

This folder contains source code for finding the optimal orientation of nine objects in Stokes flow, where drag is minimized.

### Files

- In the output folder, a **Paraview** for creating the figures from the second submission can be found. Note that **pvpython** is not included in the Docker-image.
- **scipysolver.py** contains the optimization problem used in the article. This folder uses classes from **StokesSolver.py**, which defines the MultiMesh Stokes problem for 9 objects in a channel with two inlets and one outlet.
- **ipoptsolver.py** contains a similar IPOPT implementation of the optimization problem. This is not used in the article.