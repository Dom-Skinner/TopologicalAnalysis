# TopologicalAnalysis
This repository provides a julia package to perform topological analysis using the local structure of the Delaunay tessellation.

The inputs are csv files of 2D or 3D point clouds, the outputs are distances between those point clouds.

Demonstrations of the method are provided in /examples/

# Installation
First install julia and enter package mode by pressing `]`

Next, enter `add https://github.com/Dom-Skinner/TopologicalAnalysis` to install the package. Once installed, go back to normal mode by pressing backspace. 

To use enter `using LocalCellularStructure`. The first time you enter this it will install the required scipy python package through Conda.jl which may take some time.

To use the optimal transport distance or compute a curvature you must install a Gurobi license and verify that the Gurobi.jl package is correctly installed. Check the documentation at [github.com/jump-dev/Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) for instructions.

# Reference
If you make use of this code, consider citing the following publications

1. The 2D framework was introduced in [Skinner, Song, Jeckel, Jelli, Drescher, Dunkel, Phys. Rev. Lett. (2021)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.048101) 
2. The 3D framework is introduced in [Skinner, Jeckel, Martin, Drescher, Dunkel, Sci. Adv. (2023)](https://www.science.org/doi/full/10.1126/sciadv.adg1261)
