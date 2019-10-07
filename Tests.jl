using CSV

include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

Data_dir = "/Users/Dominic/Documents/2d Cells/LocalCellularStructure/"

Positions = [[0.5*(1+sin(k*0.153214)); 0.5(1+sin(1 + k*0.253214))] for k in 1:250]
weinberg2D_wrap(Positions, Data_dir*"Test",false)
weinberg2D_wrap(Positions, Data_dir*"Test",true)
