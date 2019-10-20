using CSV

include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

Data_dir = "/Users/Dominic/Documents/2d Cells/LocalCellularStructure/"

Positions = [[0.5*(1+sin(k*0.153214)); 0.5(1+sin(1 + k*0.253214))] for k in 1:250]
weinberg2D_wrap(Positions, Data_dir*"Test_np_1",false)
weinberg2D_wrap(Positions, Data_dir*"Test_p_1",true)

w_tot = readin(Data_dir*"Test_p_",1)
append!(w_tot,readin(Data_dir*"Test_np_",1))
code_amalg = amalg2(w_tot)
network_save_str = Data_dir*"w_network_t"
compute_flip_graph(code_amalg,network_save_str)


d = calculate_distance_matrix(network_save_str,w_tot)
