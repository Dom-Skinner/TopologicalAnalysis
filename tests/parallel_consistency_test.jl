using Distributed
#@everywhere using LocalCellularStructure
using Revise
using LocalCellularStructure
using LightGraphs
using Random
using CSV,DataFrames

testing_dir = homedir()*"/.julia/dev/LocalCellularStructure/tests/data/"


in_file_3D = testing_dir*"rand_pts_3D.csv"

path_out = testing_dir*"3D_test/testing_3D_non_periodic_custom_alpha"
#find_delaunay_network(in_file_3D,path_out, periodic=false, α = 1,edge_keep=true)
#compute_motifs(path_out,path_out)


#addprocs(4)
compute_flip(testing_dir*"3D_test/", path_out, restrict = 0, edge_keep = true, dimension = 3,thresh=0.5)


g1 = loadgraph(path_out*".lgz")
g2 = loadgraph(testing_dir*"3D_test/w_network.lgz")

if g1 == g2
    println("3D network test passed")
else
    error("Failed 3D network test")
end


in_file_2D = testing_dir*"rand_pts_2D.csv"

path_out = testing_dir*"2D_test/testing_2D_alpha"
compute_flip(testing_dir*"2D_test/", path_out, restrict = 0, dimension = 2)

g1 = loadgraph(path_out*".lgz")
g2 = loadgraph(testing_dir*"2D_test/w_network.lgz")

if g1 == g2
    println("2D network test passed")
else
    error("Failed 2D network test")
end



path_out = testing_dir*"2D_test/testing_2D_alpha"
network_save_file = testing_dir*"2D_test/w_network"
weight_original = readin(path_out*"_old_avg.txt",0)
Random.seed!(1234)
weight_arr = [Dict(keys(weight_original) .=> [rand(1:100) for k in keys(weight_original)]) for i = 1:3]


calculate_distance_matrix(network_save_file,path_out*"_OT_",weight_arr,optimal_transport=true)
res_old = Matrix(CSV.read(path_out*"_old_OT_distance.txt",DataFrame))
res = Matrix(CSV.read(path_out*"_OT_distance_matrix.txt",DataFrame))
if res == res_old
    println("OT distance test passed")
else
    println(maximum(res .- res_old))
    println("OT distance test maybe failed")
end

calculate_distance_matrix(network_save_file,path_out*"_Diff_",weight_arr,optimal_transport=false)
res_old = Matrix(CSV.read(path_out*"_old_Diff_distance.txt",DataFrame))
res = Matrix(CSV.read(path_out*"_Diff_distance_matrix.txt",DataFrame))
if res == res_old
    println("Diffusion distance test passed")
else
    error("Diffusion distance test failed")
end
