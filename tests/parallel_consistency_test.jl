using Distributed
#@everywhere using LocalCellularStructure
using LocalCellularStructure
using Revise


testing_dir = homedir()*"/.julia/dev/LocalCellularStructure/tests/data/"


in_file_3D = testing_dir*"rand_pts_3D.csv"

path_out = testing_dir*"3D_test/testing_3D_non_periodic_custom_alpha"
#find_delaunay_network(in_file_3D,path_out, periodic=false, α = 1,edge_keep=true)
#compute_motifs(path_out,path_out)


#addprocs(4)
compute_flip(testing_dir*"3D_test/", path_out, restrict = 0, edge_keep = true, dimension = 3,thresh=0.5)

using LightGraphs
g1 = loadgraph(path_out*".lgz")
g2 = loadgraph(testing_dir*"3D_test/w_network.lgz")

if g1 == g2
    println("3D network test passed")
else
    error("Failed 3D network test")
end
