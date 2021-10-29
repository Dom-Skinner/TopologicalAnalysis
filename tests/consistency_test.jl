using Revise
using LocalCellularStructure
using CSV
using DataFrames
using Random

# This script is designed to ensure consistency between the old code and the
# refactored code

testing_dir = homedir()*"/.julia/dev/LocalCellularStructure/tests/data/"
in_file_2D = testing_dir*"rand_pts_2D.csv"
# New code - test periodic
r = 2
path_out = testing_dir*"testing_2D_periodic"
find_delaunay_network(in_file_2D,path_out, periodic=true, α = 0,tol=0.6)
compute_motifs(path_out,path_out,r)

# New code - test non-periodic
path_out = testing_dir*"testing_2D_non_periodic"
find_delaunay_network(in_file_2D,path_out, periodic=false, α = 0,tol=0.6)
compute_motifs(path_out,path_out,r)

# New code - test non-periodic custom alpha
path_out = testing_dir*"testing_2D_alpha"
find_delaunay_network(in_file_2D,path_out, periodic=false, α = 0.5,tol=0.6)
compute_motifs(path_out,path_out,r)


w1  =readin(testing_dir*"testing_2D_non_periodic_old_avg.txt",0)
w2  =readin(testing_dir*"testing_2D_non_periodic_avg.txt",0)

if w1 == w2
    println("Passed non-periodic test")
else
    error("Failed non-periodic test")
end

w3  =readin(testing_dir*"testing_2D_alpha_old_avg.txt",0)
w4  =readin(testing_dir*"testing_2D_alpha_avg.txt",0)

if w3 == w4
    println("Passed non-periodic test")
else
    error("Failed non-periodic test")
end

println("Passed 2D tests")

in_file_3D = testing_dir*"rand_pts_3D.csv"

path_out = testing_dir*"testing_3D_non_periodic_no_edge_keep"
find_delaunay_network(in_file_3D,path_out, periodic=false, α = 0,edge_keep=false)
compute_motifs(path_out,path_out)

path_out = testing_dir*"testing_3D_non_periodic_edge_keep"
find_delaunay_network(in_file_3D,path_out, periodic=false, α = 0,edge_keep=true)
compute_motifs(path_out,path_out)

path_out = testing_dir*"testing_3D_non_periodic_custom_alpha"
find_delaunay_network(in_file_3D,path_out, periodic=false, α = 1,edge_keep=true)
compute_motifs(path_out,path_out)



w1  =readin(testing_dir*"testing_3D_non_periodic_no_edge_keep_old_avg.txt",0)
w2  =readin(testing_dir*"testing_3D_non_periodic_no_edge_keep_avg.txt",0)

if w1 == w2
    println("Passed 3D no-edge keep test")
else
    error("Failed 3D no-edge keep test")
end


w1  =readin(testing_dir*"testing_3D_non_periodic_custom_alpha_old_avg.txt",0)
w2  =readin(testing_dir*"testing_3D_non_periodic_custom_alpha_avg.txt",0)

if w1 == w2
    println("Passed 3D alpha test")
else
    error("Failed 3D alpha test")
end



w1  =readin(testing_dir*"testing_3D_non_periodic_edge_keep_old_avg.txt",0)
w2  =readin(testing_dir*"testing_3D_non_periodic_edge_keep_avg.txt",0)

if w1 == w2
    println("Passed 3D edge keep test")
else
    error("Failed 3D edge keep test")
end
