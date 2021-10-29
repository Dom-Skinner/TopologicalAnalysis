using Revise
using LocalCellularStructure
using CSV
using DataFrames
using Random

# This script is designed to ensure consistency between the old code and the
# refactored code

testing_dir = homedir()*"/Documents/TopologicalAnalysis/data/temp/"
# Create the testing data
Random.seed!(1234)
Positions = [randn(2) for k in 1:5000]
Positions_array = vcat(Positions'...)
in_file = homedir()*"/Documents/TopologicalAnalysis/data/temp/rand_pts.csv"
CSV.write(in_file,DataFrame(Positions_array))



# New code - test periodic
r = 2
path_out = homedir()*"/Documents/TopologicalAnalysis/data/temp/testing_periodic"
find_delaunay_network(in_file,path_out, periodic=true, α = 0,tol=0.6)
compute_motifs(path_out,path_out,r)

# New code - test non-periodic
path_out = homedir()*"/Documents/TopologicalAnalysis/data/temp/testing_non_periodic"
find_delaunay_network(in_file,path_out, periodic=false, α = 0,tol=0.6)
compute_motifs(path_out,path_out,r)

# New code - test non-periodic custom alpha
path_out = homedir()*"//Documents/TopologicalAnalysis/data/temp/testing_alpha"
find_delaunay_network(in_file,path_out, periodic=false, α = 0.5,tol=0.6)
compute_motifs(path_out,path_out,r)

## Now test against the old code (must restart Julia)
include(homedir()*"/Dropbox (MIT)/3DTopology/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure
using Random

r=2
Random.seed!(1234)
Positions = [randn(2) for k in 1:5000]

# This is straight up broken in the old code anyway...
#path_out_compare = "Documents/TopologicalAnalysis/data/temp/testing_periodic_old"
#weinberg2D_wrap(Positions, path_out_compare, periodic=true, r=r, α = 0)

path_out_compare = homedir()*"/Documents/TopologicalAnalysis/data/temp/testing_non_periodic_old"
weinberg2D_wrap(Positions, path_out_compare, periodic=false, r=r, α = 0)

path_out_compare = homedir()*"/Documents/TopologicalAnalysis/data/temp/testing_alpha_old"
weinberg2D_wrap(Positions, path_out_compare, periodic=false, r=r, α = 0.5)

## Compare the two
using LocalCellularStructure
testing_dir = homedir()*"/Documents/TopologicalAnalysis/data/temp/"

w1  =readin(testing_dir*"testing_non_periodic_old_avg.txt",0)
w2  =readin(testing_dir*"testing_non_periodic_avg.txt",0)

if w1 == w2
    println("Passed non-periodic test")
else
    error("Failed non-periodic test")
end

w3  =readin(testing_dir*"testing_alpha_old_avg.txt",0)
w4  =readin(testing_dir*"testing_alpha_avg.txt",0)

if w3 == w4
    println("Passed non-periodic test")
else
    error("Failed non-periodic test")
end
