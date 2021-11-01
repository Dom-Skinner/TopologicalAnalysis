include(homedir()*"/Dropbox (MIT)/3DTopology/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

using CSV
using DataFrames
using Random

# This script is designed to ensure consistency between the old code and the
# refactored code

testing_dir = homedir()*"/.julia/dev/LocalCellularStructure/tests/data/"

# Create the testing data
Random.seed!(1234)
Positions = [randn(2) for k in 1:5000]
Positions_array = vcat(Positions'...)
CSV.write(testing_dir*"rand_pts_2D.csv",DataFrame(Positions_array))

r=2
weinberg2D_wrap(Positions, testing_dir*"testing_2D_non_periodic_old", periodic=false, r=r, α = 0)
weinberg2D_wrap(Positions, testing_dir*"testing_2D_alpha_old", periodic=false, r=r, α = 0.5)

# Create the testing data in 3D
Random.seed!(1234)
Positions = [randn(3) for k in 1:5000]
Positions_array = vcat(Positions'...)
CSV.write(testing_dir*"rand_pts_3D.csv",DataFrame(Positions_array))

label_3D(Positions, testing_dir*"testing_3D_non_periodic_no_edge_keep_old",
                edge_keep = false,α=0)
label_3D(Positions, testing_dir*"testing_3D_non_periodic_edge_keep_old",
                edge_keep = true,α=0)
label_3D(Positions, testing_dir*"testing_3D_non_periodic_custom_alpha_old",
                edge_keep = true,α=1)