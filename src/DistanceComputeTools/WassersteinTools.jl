module DistanceComputeTools


include("../PointCloudTools/PointCloudTools.jl")
using .PointCloudTools

include("../ReadWriteTools.jl")

#include("LocalCellularStructure.jl")
#using .LocalCellularStructure

include("../FlipGraphTools/NetworkDiscover.jl")
include("Wasserstein.jl")

export compute_flip_graph, calculate_distance_matrix,
        calculate_distance_matrix_parallel,W_dist,geodesic_reg
end
