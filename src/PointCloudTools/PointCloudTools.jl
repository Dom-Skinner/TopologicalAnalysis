module PointCloudTools
using CSV, DataFrames

include("Voronoi.jl")
#include("Weinberg.jl")
include("../GraphConstruct.jl")

function find_delaunay_network(path_to_csv_in, path_out; periodic=false, alpha = 0,
                                tol=0.6, edge_keep=false)

    dat_in = Matrix(CSV.read(path_to_csv_in,DataFrame))
    Positions = unique([dat_in[k,:] for k in 1:size(dat_in,1)])
    dim = length(Positions[1])

    if dim == 2
        find_delaunay_network_2D(Positions, path_out, periodic, alpha, tol)
    elseif dim == 3
        find_delaunay_network_3D(Positions, path_out, periodic, alpha, tol, edge_keep)
    else
        error("Higher than 3D Delaunay not currently supported")
    end
end

export find_delaunay_network, Delaunay_find, Voronoi_find_shape, periodic_extend!,
       graph_construct,vertex_num, shape_find, motif_size_find
end
