module PointCloudTools
using CSV, DataFrames
using HDF5

include("Voronoi.jl")
#include("Weinberg.jl")
include("../GraphConstruct.jl")

struct TopologicalNetwork
           simplices::Array
           not_edge::Array
           dim::Int
           periodic::Bool
           alpha::Number
           edge_keep::Bool
           tolerance::Number
           original_vertex_number::Int
end


function save(save_str,top::TopologicalNetwork)
    h5open(save_str, "w") do file
        write(file, "Type", "TopologicalNetwork")
        write(file, "simplices", top.simplices)
        write(file, "not_edge", top.not_edge)
        write(file, "dim", top.dim)
        write(file, "periodic", top.periodic)
        write(file, "alpha", top.alpha)
        write(file, "edge_keep", top.edge_keep)
        write(file, "tolerance", top.tolerance)
        write(file, "original_vertex_number", top.original_vertex_number)
    end
end

function load_topological_network(save_str)

    simplices = h5read(save_str,"simplices")
    not_edge = h5read(save_str,"not_edge")
    dim = h5read(save_str,"dim")
    periodic = h5read(save_str,"periodic")
    alpha = h5read(save_str,"alpha")
    edge_keep = h5read(save_str,"edge_keep")
    tolerance = h5read(save_str,"tolerance")
    original_vertex_number = h5read(save_str,"original_vertex_number")
    return TopologicalNetwork(simplices, not_edge, dim, periodic, alpha,
        edge_keep, tolerance, original_vertex_number)
end

function find_delaunay_network(path_to_csv_in, path_out; periodic=false, alpha = 0,
                                tol=0.6, edge_keep=false)

    dat_in = Matrix(CSV.read(path_to_csv_in,DataFrame))
    Positions = unique([dat_in[k,:] for k in 1:size(dat_in,1)])

    if length(Positions) != size(dat_in,1)
        println("Warning: Input file contains duplicate points.\n Recommended to rerun with duplicates removed")
    end

    top_net = find_delaunay_network_core(Positions, periodic, alpha, tol, edge_keep)
    save(path_out,top_net)

end

export find_delaunay_network, Delaunay_find, Voronoi_find_shape, periodic_extend!,
       graph_construct,vertex_num, shape_find, motif_size_find,
       find_delaunay_network_core, save,load_topological_network
end
