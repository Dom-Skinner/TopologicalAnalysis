module FlipGraphTools

using CSV, DataFrames
using DelimitedFiles
using LightGraphs
using HDF5

include("../MotifLabelTools/MotifLabelTools.jl")
using .MotifLabelTools

include("NetworkDiscover.jl")
include("NetworkDiscover3D.jl")
include("../ReadWriteTools.jl")
include("../GraphConstruct.jl")


struct FlipGraph
           g::Graph
           motif_code::Dict
end

function save(save_str,fg::FlipGraph)
    E = collect(edges(fg.g))
    h5open(save_str, "w") do file
        write(file, "Type", "FlipGraph")
        write(file, "I", [e.src for e in E])
        write(file, "J", [e.dst for e in E])
        write(file, "codes", collect(keys(fg.motif_code)))
        write(file, "values", collect(values(fg.motif_code)))
    end
end

function load_flip_graph(save_str)

    I = h5read(save_str,"I")
    J = h5read(save_str,"J")
    # TODO: what is the vectorized way to do this? Not so important as this
    # is fast for unweighted graphs
    g = SimpleGraph(maximum(J))
    for i = 1:length(J)
        add_edge!(g,I[i],J[i])
    end

    codes = h5read(save_str,"codes")
    values = h5read(save_str,"values")
    code_to_idx = Dict(codes .=> values)
    return FlipGraph(g,code_to_idx)
end


function compute_flip(data_dir_in::String, path_out; restrict = 0,
                        edge_keep = false, dimension = 2, thresh = 1.5)

    str_arr = filter(x->occursin("avg.txt",x), readdir(data_dir_in))
    compute_flip(data_dir_in.*str_arr, path_out; restrict = restrict,
                    edge_keep = edge_keep, dimension = dimension, thresh = thresh)
end

function compute_flip(str_arr, path_out; restrict = 0, edge_keep = false,
                        dimension = 2,thresh=1.5)

    weight = amalg2(readin.(str_arr))
    keep = [w[2]>restrict for w in weight]
    if dimension == 2
        fg = compute_flip_graph(weight[keep])
    else
        fg = compute_flip_graph3D(weight[keep],edge_keep,thresh)
    end
    save(path_out,fg)
end


export compute_flip,compute_flip_graph,compute_flip_graph3D,load_flip_graph

end
