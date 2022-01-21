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
