module FlipGraphTools

using CSV, DataFrames
using DelimitedFiles
using LightGraphs

include("../MotifLabelTools/MotifLabelTools.jl")
using .MotifLabelTools

include("NetworkDiscover.jl")
include("NetworkDiscover3D.jl")
include("../ReadWriteTools.jl")
include("../GraphConstruct.jl")

function compute_flip(data_dir_in, path_out; restrict = 0, edge_keep = false, dimension = 2,thresh=1.5)
    str_arr = filter(x->occursin("avg.txt",x), readdir(data_dir_in))
    weight = amalg2([readin(data_dir_in*s,0) for s in str_arr])
    keep = [w[2]>restrict for w in weight]
    if dimension == 2
        compute_flip_graph(weight[keep],path_out)
    else
        compute_flip_graph3D(weight[keep],path_out,edge_keep,thresh)
    end
end

export compute_flip

end
