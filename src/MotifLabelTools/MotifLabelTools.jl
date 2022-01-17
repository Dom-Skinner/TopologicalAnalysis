module MotifLabelTools

using CSV, DataFrames
using DelimitedFiles
using LightGraphs

include("../PointCloudTools/PointCloudTools.jl")
using .PointCloudTools

include("Labelling3D.jl")
include("Weinberg.jl")
include("../GraphConstruct.jl")
include("../ReadWriteTools.jl")


function compute_motifs(path_to_file,path_out,r=-1)

    delaunay_in = load_topological_network(path_to_file)

    if delaunay_in.dim == 2
        if r < 0; r = 2; end
        idx, tvec  = weinberg2D(delaunay_in,r)

    elseif delaunay_in.dim == 3
        if r < 0; r = 1; end
        idx, tvec = simplicial_3D(delaunay_in)
    end
    write_total(idx, tvec,path_out)
    write_avg(countmap(tvec),path_out)

end

export compute_motifs, motif_size_find, t_vec_to_simplex!,
       topological_vec_mem_save,topological_vec, weinberg_find!,weinberg2D_core

end
