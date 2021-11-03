module MotifLabelTools

using CSV, DataFrames
using DelimitedFiles
using LightGraphs

include("Labelling3D.jl")
include("Weinberg.jl")
include("../GraphConstruct.jl")
include("../ReadWriteTools.jl")


function compute_motifs(path_to_dir_in,path_out,r=-1)

    params_in = CSV.read(path_to_dir_in*".info", delim=", ",DataFrame)
    params_in = Dict(params_in.Parameter .=> params_in.value)
    if params_in["Graph type"] == "Delaunay"
        if params_in["Dimension"] == "2"
            if r < 0; r = 2; end
            idx, tvec  = weinberg2D(path_to_dir_in,params_in,r)

        elseif params_in["Dimension"] == "3"
            if r < 0; r = 1; end
            idx, tvec = simplicial_3D(path_to_dir_in,params_in)
        end
        write_total(idx, tvec,path_out)
        write_avg(countmap(tvec),path_out)
    end
end

export compute_motifs, motif_size_find, t_vec_to_simplex!,
       topological_vec_mem_save,topological_vec, weinberg_find!

end
