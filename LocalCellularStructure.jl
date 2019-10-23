module LocalCellularStructure

include("ReadWriteTools.jl")

include("VoronoiTools.jl")
using .VoronoiTools

include("WassersteinTools.jl")
using .WassersteinTools

function weinberg2D_wrap(Positions, Data_dir_str,periodic=false)
    if periodic
        Weinberg,S = weinberg2D(deepcopy(Positions),periodic)
    else
        Weinberg,S,idx = weinberg2D(deepcopy(Positions),periodic)
        Positions = Positions[idx]
    end
    write_total(Positions, Weinberg,S,Data_dir_str)
    write_avg(countmap(Weinberg),Data_dir_str)
end

export readin!, readin, amalg2, weinberg2D_wrap, compute_flip_graph,
        calculate_distance_matrix,calculate_distance_matrix_parallel,W_dist
end
