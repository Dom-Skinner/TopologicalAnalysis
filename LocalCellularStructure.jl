module LocalCellularStructure

include("ReadWriteTools.jl")

include("VoronoiTools.jl")
using .VoronoiTools

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

export readin!, readin, amalg2, weinberg2D_wrap
end
