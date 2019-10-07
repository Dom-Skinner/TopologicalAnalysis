module CommonTools

include("ReadWriteTools.jl")

include("VoronoiTools.jl")
using .VoronoiTools

function weinberg2D(Positions, Data_dir_str,periodic=false)
    Weinberg,S = weinberg2D(Positions,periodic)
    write_total(Positions, Weinberg,S,Data_dir_str)
    write_avg(Weinberg,Data_dir_str)
end

export readin!, readin, amalg2, weinberg2D
end
