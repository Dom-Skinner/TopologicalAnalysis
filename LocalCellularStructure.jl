module LocalCellularStructure

include("ReadWriteTools.jl")
include("DistanceTools.jl")

include("VoronoiTools.jl")
using .VoronoiTools

include("WassersteinTools.jl")
using .WassersteinTools

include("Tools3D.jl")
using .Tools3D

function weinberg2D_wrap(Positions, Data_dir_str; periodic=false, r=2)
    if periodic
        Weinberg,S = weinberg2D(deepcopy(Positions),periodic,r)
    else
        Weinberg,S,idx = weinberg2D(deepcopy(Positions),periodic,r)
        Positions = Positions[idx]
    end
    write_total(Positions, Weinberg,S,Data_dir_str)
    write_avg(countmap(Weinberg),Data_dir_str)
end

function label_3D(Positions,Data_dir_str; r=1)
    p, simplices, neighbrs, edge_index = Delaunay_find(Positions)
    if r != 1
        error("Todo")
        # need to have better function for finding edge points if r > 1
    end
    not_edge = setdiff(1:size(p,1), edge_index)
    tvec_tot = Array{Int64}[]
    for i in not_edge
        k_nbhd = find_nbhd(simplices,i)
        push!(tvec_tot,topological_vec(k_nbhd,i))
    end
    write_total(Positions[not_edge], tvec_tot,Data_dir_str)
    write_avg(countmap(tvec_tot),Data_dir_str)
end

function compute_flip(Data_dir; restrict = 0, dimension = 2)
    str_arr = filter(x->occursin("avg.txt",x), readdir(Data_dir))
    weight = amalg2([readin(Data_dir*s,0) for s in str_arr])
    keep = [w[2]>restrict for w in weight]
    if dimension == 2
        compute_flip_graph(weight[keep],Data_dir*"w_network")
    else
        compute_flip_graph3D(weight[keep],Data_dir*"w_network")
    end
end

export
        # For reading/writing
        readin!, readin, amalg2, get_files_dir, write_avg, weight_in,get_weights_in_dir,

        # For computing topological types
        weinberg2D_wrap, label_3D,

        # For computing flip graph
        compute_flip_graph, compute_flip,

        # For calculating distances
        calculate_distance_matrix,calculate_distance_matrix_parallel,
        W_dist, fill_W_distance_mat,fill_JS_distance_mat,
        calculate_distance_matrix_lap,distance_mat_lap,

        fill_SN_distance_mat, # depricated

        # Misc.
        subsample_dist
end
