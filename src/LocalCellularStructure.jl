module LocalCellularStructure

include("DistributionTools.jl")
include("ReadWriteTools.jl")
include("DistanceTools.jl")

include("VoronoiTools.jl")
using .VoronoiTools

include("WassersteinTools.jl")
using .WassersteinTools

include("Tools3D.jl")
using .Tools3D

function find_delaunay_network(path_to_csv_in,path_out; periodic=false, α = 0,tol=0.6)

    dat_in = Matrix(CSV.read(path_to_csv_in))
    Positions = unique([dat_in[k,:] for k in 1:size(dat_in,1)])
    dim = length(Positions[1])

    if dim == 2
        find_delaunay_network_2D(Positions, path_out, periodic, α, tol)
    elseif dim == 3
        find_delaunay_network_3D(Positions,path_out, periodic=periodic, α=α,tol=tol)
    else
        error("Higher than 3D Delanay not currently supported")
    end
end


function compute_motifs(path_to_dir_in,path_out,r)

    params_in = CSV.read(path_to_dir_in*".info", delim=", ")
    params_in = Dict(params_in.Parameter .=> params_in.value)
    if params_in["Graph type"] == "Delaunay"

        idx, Weinberg, S_tot = weinberg2D(path_to_dir_in,params_in,r)
        #write_total(idx, Weinberg,S,path_out)
        write_avg(countmap(Weinberg),path_out)
    end
end


function weinberg2D_wrap(Positions, Data_dir_str; periodic=false, r=2, α = 0)
    if periodic
        Weinberg,S = weinberg2D_old(deepcopy(Positions),periodic,r)
    else
        Weinberg,S,idx = weinberg2D_old(deepcopy(Positions),periodic,r,α=α)
        Positions = Positions[idx]
    end
    write_total(Positions, Weinberg,S,Data_dir_str)
    write_avg(countmap(Weinberg),Data_dir_str)
end

function label_3D(Positions,Data_dir_str; r=1, edge_keep = false,α = 0)
    p, simplices, neighbrs, edge_index,α_val = Delaunay_find(Positions, α = α)
    if r != 1
        error("Todo")
        # need to have better function for finding edge points if r > 1
    end

    if edge_keep
        not_edge = [i for i in 1:size(p,1)] # keep all indices
        simplices = simplices[α_val .== 1,:] # but only keep simplices that have alpha val < α
    else
        not_edge = setdiff(1:size(p,1), edge_index)
    end
    tvec_tot = Array{Int64}[]
    for i in 1:length(not_edge)

            k_nbhd = find_nbhd(simplices,not_edge[i])
            if length(k_nbhd) > 0
                push!(tvec_tot,topological_vec(k_nbhd,not_edge[i]))
            else
                not_edge[i] = -1
            end


    end
    not_edge = not_edge[not_edge.>0]
    write_total(Positions[not_edge], tvec_tot,Data_dir_str)
    write_avg(countmap(tvec_tot),Data_dir_str)
end

function compute_flip(Data_dir; restrict = 0, edge_keep = false, dimension = 2,thresh=1.5)
    str_arr = filter(x->occursin("avg.txt",x), readdir(Data_dir))
    weight = amalg2([readin(Data_dir*s,0) for s in str_arr])
    keep = [w[2]>restrict for w in weight]
    if dimension == 2
        compute_flip_graph(weight[keep],Data_dir*"w_network")
    else
        compute_flip_graph3D(weight[keep],Data_dir*"w_network",edge_keep,thresh)
    end
end

export
        # For reading/writing
        readin!, readin, amalg2, get_files_dir, write_avg, weight_in,get_weights_in_dir,

        # For computing topological types
        weinberg2D_wrap, weinberg_find!, label_3D,find_delaunay_network,

        # For computing flip graph
        compute_flip_graph, compute_flip, compute_motifs

        # For calculating distances
        calculate_distance_matrix,calculate_distance_matrix_parallel,
        W_dist, fill_W_distance_mat,fill_JS_distance_mat,
        calculate_distance_matrix_lap,distance_mat_lap,geodesic_reg,

    #    fill_SN_distance_mat, # depricated

        # Distribution tools
        tvec_dist,moments_find,find_dist_props,
        # Misc.
        subsample_dist, motif_size_find
end
