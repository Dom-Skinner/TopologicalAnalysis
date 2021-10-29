module LocalCellularStructure

include("./DistributionTools/DistributionTools.jl")
include("ReadWriteTools.jl")
include("./DistanceComputeTools/DistanceTools.jl")

include("./PointCloudTools/PointCloudTools.jl")
using .PointCloudTools

include("./DistanceComputeTools/WassersteinTools.jl")
using .DistanceComputeTools

include("./MotifLabelTools/MotifLabelTools.jl")
using .MotifLabelTools





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
        weinberg_find!, label_3D,find_delaunay_network,

        # For computing flip graph
        compute_flip_graph, compute_flip, compute_motifs

        # For calculating distances
        calculate_distance_matrix,calculate_distance_matrix_parallel,
        W_dist, fill_W_distance_mat,fill_JS_distance_mat,
        #calculate_distance_matrix_lap,distance_mat_lap,geodesic_reg,

    #    fill_SN_distance_mat, # depricated

        # Distribution tools
        tvec_dist,moments_find,find_dist_props,
        # Misc.
        subsample_dist, motif_size_find
end
