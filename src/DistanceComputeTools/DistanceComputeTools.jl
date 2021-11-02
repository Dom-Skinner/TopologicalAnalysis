module DistanceComputeTools
using Distributed

#include("../PointCloudTools/PointCloudTools.jl")
#using .PointCloudTools

include("../ReadWriteTools.jl")

#include("../FlipGraphTools/FlipGraphTools.jl")
#using .FlipGraphTools

include("Wasserstein.jl")
include("DiffusionDist.jl")


function calculate_distance_matrix(network_save_file,w_vec_in; optimal_transport= true)
    # This function is a wrapper for all other functions in this file
    # from n dictionaries in and the path to the load_graph file it returns the
    # n by n distance matrix
    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]
    if optimal_transport
        return distance_mat(g,weight)
    else
        return distance_mat_lap(g,weight)
    end
end



function triangle_index(k)
	# Helps us to iterate through upper triangular part of matrix
	T = n ->Int(0.5*n*(n+1))
	kk = Int(ceil( sqrt(2*k + 0.25) - 0.5))
	return k - T(kk-1),kk + 1
end


function distance_mat(g,weight)
	n_needed = Int(ceil(0.5*length(weight)*(length(weight)-1)))
    d = zeros(length(weight),length(weight))
	W = [weight[triangle_index(i)[1]] .- weight[triangle_index(i)[2]]  for i = 1:n_needed]
	f = x -> distance_OT(g,x)
	d_flat = pmap(f,W)
	for i  = 1:n_needed
		d[triangle_index(i)[1],triangle_index(i)[2]] = d_flat[i]
		d[triangle_index(i)[2],triangle_index(i)[1]] = d_flat[i]
	end
	return d
end

export  calculate_distance_matrix, geodesic_reg
end
