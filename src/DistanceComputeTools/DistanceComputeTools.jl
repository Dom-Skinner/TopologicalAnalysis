module DistanceComputeTools
using Distributed
using LightGraphs
using IterativeSolvers
using SparseArrays
using LinearAlgebra

include("../FlipGraphTools/FlipGraphTools.jl")
using .FlipGraphTools

include("../ReadWriteTools.jl")
include("Wasserstein.jl")



function calculate_distance_matrix(network_save_file, decode_save_file,file_out,
		str_arr; optimal_transport= true)
    # This function is a wrapper for all other functions in this file
    # from n dictionaries in and the path to the load_graph file it returns the
    # n by n distance matrix
	w_vec_in = readin.(str_arr)
	g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file, decode_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]

    d =  distance_mat(g,weight,optimal_transport)
	CSV.write(file_out,DataFrame(d,str_arr))

end


function calculate_distance_matrix(network_save_file,path_out,w_vec_in; optimal_transport= true)
    # This function is a wrapper for all other functions in this file
    # from n dictionaries in and the path to the load_graph file it returns the
    # n by n distance matrix

    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]
    d =  distance_mat(g,weight,optimal_transport)
	CSV.write(path_out,DataFrame(d,:auto))
end

function triangle_index(k)
	# Helps us to iterate through upper triangular part of matrix
	T = n ->Int(0.5*n*(n+1))
	kk = Int(ceil( sqrt(2*k + 0.25) - 0.5))
	return k - T(kk-1),kk + 1
end


function distance_mat(g,weight,optimal_transport)
	n_needed = Int(ceil(0.5*length(weight)*(length(weight)-1)))
    d = zeros(length(weight),length(weight))
	W = [weight[triangle_index(i)[1]] .- weight[triangle_index(i)[2]]  for i = 1:n_needed]
	rem_self_edges!(g)


	if optimal_transport
		f = x -> W_dist(g,x)
	else
		L = float.(laplacian_matrix(g))
		D = float.(incidence_matrix(g,oriented=true))
		f = x -> sum(abs.(D' * minres(L , x)))
	end

	d_flat = pmap(f,W)
	#d_flat = f.(W)
	for i  = 1:n_needed
		d[triangle_index(i)[1],triangle_index(i)[2]] = d_flat[i]
		d[triangle_index(i)[2],triangle_index(i)[1]] = d_flat[i]
	end
	return d
end

export  calculate_distance_matrix, geodesic_reg,load_w_graph
end
