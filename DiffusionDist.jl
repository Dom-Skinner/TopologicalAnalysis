using LightGraphs
using SparseArrays
using LinearAlgebra


function rem_self_edges!(g)
    # Removes self loops
    for e in edges(g)
        if src(e) == dst(e)
            rem_edge!(g,e)
        end
    end
end

function distance_mat_lap(g,weight)
    rem_self_edges!(g)
    L = float.(laplacian_matrix(g))
    D = float.(incidence_matrix(g,oriented=true))

    d = zeros(length(weight),length(weight))
    for i = 1:size(d,1)
        for j=(i+1):size(d,2)
                    W = weight[i] .- weight[j]
                    J = D' * (L \ W)
                    d[i,j] = sum(abs.(J))
                    d[j,i] = d[i,j]
            end
            println(100*i/size(d,1))
    end
    return d # the distnace matrix
end

function calculate_distance_matrix_lap(network_save_file,w_vec_in)
    # Given n dictionaries in and the path to the load_graph file it returns the
    # n by n distance matrix under the diffusion approximation
    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]
    return distance_mat_lap(g,weight)
end

# No need for calculate_distance_matrix_parallel here as we can reuse old def.

#================================ Test example =================================
using Random,LightGraphs
g = barabasi_albert(1000, 3, seed=123)
is_connected(g)
Random.seed!(123)
ρ1 = rand(nv(g))
ρ2 = rand(nv(g))
ρ1 = ρ1/sum(ρ1)
ρ2 = ρ2/sum(ρ2)
approx = distance_mat_lap(g,[ρ1, ρ2])[1,2]
exact = W_dist(g,Δρ)
=## ============================================================================
