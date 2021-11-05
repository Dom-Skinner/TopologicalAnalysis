using LightGraphs
using Base.Threads
using DataFrames,CSV
using Distributed


function circ_insert!(M,pair_set,val)
    idx1 = findfirst(x-> x ∈ pair_set , M)
    idx2 = findlast(x-> x ∈ pair_set , M)

    if idx2 - idx1 == 1
        insert!(M,idx2,val)
    else
        push!(M,val)
    end
end

function weinberg_flip(g,cent_node,order_mat,d,s,v1,v2,r)

    # First we perform the flipping procedure on vertices d,s,v1,v2, updating
    # the graph and the order matrices
    g2 = deepcopy(g)
    rem_edge!(g2,d,s)
    add_edge!(g2,v1,v2)
    order_copy = deepcopy(order_mat)
    order_copy[d] = order_mat[d][order_mat[d] .!= s]
    order_copy[s] = order_mat[s][order_mat[s] .!= d]
    circ_insert!(order_copy[v1],[d;s],v2)
    circ_insert!(order_copy[v2],[d;s],v1)

    # Now we use the induced induced_subgraph to compute the weinberg
    code_tot = Vector{Array{Int64}}(undef,1)
    S_tot = Array{Int64}(undef,1)
    g_ego, vmap = induced_subgraph(g2,neighborhood(g2,cent_node,r))
    vmap_inv = Dict(vmap[k] => k for k in 1:length(vmap))
    order_local = order_copy[vmap]
    for i = 1:length(vmap)
        total_order = map.(x -> get(vmap_inv, x, -1), order_copy[vmap[i]])
        order_local[i] = total_order[total_order .> 0]
    end
    if (minimum([length(neighbors(g_ego,k)) for k = 1:nv(g_ego)]) < 3) || (has_edge(g,v1,v2))
        # don't flip!!
        return [-1]
    end
    weinberg_find!(code_tot,S_tot,1,g_ego,order_local,vmap_inv[cent_node])
    if S_tot[1] == -1
        # This is bad
        savegraph( "debug2.lgz", g)
        error("..")
    end
    return code_tot[1]
end

function return_nbh_vertex(order_arr,v)
    next_ind = (i,m) -> (i+1)*(i<m) +(i == m)
    prev_ind = (i,m) -> (i-1)*(i>1) +m*(i == 1)

    d1 = findfirst(x->x==v, order_arr)
    v1 = order_arr[ next_ind(d1,length(order_arr))]
    v2 = order_arr[ prev_ind(d1,length(order_arr))]
    return v1,v2
end


function w_vec_neighbors(w_in,r)
    # This function takes a weinberg vector and creates finds all other weinberg
    # vectors that are 1 flip away while also increasing the disntance from central node
    w_neighbors = []
    g = SimpleGraph(maximum(w_in))
    for i = 1:(length(w_in)-1)
        add_edge!(g,w_in[i],w_in[i+1])
    end

    # need to do the embedding to find the order mat
    x,y,fixed_vecs = tutte_embedding(g)
    order_mat = order_mat_find(g,x,y)

    cent_node = 1
    dist_to_cent = (dijkstra_shortest_paths(g,cent_node)).dists
    for e in edges(g)
        d = dst(e)
        s = src(e)

        v11,v21 = return_nbh_vertex(order_mat[s],d)
        v22,v12 = return_nbh_vertex(order_mat[d],s)

        # We have picked something that looks flippable. We test if it is and
        # then call weinberg_flip to find the w-vector after the flip. NB!!
        # weinberg_flip might decide to not flip and will return [-1]
        if (v11 == v12) && (v21 == v22)
            #println("A Flippable triangle")
            if (dist_to_cent[d] == dist_to_cent[s]) &&
                (abs(dist_to_cent[v11]-dist_to_cent[v21]) == 2)
                #println("Don't flip!")
            else
                push!(w_neighbors, weinberg_flip(g,cent_node,order_mat,d,s,v11,v21,r))
            end
        end
    end
    return w_neighbors
end

function flip_core_parallel(x,r)
    return w_vec_neighbors(Int.(Meta.parse(x).args),r)
end

function compute_flip_graph(code_amalg,save_str;r=2)
    vector_to_idx = Dict(code_amalg[k][1] => k for k in 1:size(code_amalg,1))

    w_network = SimpleGraph(size(code_amalg,1))
    code_amalg = [code_amalg[i][1] for i = 1:size(code_amalg,1)]

    block_len = 1000 # to save memory
    nloop = Int(round(length(code_amalg)/block_len))


    for i = 1:nloop
        k_start = (i-1)*block_len + 1
        k_end = minimum([i*block_len ;length(code_amalg)])
        w_vec_nb = pmap(x->flip_core_parallel(x,r),code_amalg[k_start:k_end])
        for k = k_start:k_end
            for nb in w_vec_nb[k-k_start+1]
                if haskey(vector_to_idx,string(nb))
                    add_edge!(w_network,vector_to_idx[string(nb)],k)
                end
            end
        end
        println("Done ", i*block_len, "out of ", length(code_amalg))
    end


    savegraph( save_str*".lgz", w_network)
    df = DataFrame(codes = [code_amalg[k] for k in 1:size(code_amalg,1)],
        index = [k for k in 1:size(code_amalg,1)])
    CSV.write(save_str*".txt",  df)
end
