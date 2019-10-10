using LightGraphs
using GraphPlot

include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/VoronoiTools.jl")
using .VoronoiTools

function central_node(g)
    N = nv(g)
    ret = []
    for i = 1:N
        if length(neighborhood(g,i,2)) == N
            push!(ret,i)
        end
    end
    if length(ret) == 1
        return ret[1]
    else
        println("Winging it on the central node")
        # TODO fix this...
        x = eigenvector_centrality(g)
        v,idx = findmax(x)
        return idx
    end
end

function tutte_embedding(g)
    # This function creates a Tutte embedding of the graph.
    fixed_vecs = []
    for e in edges(g)
        if length(intersect(neighbors(g,src(e)),neighbors(g,dst(e)))) == 1
            push!(fixed_vecs,src(e))
            push!(fixed_vecs,dst(e))
            push!(fixed_vecs,intersect(neighbors(g,src(e)),neighbors(g,dst(e)))[1])
            break
        end
    end
    if length(fixed_vecs) != 3
        error("TODO: code the degenerate case...")
        # degenerate case being when there are no edges attached to only a single triangle
    end
    A = float(adjacency_matrix(g))
    D = [sum(A[i,:]) for i in 1:size(A,1)]
    for i = 1:size(A,1)
       for j in 1:size(A,2)
           if A[i,j] != 0
               A[i,j] /= D[i]
           end
       end
    end
    for i = 1:size(A,1)
        A[i,i] -= 1
    end

    x = zeros(nv(g)) ;      y = zeros(nv(g))
    x[fixed_vecs[1]] = 0.;  y[fixed_vecs[1]] = 0.;
    x[fixed_vecs[2]] = 0.;  y[fixed_vecs[2]] = 1.;
    x[fixed_vecs[3]] = 1.;  y[fixed_vecs[3]] = 0.;
    bx = A*x ; by = A*y
    keep_idx = [k ∉ fixed_vecs for k in 1:nv(g)]
    bx = bx[keep_idx] ;  by = by[keep_idx]
    A = A[keep_idx,keep_idx]

    x[keep_idx] = - A \ bx
    y[keep_idx] = - A \ by
    return x,y,fixed_vecs
end

function circ_insert!(M,pair_set,val)
    idx1 = findfirst(x-> x ∈ pair_set , M)
    idx2 = findlast(x-> x ∈ pair_set , M)

    if idx2 - idx1 == 1
        insert!(M,idx2,val)
    else
        push!(M,val)
    end
end

function weinberg_flip(g,cent_node,order_mat,d,s,v1,v2)

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
    nbh  = neighborhood(g2,cent_node,2)
    g_ego, vmap = induced_subgraph(g2,nbh)
    vmap_inv = Dict(vmap[k] => k for k in 1:length(vmap))
    order_local = order_copy[vmap]
    for i = 1:length(vmap)
        total_order = map.(x -> get(vmap_inv, x, -1), order_copy[vmap[i]])
        order_local[i] = total_order[total_order .> 0]
    end
    weinberg_find!(code_tot,S_tot,1,g_ego,order_local)
    return code_tot[1]
end

function order_mat_find(g,x,y)
    order_mat = Vector{Array{Int64}}(undef,nv(g))
    for kk = 1:nv(g)
        N_list = neighbors(g,kk)
        theta = [atan(y[s]-y[kk], x[s]-x[kk]) for s in N_list]
        order_mat[kk] = N_list[sortperm(theta)]
    end
    return order_mat
end

function return_nbh_vertex(order_arr,v)
    next_ind = (i,m) -> (i+1)*(i<m) +(i == m)
    prev_ind = (i,m) -> (i-1)*(i>1) +m*(i == 1)

    d1 = findfirst(x->x==v, order_arr)
    v1 = order_arr[ next_ind(d1,length(order_arr))]
    v2 = order_arr[ prev_ind(d1,length(order_arr))]
    return v1,v2
end


w1 = [1, 2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 6, 1, 6, 2, 6, 7, 2, 7, 8, 2, 8, 3, 8, 9, 3, 9, 10, 3, 10, 4, 10, 11, 4, 11, 12, 4, 12, 5, 12, 11, 13, 14, 15, 16, 17, 18, 7, 18, 8, 18, 17, 8, 17, 9, 17, 16, 9, 16, 15, 9, 15, 10, 15, 14, 10, 14, 13, 10, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
g = SimpleGraph(maximum(w1))
for i = 1:(length(w1)-1)
    add_edge!(g,w1[i],w1[i+1])
end

#=
function plot_a(g,x,y)
    p = scatter(x,y)

    for e in edges(g)
        plot!(p,[x[src(e)];x[dst(e)]], [y[src(e)];y[dst(e)]],leg=false)
    end
    return p
end
plotly()
p = plot_a(g,x,y)
=#

x,y,fixed_vecs = tutte_embedding(g)
order_mat = order_mat_find(g,x,y)



cent_node = central_node(g)
ds  = dijkstra_shortest_paths(g,cent_node)
dist_to_cent = ds.dists
for e in edges(g)
    d = dst(e)
    s = src(e)

    v11,v21 = return_nbh_vertex(order_mat[s],d)
    v22,v12 = return_nbh_vertex(order_mat[d],s)

    if (v11 == v12) && (v21 == v22)
        #println("A Flippable triangle")
        if (dist_to_cent[d] == dist_to_cent[s]) &&
            (abs(dist_to_cent[v11]-dist_to_cent[v21]) == 2)
            #println("Don't flip!")
        else
            w_vec = weinberg_flip(g,cent_node,order_mat,d,s,v11,v21)
            println(w_vec)
        end
    end
end
