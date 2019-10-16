using LightGraphs
using GraphPlot
using DataFrames,CSV
include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/VoronoiTools.jl")
using .VoronoiTools
include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure


function center_test(g,v)
    g_ego1,vmap = induced_subgraph(g,neighborhood(g,v,1))

    for e in edges(g_ego1)
        if length(intersect(neighbors(g,vmap[dst(e)]),neighbors(g,vmap[src(e)]))) == 1
            return false
        end
    end
    return true
end




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
        ret_idx = findall([center_test(g,r) for r in ret])
        if length(ret_idx) == 1
            return ret[ret_idx[1]]
        end

        println("Trying new method")
        edge_pts = []
        for e in edges(g)
            if length(intersect(neighbors(g,dst(e)),neighbors(g,src(e)))) == 1
                push!(edge_pts,src(e))
                push!(edge_pts,dst(e))
            end
        end
        unique!(edge_pts)
        rret = [];
        for r in ret
            if length(intersect(neighbors(g,r),edge_pts)) == 0
                 push!(rret,r)
            end
        end
        if length(rret) != 1
            savegraph("debug.lgz", g)
            println("Couldn't find central node (so guessed)")
        end
        return rret[1]

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
    if (minimum([length(neighbors(g_ego,k)) for k = 1:nv(g_ego)]) < 3) || (has_edge(g,v1,v2))
        # don't flip!!
        return [-1]
    end
    weinberg_find!(code_tot,S_tot,1,g_ego,order_local)
    if S_tot[1] == -1
        println(order_mat)
        println(order_copy)
        println(order_local)
        savegraph( "debug2.lgz", g)
        error("..")
    end
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


function w_vec_neighbors(w_in)
    w_neighbors = []
    g = SimpleGraph(maximum(w_in))
    for i = 1:(length(w_in)-1)
        add_edge!(g,w_in[i],w_in[i+1])
    end
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
                push!(w_neighbors, weinberg_flip(g,cent_node,order_mat,d,s,v11,v21))
            end
        end
    end
    return w_neighbors
end


function link_w!(g_weinberg,d_weinberg,w1,w2)
    if !haskey(d_weinberg,w1)
          d_weinberg[w1] = nv(g_weinberg)+1
          add_vertex!(g_weinberg)
    end
    if !haskey(d_weinberg,w2)
          d_weinberg[w2] = nv(g_weinberg)+1
          add_vertex!(g_weinberg)
    end
    add_edge!(g_weinberg,d_weinberg[w1],d_weinberg[w2])
end


Data_dir = "/Users/Dominic/Documents/2d Cells/Data/"
w_tot = readin(Data_dir*"Ells/Ells_",10)
append!(w_tot,readin(Data_dir*"PV/PoissonVoronoi_",10))
append!(w_tot,readin(Data_dir*"Spheres/Spheres_",10))
append!(w_tot,readin(Data_dir*"ExpData/Exp_",32))

code_amalg = amalg2(w_tot)
vector_to_idx = Dict(code_amalg[k][1] => k for k in 1:size(code_amalg,1))

w_network = SimpleGraph(size(code_amalg,1))
for i = 1:size(code_amalg,1)
    w = code_amalg[i][1]
    w_num = Meta.parse(w)
    w_nb = w_vec_neighbors(Int.(w_num.args))
    for nb in w_nb
        if haskey(vector_to_idx,string(nb))
            add_edge!(w_network,vector_to_idx[string(nb)],i)
        end
    end
    if mod(i,100) == 0
        println(i/size(code_amalg,1) * 100)
    end
end

savegraph( Data_dir*"w_network.lgz", w_network)
df = DataFrame(codes = [code_amalg[k][1] for k in 1:size(code_amalg,1)],
    index = [k for k in 1:size(code_amalg,1)])
CSV.write(Data_dir*"w_network.txt",  df)

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


#=
g = loadgraph(Data_dir*"w_network.lgz")
gplot(g)
=#
