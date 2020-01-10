using CSV,DataFrames
using LightGraphsFlows
import LightGraphs
const lg = LightGraphs
using SparseArrays
using Base.Threads
using Clp: ClpSolver # use your favorite LP solver here



function ret_weights(dict_w,N,W_code_to_idx,vmap)
    # Find the distribution of networks in the connected comp of the flip graph
    weight_arr = zeros(Float64,N)
    key_arr = collect(keys(dict_w))
    for k in key_arr
        if haskey(W_code_to_idx,k)
            weight_arr[W_code_to_idx[k]] = dict_w[k]
        end
    end
    weight_arr = weight_arr[vmap]
    return weight_arr./sum(weight_arr)
end

function load_w_graph(network_save_file)
    # Wrapper to load the flip graph, and to reduce to the largest connected comp
    g = lg.loadgraph(network_save_file*".lgz")
    #gplot(g)
    dat_in = CSV.read(network_save_file*".txt")
    W_code_to_idx = Dict(dat_in.codes .=> dat_in.index)
    W_idx_to_code = Dict(dat_in.index .=> dat_in.codes)

    ccomp = lg.connected_components(g)
    idx = argmax([length(c) for c in ccomp])
    g_con, vmap = lg.induced_subgraph(g, ccomp[idx])

    return g_con,vmap,lg.nv(g),W_code_to_idx,W_idx_to_code
end

function W_dist(g_undirected,S)

    # First create a directed version of g_undirected with edges in both
    # directions and capacities/costs equal to 1
    Nv = lg.nv(g_undirected)
    g = lg.DiGraph(Nv+2) # Create a flow-graph
    w = spzeros(Nv+2,Nv+2)
    capacity = spzeros(Nv+2,Nv+2)

    for e in lg.edges(g_undirected)
        lg.add_edge!(g,lg.src(e),lg.dst(e))
        lg.add_edge!(g,lg.dst(e),lg.src(e))
        capacity[lg.src(e),lg.dst(e)] = 1.
        capacity[lg.dst(e),lg.src(e)] = 1.
        w[lg.src(e),lg.dst(e)] = 1.
        w[lg.dst(e),lg.src(e)] = 1.
    end

    # Take the Nv+1th vertex as a source and Nv+2th as a sink.
    demand = spzeros(Nv+2,Nv+2)
    for i in 1:Nv
        if S[i] > 0
            lg.add_edge!(g,Nv+1,i)
            capacity[Nv+1,i] = S[i]*(1+5e-5)
        elseif S[i] < - 5e-6
            lg.add_edge!(g,i,Nv+2)
            demand[i,Nv+2] = -S[i]*(1-5e-5)
            capacity[i,Nv+2] = 1.

        end
    end
    # call min cost flow
    flow = mincost_flow(g,spzeros(lg.nv(g)), capacity , w, ClpSolver(),
                edge_demand=demand, source_nodes=[Nv+1], sink_nodes=[Nv+2])

    return sum(abs.(flow[1:Nv,1:Nv]))
end


function distance_mat(g,weight)
    d = zeros(length(weight),length(weight))
    for i = 1:size(d,1)
    Threads.@threads for j=(i+1):size(d,2)
                    W = weight[i] .- weight[j]
                    d[i,j] = W_dist(g,sign(sum(W))*W + W*(sign(sum(W)) == 0))
                    d[j,i] = d[i,j]
                    #println("i= ",i, " j = ", j)
                    #println(d[i,j])
                end
        println(100*i/size(d,1))
    end
    println("Catching exceptions")
    #This is a HACK
    for i = 1:size(d,1)
    Threads.@threads for j=(i+1):size(d,2)
                    if d[i,j] == 0
                        W = weight[i] .- weight[j]
                        d[i,j] = W_dist(g,-sign(sum(W))*W + W*(sign(sum(W)) == 0))
                        d[j,i] = d[i,j]
                    end
                end
        println(100*i/size(d,1))
    end

    return d
end

function calculate_distance_matrix(network_save_file,w_vec_in)
    # This function is a wrapper for all other functions in this file
    # from n dictionaries in and the path to the load_graph file it returns the
    # n by n distance matrix
    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]

    return distance_mat(g,weight)
end

function calculate_distance_matrix_parallel(network_save_file,w_vec_in)
    # If we want to write our own parallel distance matrix function
    g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph(network_save_file)
    weight = [ret_weights(w_vec_in[i],N,W_code_to_idx,vmap) for i in 1:length(w_vec_in)]

    return g,weight
end

#=
weight = []
for i = 1:7
    push!(weight,ret_weights(w_ells[i],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_spheres[i],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_PV[i],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_exp[i],N,W_code_to_idx,vmap))
end
=#

#=
using MultivariateStats
using Plots

Out = permutedims(transform(fit(MDS, d, maxoutdim=2, distances=true)))
savefig(plot_res(Out),"MDS.pdf")

function make_heatmap(d)
    idx = [collect(1:4:28);collect(2:4:28);collect(3:4:28);collect(4:4:28)]
    d2 = similar(d)
    for i=1:28
        for j = 1:28
             d2[i,j] = d[idx[i],idx[j]]
         end
     end
    p = heatmap(d2,colorbar_title="W distance",xlabel="Sample number",
            ylabel="Sample number", dpi = 350,c=:BuGn)
    savefig(p, "Corellation.png")
end
function plot_res(Out)
    names = ["Ellipsoids","Spheres","Poisson-Voronoi","Experiments"]
    p = scatter()
    for i = 1:4
        idx = 0:4:(size(Out,1)-4)
        scatter!(p,Out[i.+idx,1],Out[i.+idx,2],color=i,label=names[i])
    end
    return p
end
=#
