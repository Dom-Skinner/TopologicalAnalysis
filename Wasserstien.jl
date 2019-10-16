using GraphPlot
using CSV,DataFrames
using LightGraphsFlows
import LightGraphs
const lg = LightGraphs
using SparseArrays
using Base.Threads
using Clp: ClpSolver # use your favorite LP solver here


include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/VoronoiTools.jl")
using .VoronoiTools
include("/Users/Dominic/Documents/2d Cells/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure


function ret_weights(dict_w,N,W_code_to_idx,vmap)
    weight_arr = zeros(Float64,N)
    key_arr = collect(keys(dict_w))
    for k in key_arr
        weight_arr[W_code_to_idx[k]] = dict_w[k]
    end
    println(sum(weight_arr))
    weight_arr = weight_arr[vmap]
    println(sum(weight_arr))
    return weight_arr./sum(weight_arr)
end

function load_w_graph()
    Data_dir = "/Users/Dominic/Documents/2d Cells/Data/"
    g = lg.loadgraph(Data_dir*"w_network.lgz")
    #gplot(g)
    dat_in = CSV.read(Data_dir*"w_network.txt")
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
            capacity[Nv+1,i] = S[i]
        elseif S[i] < - 5e-6
            lg.add_edge!(g,i,Nv+2)
            demand[i,Nv+2] = -S[i]
            capacity[i,Nv+2] = 1.

        end
    end

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
    return d
end
g,vmap,N,W_code_to_idx,W_idx_to_code = load_w_graph()

Data_dir = "/Users/Dominic/Documents/2d Cells/Data/"
w_ells = readin(Data_dir*"Ells/Ells_",10)
w_spheres = readin(Data_dir*"Spheres/Spheres_",10)
w_PV = readin(Data_dir*"PV/PoissonVoronoi_",10)
w_exp = readin(Data_dir*"ExpData/Exp_",32)

weight = []
color = []
for i = 1:1
    push!(weight,ret_weights(w_ells[i],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_spheres[i],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_PV[1],N,W_code_to_idx,vmap))
    push!(weight,ret_weights(w_exp[1],N,W_code_to_idx,vmap))
    append!(color,[1;2;3;4])
end

d = distance_mat(g,weight)

using MultivariateStats
using Plots

Out = permutedims(classical_mds(d,2))
scatter(Out[:,1],Out[:,2],color=[1;1;2;2;3;3;4;4],leg=false)
