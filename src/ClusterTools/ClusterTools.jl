module ClusterTools
using LightGraphs
using LinearAlgebra
using Arpack
using Clustering: kmeans
using StatsBase: var,mean,countmap

include("../PointCloudTools/PointCloudTools.jl")
using .PointCloudTools

function topological_cluster(Positions::Array;k=2,r=3,α=0)
    p, simplices, neighbours, edge_index, α_val,α = Delaunay_find(Positions, α)
    g_full = graph_construct(simplices,size(p,1))

    if is_connected(g_full)
        return topological_cluster(g_full,k=k,r=r)
    end

    # If we reach here, the graph is not connected.
    g_con_idx = connected_components(g_full)
    g_con_idx = g_con_idx[argmax(length.(g_con_idx))] # work with largest connected component

    g_part,vmap = induced_subgraph(g_full,g_con_idx)
    labels_part = topological_cluster(g_part,k=k,r=r)
    labels_full = zeros(Int64,nv(g_full))
    labels_full[vmap] .= labels_part


    # Now take a graph that is connected, and assign the remaining points to their
    # most common neighboring labels
    g_full = minimally_connected_graph(Positions,α)

    for j = 1:length(labels_full)
        if labels_full[j] > 0
            continue
        end

        for r_val = 2:10
            N = unique(neighborhood(g_full,j,r_val))
            local_labs = filter(x->x>0,labels_full[N])
            if length(local_labs) > 0
                N = countmap(local_labs)
                v = collect(values(N))
                idx = collect(keys(N))
                labels_full[j] = idx[argmax(v)]
                break
            end
        end


    end

    return labels_full
end


function minimally_connected_graph(Positions,α;iter_max_val=10)
    # If the alpha complex not connected, try again with a larger alpha
    # Ideally don't have to use α=Inf since this may contain many spurious connections

    # try up to (2^iter_max_val)*α before going to Inf
    for iter = 1:iter_max_val

        α = 2*α
        p, simplices, neighbours, edge_index, α_val,α = Delaunay_find(Positions, α)
        g_full = graph_construct(simplices,size(p,1))

        if is_connected(g_full)
            return g_full
        end
    end

    p, simplices, neighbours, edge_index, α_val,α = Delaunay_find(Positions, Inf)
    return graph_construct(simplices,size(p,1))

end

function topological_cluster(g_full::Graph;k=2,r=3)

    renorm = x-> (x .- mean(x)) ./ norm(x .- mean(x))

    # compute topological statistics
    nbh = [neighborhood(g_full,i,r) for i = 1:nv(g_full)]
    n = length.(nbh)
    r1 = renorm([var(n[nb]) for nb in nbh])
    r2 = renorm([mean(n[nb]) for nb in nbh])

    # compute diffusion statistics
    λ,evecs = eigs(laplacian_matrix(g_full),nev=minimum([6;k]), which=:SR,maxiter=750)

    # combine and take kmeans
    vec_tot = hcat([r1,r2,evecs[:,2:end]]...)
    labels =  kmeans(collect(transpose(vec_tot)), k).assignments
    make_connected!(labels,g_full)

    return labels
end
#=
function edge_reassign!(labels,g_full;iter=10)
    r = 2
    code_tot,_ = weinberg2D_core(g_full,nv(g_full),r)

    for count = 1:iter
        cluster_top = [countmap(code_tot[labels .== i]) for i = 1:k]

        boundary_assign = [labels[neighborhood(g_full,i,1)] for i = 1:nv(g_full)]
        nbh_counts =  countmap.(boundary_assign)
        boundary_assign = [nbh_counts[i][labels[i]]/ sum(collect(values(nbh_counts[i]))) for i = 1:length(nbh_counts)]


        for i = 1:length(labels)
            if boundary_assign[i]  < 0.5
                n = [haskey(clust,code_tot[i]) ? clust[code_tot[i]] : 0 for clust in cluster_top]
                nmax = maximum(n)
                nmax = findall(n .== nmax)
                labels[i] = nmax[argmax([nbh_counts[i][j] for j in nmax])]
            end
        end
    end
end
=#

function make_connected!(labels,g_full)
    # takes the clustering done by topological_cluster, and reassigns points that
    # are isolated, i.e. not in a major connected component.
    for l in unique(labels)
        # Find the connected components of the subgraph for a particular label
        g_sub,vmap = induced_subgraph(g_full,findall(l .== labels))
        comp = connected_components(g_sub)
        comp_len = length.(comp)

        # If not in a major connected component reassign
        for i = 1: length(comp)
            if comp_len[i] < 0.5*maximum(comp_len)
                # Find the labels of the neighbors of this component in the full
                # graph, and assign the component to be equal to the most common one.
                orig_index = vmap[comp[i]]
                N = unique(vcat([neighborhood(g_full,o,1) for o in orig_index]...))
                N = countmap(labels[setdiff(N,orig_index)])
                v = collect(values(N))
                idx = collect(keys(N))
                labels[orig_index] .= idx[argmax(v)]
            end
        end
    end
end

export  topological_cluster
end
