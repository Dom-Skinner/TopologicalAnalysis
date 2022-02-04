using LightGraphs
using LinearAlgebra
using Arpack
using Clustering: kmeans
using StatsBase: var,mean,countmap,median

function custom_NN(vec,point)
    r2 = sum( (vec .- permutedims(point)).^2,dims=2)
    return argmin(r2)[1]
end

function topological_cluster(mot,fg,pos;k=2,r=3,alpha=0,idx_in=[])


    L = laplacian_matrix(fg.g)
    D = 0.04*ones(nv(fg.g))
    for m in mot.tvec
        if haskey(fg.motif_code,m)
            D[fg.motif_code[m]] = D[fg.motif_code[m]] + 1
        end
    end
    DI = spdiagm(0=>1 ./D)

    λ,evecs = eigs(L*DI, which=:SR,maxiter=5000)
    t1 = [haskey(fg.motif_code,m) ? real(evecs[fg.motif_code[m],1]) : NaN for m in mot.tvec]
    t2 = [haskey(fg.motif_code,m) ? real(evecs[fg.motif_code[m],2]) : NaN for m in mot.tvec]
    t3 = [haskey(fg.motif_code,m) ? real(evecs[fg.motif_code[m],3]) : NaN for m in mot.tvec]
    tvec = hcat([t1,t2,t3]...)
    println(size(tvec))
    return remove_nans(tvec,pos,alpha)

    g = minimally_connected_graph(pos,alpha/2)
    nev=minimum([6,k])
    return topological_evec(g,r=r,nev=nev)


    # First take connected component of alpha complex
    g_part,vmap = connected_subgraph(pos,alpha)
    # Then find the topological + diffusion coords
    nev=minimum([6,k])
    vec_tot = topological_evec(g_part,r=r,nev=nev)
    # Apply k-means
    if length(idx_in)>0
        kmeans_res = kmeans(collect(transpose(vec_tot)), k,init=idx_in,maxiter=1)
    else
        kmeans_res = kmeans(collect(transpose(vec_tot)), k)
    end


    idx = [custom_NN(vec_tot,kmeans_res.centers[:,i]) for i = 1:k]
    idx = vmap[idx]

    labels =  kmeans_res.assignments
    # Make the k-means clustering connnected in physical space
    labels = make_connected(labels,g_part)
    # Finally, assign the remaining points not connected to the main component
    labels = assign_remainders(vmap,labels,pos,alpha)
    return labels,idx

end

function remove_nans(tvec::Array,pos::Array,α)

    g_full = minimally_connected_graph(pos,α)

    for j = 1:size(tvec,1)
        if !isnan(tvec[j,1])
            continue
        end

        N = unique(neighborhood(g_full,j,3))
        N = filter(x->!isnan(tvec[x,1]),N)
        for k = 1:size(tvec,2)
            tvec[j,k] = median(tvec[N,k])
        end
    end

    return tvec
end


function connected_subgraph(pos,alpha)
    # return the larges connected component of the alpha complex
    p, simplices, neighbours, edge_index, α_val,α_used = Delaunay_find(pos, alpha)
    g_full = graph_construct(simplices,size(p,1))

    g_con_idx = connected_components(g_full)
    g_con_idx = g_con_idx[argmax(length.(g_con_idx))] # work with largest connected component

    return induced_subgraph(g_full,g_con_idx)
end


function assign_remainders(vmap::Array,labels_part::Array,pos::Array,α)
    # Now take a graph that is connected, and assign the remaining points to their
    # most common neighboring labels

    g_full = minimally_connected_graph(pos,α)

    if nv(g_full) == length(labels_part)
        return labels_part
    end
    labels_full = zeros(Int64,nv(g_full))
    labels_full[vmap] .= labels_part

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

function topological_evec(g_full::Graph;r=3,nev=6)
    # For a given graph, compute the topological + spectral vector which we will
    # apply k-means to.
    @assert is_connected(g_full) # If g_full not connected, e-vecs will just be
    # connected components
    renorm = x-> (x .- mean(x)) ./ norm(x .- mean(x))

    # compute topological statistics
    nbh = [neighborhood(g_full,i,r) for i = 1:nv(g_full)]
    n = length.(nbh)
    r1 = renorm([var(n[nb]) for nb in nbh])
    r2 = renorm([mean(n[nb]) for nb in nbh])

    #scale = norm(pos[:,1])
    # compute diffusion statistics
    λ,evecs = eigs(laplacian_matrix(g_full),nev=nev, which=:SR,maxiter=750)

    # combine and take kmeans
    return hcat([r1,r2,evecs[:,2:end]]...)#hcat([r1,r2,pos./scale]...)#
end


function make_connected(l_orig,g_full)
    # takes the clustering done by topological_cluster, and reassigns points that
    # are isolated, i.e. not in a major connected component.
    labels = copy(l_orig)
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
    return labels
end
