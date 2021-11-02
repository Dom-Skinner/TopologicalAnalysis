using LightGraphs
using SparseArrays
using LinearAlgebra
using IterativeSolvers

function rem_self_edges!(g)
    # Removes self loops
    for e in collect(edges(g))
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
                    #J = D' * (L \ W) # Gaussian elimination too memory intensive
                    J = D' * minres(L , W)
                    d[i,j] = sum(abs.(J))
                    d[j,i] = d[i,j]
            end
            println(100*i/size(d,1))
    end
    return d # the distnace matrix
end
