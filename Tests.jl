using Random

include("/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

Data_dir = "/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/"

Random.seed!(1234)
Positions = [randn(3) for k in 1:1000]
p, simplices, neighbours, edge_index = Delaunay_find(Positions)


function add_labeled!(Orig::Array{Int64,2},New::Array{Int64,2},
    dict_ll::Int64,ll::Int64,Label_arr::BitArray{1},
    Nzero::Array{Int64,1},s_code::Array{Int64,1})

    @inbounds for j = 1:size(Orig,2),i=1:size(Orig,1)
        if dict_ll == Orig[i,j]
            New[i,j] = ll
            Nzero[i] -= 1
            Label_arr[i] = (Nzero[i] == 1)
            s_code[i] += 1 << ll #2^ll
        end
    end
end

function fill_labels!(New::Array{Int64,2},idx_saved::Array{Int64,1},
            Orig::Array{Int64,2},Label_arr::BitArray{1},
            num_zero::Array{Int64,1},s_code::Array{Int64,1})

    @inbounds for n = 5:length(idx_saved)
        # Take the minimum val of sum_t for all potential choices
        min_val = typemax(Int64)
        idx = 0
        @inbounds for i in (1:size(Label_arr,1))
            if Label_arr[i] && (s_code[i] < min_val)
                min_val = s_code[i]
                idx = i
            end
        end

        idx_saved[n] = Orig[idx,findfirst(iszero,New[idx,:])]
        add_labeled!(Orig,New,idx_saved[n],n,Label_arr,num_zero,s_code)
    end

end

function compare_fun!(A,B)
    @inbounds for i = 1:length(A)
        if  A[i] < B[i]
            return
        elseif A[i] > B[i]
            A .= B
            return
        end
    end
    return
end

function find_to_perm!(to_perm,nbhd,M)
    @inbounds for i = 1:M
        kk = 1
        for j = 1:4
            if k != nbhd[i,j]
                to_perm[i,kk] = j
                kk += 1
            end
        end
    end
end

k = 6
k_nbhd = simplices[[ k ∈ simplices[i,:] for i in 1:size(simplices,1) ],:]
function topological_vec(nbhd,central_vertex;r=1)
    # This function takes a neighborhood of simplices around a central vertex
    # and gives it a cannonical labeling.

    N = length(unique(nbhd))
    if N > 63
        error("Too many to store in int64")
    end

    if r != 1
        error("Todo")
        # To fix this, need to modify the function find_to_perm! and M, so that
        # the for loop goes over all simplices neighboring k rather than all
        # possible simplices in the local network
    end

    M = size(nbhd,1)
    simplex_n = zeros(Int64,size(nbhd))
    to_label = falses(M)
    num_zero = zeros(Int64,M)
    sim_code = zeros(Int64,M)
    Num2idx = zeros(Int64,N)
    Topological_vec = typemax(Int64)*ones(Int64,M)
    perms = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2],[3;2;1]]
    to_perm = zeros(Int64,M,3)
    find_to_perm!(to_perm,nbhd,M)

    for i = 1:M, p in 1:6
        Num2idx .= 0
        to_label .= false
        simplex_n .= 0
        num_zero .= 4
        sim_code .= 0

        Num2idx[1] = central_vertex
        add_labeled!(nbhd,simplex_n,Num2idx[1],1,to_label,num_zero,sim_code)

        for j in 1:3
            Num2idx[j+1] = nbhd[i,to_perm[i,perms[p][j]]]
            add_labeled!(nbhd,simplex_n,Num2idx[j+1],j+1,to_label,num_zero,sim_code)
        end
        fill_labels!(simplex_n,Num2idx,nbhd,to_label,num_zero,sim_code)
        sort!(sim_code)
        compare_fun!(Topological_vec,sim_code)
    end
    #println(Topological_vec)

    return

end

using BenchmarkTools
using Profile
Profile.clear()

@btime topological_vec(k_nbhd,k)
Profile.print()
Juno.profiler()
