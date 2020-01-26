
function add_labeled!(Orig::Array{Int64,2},to_label::Array{Int64,1},
    dict_ll::Int64,ll::Int64, Nzero::Array{Int64,1},s_code::Array{Int64,1})
    @inbounds for j = 1:size(Orig,2),i=1:size(Orig,1)
        if dict_ll == Orig[i,j]
            to_label[i] -= j
            Nzero[i] -= 1
            s_code[i] += 1 << ll #2^ll
        end
    end
end

function fill_labels!(to_label::Array{Int64,1}, Orig::Array{Int64,2},
            num_zero::Array{Int64,1},s_code::Array{Int64,1},N::Int64)

    @inbounds for n = 5:N
        # Take the minimum val of s_code for all potential choices
        min_val = typemax(Int64)
        idx = 0
        @inbounds for i in 1:length(num_zero)
            if (num_zero[i] == 1) && (s_code[i] < min_val)
                min_val = s_code[i]
                idx = i
            end
        end

        add_labeled!(Orig,to_label, Orig[idx,to_label[idx]],n,num_zero,s_code)
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

function find_to_perm!(to_perm,nbhd,central_vertex,M)
    @inbounds for i = 1:M
        kk = 1
        for j = 1:4
            if central_vertex != nbhd[i,j]
                to_perm[i,kk] = j
                kk += 1
            end
        end
    end
end

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

    # num_zero tracks how many of the vertices of a given simplex are zero
    # (meaning unlabeled). If num_zero = 1, then which_to_label will be the
    # index in 1:4 that needs labelling. Sim code contains the
    # labelings in the form sim_code = 2^n1 + 2^n2 + ... Topological_vec is the
    # smallest out of all ordered sim codes.
    M = size(nbhd,1)
    num_zero = zeros(Int64,M)
    which_to_label = zeros(Int64,M)
    sim_code = zeros(Int64,M)
    Topological_vec = typemax(Int64)*ones(Int64,M)
    perms = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2],[3;2;1]]
    to_perm = zeros(Int64,M,3)
    find_to_perm!(to_perm,nbhd,central_vertex,M)

    for i = 1:M, p in 1:6
        which_to_label .= 10
        num_zero .= 4
        sim_code .= 0

        add_labeled!(nbhd,which_to_label,central_vertex,1,num_zero,sim_code)

        for j in 1:3
            add_labeled!(nbhd,which_to_label, nbhd[i,to_perm[i,perms[p][j]]],
                            j+1,num_zero,sim_code)
        end
        fill_labels!(which_to_label,nbhd,num_zero,sim_code,N)
        sort!(sim_code)
        compare_fun!(Topological_vec,sim_code)
    end

    return

end
