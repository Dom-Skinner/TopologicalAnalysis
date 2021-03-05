function find_nbhd(simplices,k)
    N = size(simplices,1)
    contains_k = falses(N)
    #@inbounds for j = 1:4, i=1:N
    for j = 1:4, i=1:N
        if k == simplices[i,j]
            contains_k[i] = true
        end
    end
    simplices[contains_k,:]
end

#=
function add_labeled!(Orig::Array{Int64,2},to_label::Array{Int64,1},
    dict_ll::Int64,ll::Int64, Nzero::Array{Int64,1},s_code::Array{Int64,1})
    #@inbounds for j = 1:size(Orig,2),i=1:size(Orig,1)
    for j = 1:size(Orig,2),i=1:size(Orig,1)
        if dict_ll == Orig[i,j]
            to_label[i] -= j
            Nzero[i] -= 1
            s_code[i] += 1 << ll #2^ll
        end
    end
end
=#
function add_labeled!(Orig::Array{Int64,2},to_label::Array{Int64,1},
    coord_arr,ll::Int64, Nzero::Array{Int64,1},s_code::Array{Int64,1})
    for c in coord_arr
            to_label[c[1]] -= c[2]
            Nzero[c[1]] -= 1
            s_code[c[1]] += 1 << ll #2^ll
    end
end


function deg_case(nbhd,idx_loc,i_save,p_save,to_perm,central_vertex)

    # If the standard algorithm fails to label the simplicial complex, we need
    # to make another guess. We have to try all possible guesses from all possible
    # initial starting points in i_save, p_save.

    tvec_tot = []
    N = length(unique(nbhd))
    M = size(nbhd,1)
    num_zero = zeros(Int64,M)
    which_to_label = zeros(Int64,M)
    sim_code = zeros(Int64,M)
    Topological_vec = typemax(Int64)*ones(Int64,M)

    perms_ = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2],[3;2;1]]
    for kk = 1:length(i_save)
        i = i_save[kk]
        p = p_save[kk]
        which_to_label .= 10
        num_zero .= 4
        sim_code .= 0
        Topological_vec .= typemax(Int64)
        Topological_vec[1] = 30

        add_labeled!(nbhd,which_to_label,idx_loc[central_vertex],1,num_zero,sim_code)

        for j in 1:3
            add_labeled!(nbhd,which_to_label, idx_loc[nbhd[i,to_perm[i,perms_[p][j]]]],
                            j+1,num_zero,sim_code)
        end
        num_zero[i] = -1
        fill_labels(which_to_label,nbhd,num_zero,sim_code,N,Topological_vec,idx_loc)

        # The above gets us to the starting point of i_save,p_save, where no more
        # simplices can be labeled. We now try all possible guesses



        N = length(unique(nbhd))
        already_labeled = unique(nbhd[num_zero .== -1,:])
        to_label = setdiff(unique(nbhd),already_labeled)

        n_choice = minimum(num_zero[num_zero .> 0])
        simp_choice = findall(x-> x==n_choice, num_zero)
        if n_choice == 2
            perms = [[1;2],[2;1]]
        elseif n_choice == 3
            perms = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2],[3;2;1]]
        else
            error()
        end

        which_to_label_cp = copy(which_to_label)
        num_zero_cp = copy(num_zero)
        sim_code_cp = copy(sim_code)
        Topological_vec_cp = copy(Topological_vec)


        for s in simp_choice
            for p in perms
                which_to_label_cp .= which_to_label
                num_zero_cp .= num_zero
                sim_code_cp .= sim_code
                Topological_vec_cp .= Topological_vec

                lab_next = nbhd[s,findall(x -> x∈ to_label, nbhd[s,:])]
                add_labeled!(nbhd,which_to_label_cp,idx_loc[lab_next[p[1]]],
                            length(already_labeled)+1,num_zero_cp,sim_code_cp)
                add_labeled!(nbhd,which_to_label_cp,idx_loc[lab_next[p[2]]],
                            length(already_labeled)+2,num_zero_cp,sim_code_cp)
                            num_zero_cp[s] = -1

                Topological_vec_cp[findlast(x-> x< typemax(Int64),Topological_vec_cp) .+ 1 ] = sim_code_cp[s]

                fill_labels(which_to_label_cp,nbhd,num_zero_cp,sim_code_cp,N,Topological_vec_cp,idx_loc)

                push!(tvec_tot,copy(Topological_vec_cp))
            end
        end
    end
    Topological_vec .= minimum(tvec_tot)
end


function fill_labels(to_label::Array{Int64,1}, Orig::Array{Int64,2},
            num_zero::Array{Int64,1},s_code::Array{Int64,1},N::Int64,Topological_vec,idx_loc)
    # This is the core function of the labeling algorithm. It assumes that the
    # first simplex has been chosen and labeled. The rest of the labeling now
    # procedes deterministically, but will terminate if the labeling is not less
    # than the current value of Topological_vec
    N_sim_labeled = sum(num_zero .== -1)

    flag = 0
    n_l  = length(unique(Orig[num_zero .== -1,:])) + 1
    for n = n_l:N
    #@inbounds for n = 5:N
        # Take the minimum val of s_code for all potential choices
        min_val = typemax(Int64)
        idx = 0
        #@inbounds for i in 1:length(num_zero)
        for i in 1:length(num_zero)
            if (num_zero[i] == 1) && (s_code[i] < min_val)
                min_val = s_code[i]
                idx = i
            end
        end
        if idx == 0
            if Topological_vec[N_sim_labeled+1] < typemax(Int64)
                if flag == 1
                    Topological_vec[N_sim_labeled+1:end] .= typemax(Int64)
                    return true
                else
                    return false
                end
            else
                return true
            end
        end

        # now the next labeled point has been chosen, update the arrays
        add_labeled!(Orig,to_label, idx_loc[Orig[idx,to_label[idx]]],n,num_zero,s_code)

        # Find the simplices that have just been fully labeled
        just_labeled = Int64[]
        for i = 1:length(num_zero)
            if (num_zero[i] == 0)
                num_zero[i] = - 1
                push!(just_labeled,s_code[i])
            end
        end
        sort!(just_labeled)

        # If the new topological vector is less than the original, update so the
        # new topological vector replaces the old. If it is greater than the
        # original terminate this function early
        if flag != 1
            for i = 1:length(just_labeled)
                if just_labeled[i] < Topological_vec[N_sim_labeled+i]
                    flag = 1
                    break
                elseif just_labeled[i] > Topological_vec[N_sim_labeled+i]
                    return false
                end
            end
        end
        if flag == 1
            for i = 1:length(just_labeled)
                Topological_vec[N_sim_labeled+i] = just_labeled[i]
            end
        end

        N_sim_labeled += length(just_labeled)
    end
    return true
end

function compare_fun!(A,B)
    #@inbounds for i = 1:length(A)
    for i = 1:length(A)
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
    # This function finds which vertices are to be permuted
    @inbounds for i = 1:M
        kk = 1
        for j = 1:4
            if central_vertex != nbhd[i,j] && kk < 4
                to_perm[i,kk] = j
                kk += 1
            elseif (kk == 4) && central_vertex != nbhd[i,j]
                to_perm[i,1] = 0
            end
        end
    end
end

function topological_vec(k_nbhd,central_vertex;r=1)
    # This function takes a neighborhood of simplices around a central vertex
    # and gives it a cannonical labeling.

    # First make a copy of k_nbhd and relabel the vertices 1:N
    nbhd = copy(k_nbhd)
    nbhd_unique = unique(nbhd)
    N = length(nbhd_unique)
    if N > 62
        error("Too many to store in int64")
    end

    for i in 1:length(nbhd_unique)
        for k = 1:size(k_nbhd,2), j = 1:size(k_nbhd,1)
            if k_nbhd[j,k] == nbhd_unique[i]
                nbhd[j,k] = i
            end
        end
    end
    central_vertex = findfirst(isequal(central_vertex),nbhd_unique)


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
    Topological_vec[1] = 30
    perms = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2],[3;2;1]]
    to_perm = zeros(Int64,M,3)
    find_to_perm!(to_perm,nbhd,central_vertex,M)
    # idx_loc[i] contains an array of indices [x,y] s.t. nbhd[x,y] = i
    idx_loc = [findall(isequal(i),nbhd) for i in 1:N]
    i_save = []
    p_save = []
    for i = 1:M
        if to_perm[i,1] != 0
            for p in 1:6
                which_to_label .= 10
                num_zero .= 4
                sim_code .= 0

                add_labeled!(nbhd,which_to_label,idx_loc[central_vertex],1,num_zero,sim_code)

                for j in 1:3
                    add_labeled!(nbhd,which_to_label, idx_loc[nbhd[i,to_perm[i,perms[p][j]]]],
                                    j+1,num_zero,sim_code)
                end
                num_zero[i] = -1
                tvec_cp =  copy(Topological_vec)
                optimal = fill_labels(which_to_label,nbhd,num_zero,sim_code,N,Topological_vec,idx_loc)

                # For the degenerate case, we save all the initial configurations
                # that give rise to the best possible start
                if tvec_cp != Topological_vec
                    i_save = []
                    p_save = []
                end

                if optimal
                    push!(i_save,i)
                    push!(p_save,p)
                end


                #compare_fun!(Topological_vec,sim_code)

            end
        end
    end
    # if there are still unlabeled vertices, then we have to make another guess
    # to finish labeling. Do this with function deg_case
    if any(x-> x == typemax(Int64),Topological_vec)
        Topological_vec = deg_case(nbhd,idx_loc,i_save,p_save,to_perm,central_vertex)
    end
    return Topological_vec

end

function t_vec_to_simplex!(tvec,simplex)
    # Converts a topological vector back into the simplex form
    for i = 1:length(tvec)
        j = 4
        ttemp = tvec[i]
        for p2 = 62:-1:1
            if ttemp >= (1 << p2)
                simplex[i,j] = p2
                ttemp -= (1 << p2)
                j -= 1
            end
        end
    end
end
