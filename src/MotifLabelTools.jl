using LightGraphs
using Statistics: mean
using StatsBase: countmap
using LinearAlgebra: norm

struct MotifArray
           idx::Array
           tvec::Array
           dim::Int
           r::Int
           regions::Union{Array, Missing}
end

struct MotifDist
           map::Dict
           dim::Int
           r::Int
end

function MotifArray(idx::Array, tvec::Array, dim::Int, r::Int)
    return MotifArray(idx, tvec, dim, r, missing)
end

function MotifArray(m::MotifArray, regions::Array)
    return MotifArray(m.idx, m.tvec, m.dim, m.r, regions)
end

function compute_motifs(delaunay_in::TopologicalNetwork,r=-1)
    if delaunay_in.dim == 2
        if r < 0; r = 2; end
        idx, tvec  = weinberg2D(delaunay_in,r)

    elseif delaunay_in.dim == 3
        if r < 0; r = 1; end
        idx, tvec = simplicial_3D(delaunay_in)
    end
    return MotifArray(idx,tvec,delaunay_in.dim,r)
end

function avg_motif(motif::Vararg{MotifDist,N}) where {N}
        tvecs = unique(vcat([collect(keys(m.map)) for m in motif]...))
        count = zeros(Int64,length(tvecs))
        for i = 1:length(tvecs)
            for m in motif
                if haskey(m.map,tvecs[i])
                    count[i] = count[i] + m.map[tvecs[i]]
                end
            end
        end
        return MotifDist(Dict(tvecs .=> count),motif[1].dim,motif[1].r)
end
function avg_motif(motif::Vararg{MotifArray,N}) where {N}
        return MotifDist(countmap(vcat([m.tvec for m in motif]...)),motif[1].dim,motif[1].r)
end
################################################################################
# 2D Labeling code
################################################################################

function shape_find(index_points,vert,regions,ridge_pts,ridge_vert,ind)
    # This function takes in an index and returns the voronoi cell associated
    # with that index, and a graph of the vertices and edges of the polytope
    kk = index_points[ind]
    if any(regions[kk] .== -1)
        error("error in shape_find")
    end

    poly_pts = vert[regions[kk].+1,:] # the actual points of the polytope
    faces = ridge_vert[ (ridge_pts[:,1] .== ind) .| (ridge_pts[:,2] .== ind)]
    point_lookup = Dict(regions[kk].+1 .=> 1:length(regions[kk]) )

    g = SimpleGraph(length(regions[kk]))

    for s in 1:length(faces)
        face = faces[s].+1
        for i = 1:(length(face)-1)
            add_edge!(g, point_lookup[face[i]], point_lookup[face[i+1]])
        end
        add_edge!(g,point_lookup[face[end]],point_lookup[face[1]])
    end
    face = faces[1].+1
    face_points = [point_lookup[face[i]] for i in 1:length(face)]

    return g, poly_pts, face_points
end



function right_label(Points,g,ind)
    # This function orders the neighbouring points around a vertex, ind, according to their
    # counter clockwise orientation, where the the out of plane vector is defined
    # to be from the origin to the vertex. Shift origin so it lies within the shape
    Points = Points - repeat(mean(Points,dims=1),size(Points,1),1)
    renorm = x-> x/norm(x)
    e1 = renorm(Points[ind,:])
    nbs = neighbors(g,ind)
    a = Points[nbs[1],:]
    e2 = renorm(a - dot(a,e1)*e1) #/norm(a - dot(a,e1)*e1)
    e3 = cross(e1,e2)
    plane_angle = [angle(dot(Points[n,:],e2) + dot(Points[n,:],e3)*im) for n in nbs[2:end]]
    plane_angle = [(p < 0) ? p + 2π : p for p in plane_angle]
    order =  nbs[ vcat([1], (sortperm(plane_angle).+1))]
    return order
end

function weinberg_vect(g_aug,e,order_mat; mirror=false)
    # This function is the core part of the Weinberg algorithm. It takes the
    # augmented graph (with forward and backward edges) and a starter edge and
    # finds the Weinberg code. The mirror option either picks the right most (default) or
    # left most available edge (or right most in the mirrored representation)
    vertex_to_code = Dict(src(e) => 1 ) # convert vertex numbers to their code numbering
    code = [1]
    current_vertex = dst(e)
    previous_vertex = src(e)
    while ne(g_aug) > 0

        if haskey(vertex_to_code, current_vertex)
            push!(code,vertex_to_code[current_vertex])
            if has_edge(g_aug,current_vertex,previous_vertex)
                rem_edge!(g_aug, previous_vertex,current_vertex)
                previous_vertex,current_vertex = current_vertex,previous_vertex
                continue
            end
        else
            vertex_to_code[current_vertex] = maximum(code)+1
            push!(code,vertex_to_code[current_vertex])
        end


        edge_choices = order_mat[current_vertex]

        ind = findfirst(x->x==previous_vertex,edge_choices)
        if mirror
            edge_choices = circshift(edge_choices, 1-ind)
            next_vertex_ind = findlast(ed -> has_edge(g_aug, current_vertex,ed),edge_choices)
        else
            edge_choices = circshift(edge_choices, -ind)
            next_vertex_ind = findfirst(ed -> has_edge(g_aug, current_vertex,ed),edge_choices)
        end
         if isnothing(next_vertex_ind)
             return code
         end

        rem_edge!(g_aug, previous_vertex,current_vertex)
        previous_vertex,current_vertex = current_vertex,edge_choices[next_vertex_ind]
    end
    error("Error occured in Weinberg computation")
end

function augmented_graph(g)
    g_aug = SimpleDiGraph(size(g)[1])
    for e in edges(g)
        add_edge!(g_aug, src(e), dst(e))
        add_edge!(g_aug, dst(e), src(e))
    end
    return g_aug
end

function weinberg_find!(code_tot,S_tot,kk,g,order_mat,cent_node = -1)
    g_aug = augmented_graph(g)
    vecs = Array{Int64}[]
    for e in edges(g_aug)
        if (cent_node < 0) || (src(e) == cent_node)
            append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=true)])
            append!(vecs,[weinberg_vect(augmented_graph(g),e,order_mat;mirror=false)])
        end
    end
    sort!(vecs)
    S_tot[kk] = Int64(length(vecs)/length(unique(vecs))) # size of the symmetry group
    code_tot[kk] =  copy(vecs[1])

end


function edge_neighbors(g,edge_index)
    e_neighbors = Array{Int64}(undef,0)
    for e in edge_index
        append!(e_neighbors,neighbors(g,e))
    end
    return unique!(e_neighbors)
end


function weinberg2D_core(g,N,idx,r=2)
    code_tot = Vector{Array{Int64}}(undef,N)
    S_tot = Array{Int64}(undef,N)

    for k in idx
        nbh  = neighborhood(g,k,r)
        g_ego, vmap = induced_subgraph(g,nbh)
        vmap_inv = Dict(vmap[k] => k for k in 1:length(vmap))

        x,y,fixed_vecs = tutte_embedding(g_ego)
        order_local = order_mat_find(g_ego,x,y)
        
        weinberg_find!(code_tot,S_tot,k,g_ego,order_local,vmap_inv[k])
    end
    return code_tot[idx], S_tot[idx]
end

function motif_size_find(Pos,r=2;periodic=true)
    Positions = deepcopy(Pos)
    N = length(Positions)
    index_ref = [x for x in 1:N]
    if periodic; periodic_extend!(Positions,index_ref); end
    # All other operations are performed both periodic and non-periodic systems,
    # but for non-periodic systems should have no effect
    _, simplices, _, _ = Delaunay_find(Positions)
    g_full = graph_construct(simplices,length(Positions))
    motif_lens = [length(neighborhood(g_full,k,r)) for k in 1:N]
    return motif_lens
end


function weinberg2D(delaunay_in,r)

    simplices = delaunay_in.simplices
    g = graph_construct(simplices,maximum(simplices))
    edge_index = setdiff(1:nv(g), delaunay_in.not_edge)

    periodic = delaunay_in.periodic

    N = periodic ? delaunay_in.original_vertex_number : nv(g)

    idx = Vector(1:N)
    if !periodic
        for i = 1:r-1
            edge_index = edge_neighbors(g,edge_index)
        end
        idx = setdiff(1:nv(g), edge_index)
    end

    code_tot,S_tot = weinberg2D_core(g,N,idx,r)

    return idx, code_tot,S_tot
    
end

################################################################################
# 3D Labeling code
################################################################################

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
        num_labeled = 4#length(unique(nbhd[num_zero .== -1,:]))
        fill_labels(which_to_label,nbhd,num_zero,sim_code,N,Topological_vec,idx_loc,num_labeled)

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
                for j = 1:length(p)
                    add_labeled!(nbhd,which_to_label_cp,idx_loc[lab_next[p[j]]],
                                length(already_labeled)+j,num_zero_cp,sim_code_cp)
                end
                            num_zero_cp[s] = -1

                Topological_vec_cp[findlast(x-> x< typemax(Int64),Topological_vec_cp) .+ 1 ] = sim_code_cp[s]

                num_labeled = -1
                fill_labels(which_to_label_cp,nbhd,num_zero_cp,sim_code_cp,N,Topological_vec_cp,idx_loc,num_labeled)

                push!(tvec_tot,copy(Topological_vec_cp))
            end
        end
    end

    Topological_vec .= minimum(tvec_tot)
end


function fill_labels(to_label::Array{Int64,1}, Orig::Array{Int64,2},
            num_zero::Array{Int64,1},s_code::Array{Int64,1},N::Int64,
            Topological_vec,idx_loc,num_labeled)
    # This is the core function of the labeling algorithm. It assumes that the
    # first simplex has been chosen and labeled. The rest of the labeling now
    # procedes deterministically, but will terminate if the labeling is not less
    # than the current value of Topological_vec
    N_sim_labeled = sum(num_zero .== -1)

    flag = 0
    if num_labeled > 0
        n_l  = num_labeled + 1#length(unique(Orig[num_zero .== -1,:])) + 1
    else
        n_l = length(unique(Orig[num_zero .== -1,:])) + 1
    end
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
    nbhd = similar(k_nbhd)
    return topological_vec_mem_save(nbhd,k_nbhd,central_vertex;r=1)
end

function topological_vec_mem_save(nbhd,k_nbhd,central_vertex;r=1)
    # This function takes a neighborhood of simplices around a central vertex
    # and gives it a cannonical labeling.

    # First make a copy of k_nbhd and relabel the vertices 1:N
    #nbhd = copy(k_nbhd)
    nbhd_unique = unique(k_nbhd)
    N = length(nbhd_unique)
    if N > 62
        println("Too many to store in int64")
        println(N)
        println(nbhd_unique)
        return []
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
                num_labeled = 4
                optimal = fill_labels(which_to_label,nbhd,num_zero,sim_code,N,Topological_vec,idx_loc,num_labeled)

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
    if (maximum(Topological_vec) < typemax(Int64)) && (length(Topological_vec) != length(unique(Topological_vec)))
        println("tvec = ", Topological_vec)
        println("simplex = ",nbhd)
        println("central_vertex = ", central_vertex)
        error()
    end
    #if maximum(Topological_vec) == typemax(Int64)
        #println("tvec = ", Topological_vec)
        #println("Simplex = ", nbhd)
        #error()
    #end
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



function simplicial_3D(delaunay_in)

    simplices = delaunay_in.simplices
    not_edge = delaunay_in.not_edge
    periodic = delaunay_in.periodic

    if periodic
        N = delaunay_in.original_vertex_number
        not_edge = not_edge[not_edge .<= N] # don't compute motif for periodic copies
    end

    tvec_tot = Array{Int64}[]
    for i in 1:length(not_edge)
        k_nbhd = find_nbhd(simplices,not_edge[i])
        if length(k_nbhd) > 0
            push!(tvec_tot,topological_vec(k_nbhd,not_edge[i]))
        else
            not_edge[i] = -1
        end
    end

    not_edge = not_edge[not_edge.>0]
    return not_edge, tvec_tot
end
