using LightGraphs
using DataFrames, CSV
using Base.Threads

function three_match(A,k,l)
    # finds if rows k and l of A share a face
    count = 0
    #@inbounds for i = 1:4, j = (1+count):4
    for i = 1:4, j = (1+count):4
        if A[k,i] == A[l,j]
            count += 1
        end
    end
    return count == 3
end

# First find all edges corresponding to 3 simplices that can be flipped

function find_3_simplex(simplex)
    to_return = Vector{Int64}[]
    pairs = falses(size(simplex,1))
     for i = 1:( size(simplex,1) - 2)
        pairs .= false
        # First find each tetrahedron that shares a face with tetrahedron i
        for j = (i+1): size(simplex,1)
            pairs[j] = three_match(simplex,i,j)
        end
        # Loop through the pairs and see if a third one is paired with both
        for j = (i+1): (size(simplex,1)-1)
            if pairs[j]
                for k = (j+1) : size(simplex,1)
                    if pairs[k] && three_match(simplex,j,k)
                            push!(to_return,[i;j;k])
                    end
                end
            end
        end
    end
    return to_return
end


function flip_partition!(count1,count2,simplex,ijk)
    # Finds the common edge of the three simplices that are to be flipped.
    #@inbounds for i in 1:4
    for i in 1:4
        scan = false
        for i2 = 2:3, j2 = 1:4
            if simplex[ijk[1],i] == simplex[ijk[i2],j2]
                scan = !scan
            end
        end
        if scan
            push!(count1,simplex[ijk[1],i])
        else
            push!(count2,simplex[ijk[1],i])
        end
    end

    #@inbounds for i = 1:4
    for i = 1:4
        scan = true
        for j = 1:4
            if simplex[ijk[2],i] == simplex[ijk[1],j]
                scan = false
            end
        end
        if scan
            push!(count1,simplex[ijk[2],i])
        end
    end
end


function perform_3_flip(simplex,ijk,central_vertex; r = 1)

    keep = trues(size(simplex,1))
    for i in ijk
        keep[i] = false
    end
    new_simplex = simplex[keep,:]

    # find which indices will be affected by the flip
    count1 = Int64[]
    count2 = Int64[]
    flip_partition!(count1,count2,simplex,ijk)
    row1 = zeros(Int64,4)
    row2 = zeros(Int64,4)
    for i = 1:3
        row1[i] = count1[i]
        row2[i] = count1[i]
    end
    row1[4] = count2[1]
    row2[4] = count2[2]
    sort!(row1); sort!(row2)

    if r > 1
        error("TODO")
        # need to do proper checks to see if the new simplex is still in the
        # local network of the origin. For r=1 this is easy and is done below
    end
    if central_vertex == count2[1] # TODO search and remove points not in new local simplicial complex
        new_simplex = [new_simplex; row1']
    elseif central_vertex == count2[2]
        new_simplex = [new_simplex; row1']
    else
        new_simplex = [new_simplex; row1';row2']
    end

    tvec = topological_vec(new_simplex,central_vertex,r=r)
    for j = 2:length(tvec)
        if tvec[j-1] == tvec[j]
            # Can't flip as there is a parallel tetrahedron
            return Int64[]
        end
    end
    return tvec

end


function face_find(simplex,row,face)
    for i = 1:size(simplex,1)
        count = 0
        if i != row
            for j = 1:4, f in face
                if simplex[i,j] == f
                    count += 1
                end
            end
        end
        if count == 3
            return i
        end
    end
    return 0
end

function find_2_external(simplex)
    # Finds all pairs of simplces that are neighbours and could be neighbours
    # with some external simplex, hence could be flipped
    faces = [[1;2;3],[1;2;4],[1;3;4],[2;3;4]]
    nbhd = zeros(Int64,4)
    to_return = Vector{Int64}[]

    #@inbounds for i = 1:size(simplex,1)
    for i = 1:size(simplex,1)
        external = false
        for j = 1:4
            nbhd[j] = face_find(simplex,i,simplex[i,faces[j]])
            if nbhd[j] == 0
                external = true
            end
        end
        if external
            for j = 1:4
                if nbhd[j] > i
                    push!(to_return,[i;nbhd[j]])
                end
            end
        end
    end
    return to_return
end

function count_external_faces(simplex,i)
    # counts the number of external faces that a particular simplex has
    faces = [[1;2;3],[1;2;4],[1;3;4],[2;3;4]]
    count = 0
    for j = 1:4
        count += face_find(simplex,i,simplex[i,faces[j]]) == 0
    end
    return count
end

function find_face_external!(face_array,simplex,i)
    # counts the number of external faces that a particular simplex has
    faces = [[1;2;3],[1;2;4],[1;3;4],[2;3;4]]
    for j = 1:4
        if face_find(simplex,i,simplex[i,faces[j]]) == 0
            push!(face_array,simplex[i,faces[j]])
        end
    end
end

function find_2A(simplex,central_vertex, r = 1)
    # Finds all simplices that could be deleted in a flip of type 2A
    # Need only find simplices with at least 2 exposed faces and not containing
    # the central vertex (for r =1); these are then eligible to be flipped
    to_flip = Int64[]

    for i = 1:size(simplex,1)
        proceed = true
        if r != 1
            error("TODO code for r > 1")
            # for r = 1 just check the simplices involved aren't the central vertex
            # for r > 1 would need to check that they are all distance r away
            # from the central vertex
        end
        for j = 1:4
            if central_vertex == simplex[i,j]
                proceed = false
                continue
            end
        end
        if proceed
            if count_external_faces(simplex,i) > 1
                push!(to_flip,i)
            end
        end
    end
    return to_flip
end


function perform_flip_2A(simplex,i,central_vertex;r=1)
    # given a simplex with 2 external faces, we perform the flip of type 2A
    # and return the new topological vector
    new_simplex = copy(simplex)
    new_simplex = new_simplex[ vcat(1:(i-1),(i+1):size(simplex,1)),:]
    return topological_vec(new_simplex,central_vertex,r=r)
end

function find_2A_rev(simplex, r = 1)
    # Does the reverse of a 2A flip. Needed to better calculate the flip graph.
    # Find all pairs of exposed faces joined at an edge, and add a new simplex
    # with those 4 vertices
    to_flip = Vector{Int64}[]

    # First find all external faces
    external_faces = Vector{Int64}[]
    for i = 1:size(simplex,1)
        find_face_external!(external_faces,simplex,i)
    end

    # If 2 external faces have a common edge, then can create new simplex
    for i = 1: length(external_faces)-1
        for j = (i+1): length(external_faces)
            if faces_share_edge(external_faces,i,j)
                push!(to_flip, unique(vcat(external_faces[i],
                external_faces[j])))
            end
        end
    end
    return to_flip
end

function perform_flip_2A_rev(simplex,to_create,central_vertex;r=1)
    # given a simplex with 2 external faces, we perform the flip of type 2A in reverse
    # and return the new topological vector
    new_simplex = [simplex;to_create']
    return topological_vec(new_simplex,central_vertex,r=r)
end

function faces_share_edge(face_array,i,j)
    # Check if the faces face_array[i], face_array[j] share an edge i.e. 2 vertices
    count = 0
    for l = 1:3, s = 1:3
        if face_array[i][l] == face_array[j][s]
            count +=1
            if count == 2
                return true
            end
        end
    end
    return false
end


function find_2B(simplex)
    # Finds all simplices that could be created in a flip of type 2B
    # Need to find all vertices connected to exactly 3 exposed faces
    to_flip = Vector{Int64}[]

    # First find all external faces
    external_faces = Vector{Int64}[]
    for i = 1:size(simplex,1)
        find_face_external!(external_faces,simplex,i)
    end

    # If 3 external faces pairwise all have a common edge, then can create new simplex
    for i = 1: length(external_faces)-2
        for j = (i+1): length(external_faces)-1
            if !faces_share_edge(external_faces,i,j)
                continue
            end
            for k = (j+1): length(external_faces)
                if faces_share_edge(external_faces,i,k) &&
                     faces_share_edge(external_faces,j,k)
                     M = unique(vcat(external_faces[i],
                     external_faces[j],external_faces[k]))
                     if length(M) == 4 # exclude degenerate cases
                         push!(to_flip, M)

                     end
                end
            end
        end
    end
    return to_flip
end

function perform_flip_2B(simplex,to_create,central_vertex;r=1)
    # given a simplex with 2 external faces, we perform the flip of type 2B
    # and return the new topological vector
    new_simplex = [simplex; to_create']
    return topological_vec(new_simplex,central_vertex,r=r)
end
#=
function find_idx_counts!(idx_seen,n_count,simplex,ij)
    # fills out vectors idx_seen and n_count with the vertices seen and the
    # number of times they were seen
    #@inbounds for i = 1:4
    for i = 1:4
        idx_seen[i] = simplex[ij[1],i]
        n_count[i] = 1
    end
    #@inbounds for i = 1:4
    for i = 1:4
        seen = false
        for j=1:4
            if idx_seen[j] == simplex[ij[2],i]
                n_count[j] += 1
                seen = true
            end
        end
        if !seen
            idx_seen[5] = simplex[ij[2],i]
            n_count[5] = 1
        end
    end
end

function perform_2_plus_external_flip(simplex,ij,central_vertex;r=1)
    # given 2 external simplices that can be flipped, we find the post flip
    # simplices and return the new topological vector
    idx_seen = zeros(Int64,5)
    n_count = zeros(Int64,5)
    find_idx_counts!(idx_seen,n_count,simplex,ij)

    faces = [[2;3;4],[1;3;4],[1;2;4],[1;2;3]]
    fixed_v1 = Int64[]
    fixed_v2 = Int64[]
    for j = 1:4
        if face_find(simplex,ij[1],simplex[ij[1],faces[j]]) == 0
            push!(fixed_v1,simplex[ij[1], j])
        end
        if face_find(simplex,ij[2],simplex[ij[2],faces[j]]) == 0
            push!(fixed_v2,simplex[ij[2],j])
        end
    end
    fixed_v = intersect(fixed_v1,fixed_v2)[1]
    row1 = zeros(Int64,4)
    count1 = 2
    row2 = zeros(Int64,4)
    count2 = 2

    row1[1] = fixed_v
    row2[1] = fixed_v

    for i = 1:5
        if n_count[i] == 1
            row1[count1] = idx_seen[i]
            count1 += 1
            row2[count2] = idx_seen[i]
            count2 += 1
        elseif (n_count[i] == 2) && (idx_seen[i] != fixed_v)
            if count1 == count2
                row1[count1] = idx_seen[i]
                count1 += 1
            else
                row2[count2] = idx_seen[i]
                count2 += 1
            end
        end
    end

    sort!(row1)
    sort!(row2)
    new_simplex = copy(simplex)
    for i = 1:4
        new_simplex[ij[1],i] = row1[i]
        new_simplex[ij[2],i] = row2[i]
    end
    return topological_vec(new_simplex,central_vertex,r=r)
end
=#
function find_all_neighbors(attempt,simplex,central_vertex;r=1)
    # Given a local simplicial complex, this function finds all others that
    # are 1 flip away, trying type 1A,1B,2A,2B, flips along the way.
    tvec_nhbd = []

    to_ret = find_3_simplex(simplex)
    for t in to_ret
        push!(tvec_nhbd,perform_3_flip(simplex,t,central_vertex))
    end

    to_flip = find_2A(simplex,central_vertex)
    for s in to_flip
        push!(tvec_nhbd,perform_flip_2A(simplex,s,central_vertex))
    end

    if attempt == 1
        # Only do the reverse caclulation for the initial vertices to save time
        to_flip = find_2A_rev(simplex)
        for kk in 1:length(to_flip)
            push!(tvec_nhbd,perform_flip_2A_rev(simplex,to_flip[kk],central_vertex;r=r))
        end
    end

    to_flip = find_2B(simplex)
    for kk in 1:length(to_flip)
        push!(tvec_nhbd,perform_flip_2B(simplex,to_flip[kk],central_vertex;r=r))
    end
    unique!(tvec_nhbd)
    return tvec_nhbd
end

function find_flip_graph3D(tvec_tot)
    g = SimpleGraph(length(tvec_tot))
    code_to_idx_orig = Dict(tvec_tot[i] => i for i in 1:length(tvec_tot))
    code_to_idx = Dict(tvec_tot[i] => i for i in 1:length(tvec_tot))
    println(nv(g))
    tvec_new = []
    splock = SpinLock()
    Threads.@threads for i = 1:length(tvec_tot)
    #for i = 1:length(tvec_tot)
        tvec = tvec_tot[i]
        simplex = zeros(Int64,length(tvec),4)
        t_vec_to_simplex!(tvec,simplex)
        central_vertex = 1 # By convention

        tvec_nhbd = find_all_neighbors(1,simplex,central_vertex,r=1)

       lock(splock)
        for t in tvec_nhbd
            if haskey(code_to_idx,t)
                add_edge!(g,code_to_idx[t],code_to_idx[tvec])
            else
                code_to_idx[t] = nv(g) + 1
                add_vertex!(g)
                add_edge!(g,code_to_idx[t],code_to_idx[tvec])
                push!(tvec_new,t)
            end
        end
        unlock(splock)

    end
    println(ne(g))

    # Now add the edges for the newly added vecs
    splock = SpinLock()
    Threads.@threads for i = 1:length(tvec_new)
        tvec = tvec_new[i]
        simplex = zeros(Int64,length(tvec),4)
        t_vec_to_simplex!(tvec,simplex)
        central_vertex = 1 # By convention

        tvec_nhbd = find_all_neighbors(2,simplex,central_vertex,r=1)

       lock(splock)
        for t in tvec_nhbd
            # This time ignore if their neighbours haven't been seen before
            if haskey(code_to_idx,t)
                add_edge!(g,code_to_idx[t],code_to_idx[tvec])
            end
        end
        unlock(splock)

    end

    # remove vertices that won't effect the diffusion calculation
    for v = nv(g) : -1 : (length(tvec_tot) +1)
        if length(neighbors(g,v)) < 2
            rem_vertex!(g,v)
        end
    end

    # To thin graph even further, remove vertices that are not central enough
    rank = pagerank(g)
    to_keep = rank .> 1.4*mean(rank)
    to_keep[1:length(tvec_tot)] .= true
    g_new, d = induced_subgraph(g,[i for i in 1:nv(g)][to_keep])


    # Use this testing code to see what fraction of the motifs are accounted for
    # in the largest connected component of the new flip graph.
    #=
        # Load the number of unique motifs and how often they occur
        str_arr = filter(x->occursin("avg.txt",x), readdir(Dir_out))
        weight = amalg2([readin(Dir_out*s,0) for s in str_arr])
        weight = [w[2] for w in weight]
        core_v = length(weight)

        # Verify that most of the motifs are in the connected component of the new flip graph
        g_comp = connected_components(g_new)
        v_in_core = intersect([i for i = 1:core_v],g_comp[1])
        sum(weight[v_in_core])/sum(weight)
    =#
    println(nv(g_new))
    println(ne(g_new))

    return g_new
end

function compute_flip_graph3D(code_amalg,save_str)
    tvec_tot = []
    for i = 1:size(code_amalg,1)
        w = code_amalg[i][1]
        w_num = Meta.parse(w)
        push!(tvec_tot,Int.(w_num.args))
    end
    g = find_flip_graph3D(tvec_tot)
    savegraph( save_str*".lgz", g)
    df = DataFrame(codes = [code_amalg[k][1] for k in 1:size(code_amalg,1)],
        index = [k for k in 1:size(code_amalg,1)])
    CSV.write(save_str*".txt",  df)
end
#=
tvec_tot = []
for i = 1:1000
    k_nbhd = find_nbhd(simplices,i)
    push!(tvec_tot,topological_vec(k_nbhd,i))
end
    unique!(tvec_tot)



simplex = zeros(Int64,length(tvec_tot[1]),4)
t_vec_to_simplex!(tvec_tot[1],simplex)
@profiler g = find_flip_graph3D(tvec_tot)
=#
