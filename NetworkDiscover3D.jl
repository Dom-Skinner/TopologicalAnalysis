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
    if central_vertex == count2[1]
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


function find_flip_graph3D(tvec_tot)
    g = SimpleGraph(length(tvec_tot))
    code_to_idx = Dict(tvec_tot[i] => i for i in 1:length(tvec_tot))

    splock = SpinLock()
    Threads.@threads for i = 1:length(tvec_tot)
        tvec = tvec_tot[i]
        simplex = zeros(Int64,length(tvec),4)
        t_vec_to_simplex!(tvec,simplex)
        central_vertex = 1 # By convention

        tvec_nhbd = []

        to_ret = find_3_simplex(simplex)
        for t in to_ret
            push!(tvec_nhbd,perform_3_flip(simplex,t,central_vertex))
        end

        to_ret = find_2_external(simplex)
        for t in to_ret
            push!(tvec_nhbd,perform_2_plus_external_flip(simplex,t,central_vertex))
        end
        unique!(tvec_nhbd)

        lock(splock)
        for t in tvec_nhbd
            if haskey(code_to_idx,t)
                add_edge!(g,code_to_idx[t],code_to_idx[tvec])
            end
        end
        unlock(splock)

    end

    #TODO? Could include the unseen vectors into the flip graph and then find
    # their neighbours to get a larger connected component.
    return g
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
