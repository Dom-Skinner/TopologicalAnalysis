using LightGraphs

function graph_construct(simplices,N)
    g  = SimpleGraph(N)
  for i = 1:length(simplices[:,1])
    add_edge!(g, simplices[i,1], simplices[i,2])
    add_edge!(g, simplices[i,1], simplices[i,3])
    add_edge!(g, simplices[i,2], simplices[i,3])
    if length(simplices[i,:]) == 4
      add_edge!(g, simplices[i,1], simplices[i,4])
      add_edge!(g, simplices[i,2], simplices[i,4])
      add_edge!(g, simplices[i,3], simplices[i,4])
    end
  end
  return g
end


function vertex_num(g)
    N = nv(g)
    vertex_info = zeros(Int64,N)
    for i = 1:N
      vertex_info[i] = length(neighbors(g,i))
    end
    return vertex_info
end


function tutte_embedding(g)
    # This function creates a Tutte embedding of the graph.
    fixed_vecs = []
    for e in edges(g)
        if length(intersect(neighbors(g,src(e)),neighbors(g,dst(e)))) == 1
            push!(fixed_vecs,src(e))
            push!(fixed_vecs,dst(e))
            push!(fixed_vecs,intersect(neighbors(g,src(e)),neighbors(g,dst(e)))[1])
            break
        end
    end
    if length(fixed_vecs) != 3
        for e in edges(g)
            push!(fixed_vecs,src(e))
            push!(fixed_vecs,dst(e))
            push!(fixed_vecs,intersect(neighbors(g,src(e)),neighbors(g,dst(e)))[1])
            break
        end
        #println("warning: Degenerate case in Tutte embedding")
        # degenerate case being when there are no edges attached to only a single triangle
    end
    A = float(adjacency_matrix(g))
    D = [sum(A[i,:]) for i in 1:size(A,1)]
    for i = 1:size(A,1)
       for j in 1:size(A,2)
           if A[i,j] != 0
               A[i,j] /= D[i]
           end
       end
    end
    for i = 1:size(A,1)
        A[i,i] -= 1
    end

    x = zeros(nv(g)) ;      y = zeros(nv(g))
    x[fixed_vecs[1]] = 0.;  y[fixed_vecs[1]] = 0.;
    x[fixed_vecs[2]] = 0.;  y[fixed_vecs[2]] = 1.;
    x[fixed_vecs[3]] = 1.;  y[fixed_vecs[3]] = 0.;
    bx = A*x ; by = A*y
    keep_idx = [k ∉ fixed_vecs for k in 1:nv(g)]
    bx = bx[keep_idx] ;  by = by[keep_idx]
    A = A[keep_idx,keep_idx]

    x[keep_idx] = - A \ bx
    y[keep_idx] = - A \ by
    return x,y,fixed_vecs
end


function order_mat_find(g,x,y)
    # Finds the order mat given the graph and it's embedding
    order_mat = Vector{Array{Int64}}(undef,nv(g))
    for kk = 1:nv(g)
        N_list = neighbors(g,kk)
        theta = [atan(y[s]-y[kk], x[s]-x[kk]) for s in N_list]
        order_mat[kk] = N_list[sortperm(theta)]
    end
    return order_mat
end
