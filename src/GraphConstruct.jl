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
    return [length(neighbors(g,i)) for i = 1:nv(g)]
end


function get_peripheral_cycle(g,cent_node,r=2)
    # Find vertices r-1 edges from central node
    outer_vertices = setdiff(1:nv(g),neighborhood(g,cent_node,r-1))
    g_periph,vmap = induced_subgraph(g,outer_vertices)

    
    #####################
    # This block of code is not necessary but saves time for big graph_construct
    # It identifies edges that are clearly not part of the largest peripheral cycle
    # and removes them. Saves time in the expensive cycle computations.
    for i = 1:nv(g_periph)
        nn  = neighbors(g_periph,i)
        if (nv(g_periph) > 3) && (length(nn) == 2) && has_edge(g_periph,nn[1],nn[2])
            nn  = copy(neighbors(g_periph,i))
            rem_edge!(g_periph,nn[1],i)
            rem_edge!(g_periph,nn[2],i)
        end
    end
    ####################
    
    # Find cycles in the subgraph of outer edges
    
    cycles = simplecycles_limited_length(g_periph, length(outer_vertices))  
    cycles = filter(x->length(x)>2,cycles)
    # only keep cycles that have no chords
    cycles = filter(x->ne(induced_subgraph(g_periph,x)[1]) == length(x),cycles)
    
    c_special = vmap[cycles[argmax(length.(cycles))]]
    
    # Check that it is really a peripheral cycle
    test_con,_ = induced_subgraph(g,setdiff(1:nv(g),c_special))
    @assert is_connected(test_con)

    return c_special
end

function tutte_embedding(g,cent_node,r=2)
    # This function creates a Tutte embedding of the graph.
 
    
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

    
    # arrange the points making up the outer edges around a circle
    c_special = get_peripheral_cycle(g,cent_node,r)
    n_poly = length(c_special)
    
    x = zeros(nv(g)) ;      y = zeros(nv(g))
    x[c_special] .= cos.(2*pi/n_poly .* (0:n_poly-1))
    y[c_special] .= sin.(2*pi/n_poly .* (0:n_poly-1))

    bx = A*x ; by = A*y
    keep_idx = [k ∉ c_special for k in 1:nv(g)]
    bx = bx[keep_idx] ;  by = by[keep_idx]
    A = A[keep_idx,keep_idx]

    x[keep_idx] = - A \ bx
    y[keep_idx] = - A \ by
    return x,y,c_special
end



function order_mat_find(g,x,y)
    # Finds the order mat given the graph and it's embedding
    order_mat = Vector{Array{Int64}}(undef,nv(g))
    for kk = 1:nv(g)
        N_list = neighbors(g,kk)
        theta = [atan(y[s]-y[kk], x[s]-x[kk]) for s in N_list]
        tst=filter(x->x>0,abs.(theta .- theta'))
        if (length(tst) > 0) && (minimum(tst) < 1e-15)
            println("WARNING: Insufficient precision in graph embedding.")
        end
        order_mat[kk] = N_list[sortperm(theta)]
    end
    return order_mat
end
