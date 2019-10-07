using LightGraphs

#=
function graph_construct(Positions)
  p, simplices, neighbors, edge_index = Delaunay_find(Positions)
  g  = SimpleGraph(length(p[:,1]))

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
=#

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
