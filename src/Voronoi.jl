#https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html
#https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html
#https://graphics.stanford.edu/courses/cs268-11-spring/handouts/AlphaShapes/as_fisher.pdf

using PyCall
using LinearAlgebra: norm, det
using Statistics: median, std
using CSV, DataFrames
const scipy_qhull = PyNULL()


#############################################################################
struct TopologicalNetwork
           simplices::Array
           not_edge::Array
           dim::Int
           periodic::Bool
           alpha::Number
           edge_keep::Bool
           tolerance::Number
           original_vertex_number::Int
           clusters::Union{Array, Missing}
end

function TopologicalNetwork(simplices::Array, not_edge::Array, dim::Int,
        periodic::Bool, alpha::Number, edge_keep::Bool, tolerance::Number,
        original_vertex_number::Int)
    return TopologicalNetwork(simplices, not_edge, dim, periodic, alpha,
        edge_keep,tolerance,original_vertex_number,missing)
end

function TopologicalNetwork(tn::TopologicalNetwork, clusters::Array)
    return TopologicalNetwork(tn.simplices, tn.not_edge, tn.dim, tn.periodic, tn.alpha,
        tn.edge_keep,tn.tolerance,tn.original_vertex_number,clusters)
end


function find_delaunay_network(path_to_csv_in::String; periodic=false, alpha = 0,
                                tol=0.6, edge_keep=false)

    dat_in = Matrix(CSV.read(path_to_csv_in,DataFrame))
    Positions = unique([dat_in[k,:] for k in 1:size(dat_in,1)])

    if length(Positions) != size(dat_in,1)
        println("Warning: Input file contains duplicate points.\n Recommended to rerun with duplicates removed")
    end

    return find_delaunay_network_core(Positions, periodic, alpha, tol, edge_keep)

end
function find_delaunay_network(Positions::Array; periodic=false, alpha = 0,
                                tol=0.6, edge_keep=false)
    return find_delaunay_network_core(Positions, periodic, alpha, tol, edge_keep)
end
#############################################################################

function __init__()
    copy!(scipy_qhull, pyimport_conda("scipy.spatial.qhull", "scipy"))
end

function circumradius(e1,e2,e3)
    # see http://mathworld.wolfram.com/Circumradius.html
    a = sqrt(e1[1]^2+e1[2]^2)
    b = sqrt(e2[1]^2+e2[2]^2)
    c = sqrt(e3[1]^2+e3[2]^2)
    return a*b*c/sqrt( (a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))
end

function sidelengths(e1,e2,e3)
    # see http://mathworld.wolfram.com/Circumradius.html
    a = sqrt(e1[1]^2+e1[2]^2)
    b = sqrt(e2[1]^2+e2[2]^2)
    c = sqrt(e3[1]^2+e3[2]^2)
    return [a,b,c]
end

function circumradius3D(e1,e2,e3,e4)
    # see http://mathworld.wolfram.com/Circumsphere.html
    e_n = [norm(e1)^2 ; norm(e2)^2; norm(e3)^2; norm(e4)^2]
    e_mat = [e1' 1 ;e2' 1 ;e3' 1 ;e4' 1]
    e_mat = [e_n e_mat]
    a = det(e_mat[:,2:5] )
    Dx = det(e_mat[:,[1; 3:5]])
    Dy = det(e_mat[:,[1:2; 4:5]])
    Dz = det(e_mat[:,[1:3; 5]])
    c = det(e_mat[:,1:4])
    r = Dx^2+Dy^2+Dz^2-4*a*c
    if r >= 0 && abs(a) > 0
    	return sqrt(r)/(2*abs(a))
    else
    	return Inf # infinite radius
    end
end

# function alpha_shape2D(simplices,p)
#     α_val = zeros(length(simplices[:,1]))
#     for s = 1:length(simplices[:,1])
#         e1 = p[simplices[s,1],:] - p[simplices[s,2],:] #side lengths in 2D of triangle
#         e2 = p[simplices[s,1],:] - p[simplices[s,3],:]
#         e3 = p[simplices[s,2],:] - p[simplices[s,3],:]
#         α_val[s] = circumradius(e1,e2,e3)
#     end
#     return α_val
# end

function alpha_shape2D(simplices,p) #modified by Elise
    α_val = zeros(length(simplices[:,1]),3)
    for s = 1:length(simplices[:,1])
        e1 = p[simplices[s,1],:] - p[simplices[s,2],:] #side lengths in 2D of triangle
        e2 = p[simplices[s,1],:] - p[simplices[s,3],:]
        e3 = p[simplices[s,2],:] - p[simplices[s,3],:]
        α_val[s,:] = sidelengths(e1,e2,e3)
    end
    return α_val #Nx3 array of edge lengths
end

function alpha_shape3D(simplices,p)
    α_val = zeros(length(simplices[:,1]))
    for s = 1:length(simplices[:,1])
        e1 = p[simplices[s,1],:]
        e2 = p[simplices[s,2],:]
        e3 = p[simplices[s,3],:]
        e4 = p[simplices[s,4],:]
        α_val[s] = circumradius3D(e1,e2,e3,e4)
    end
    return α_val
end

function alpha_shape_eval(simplices,p)
    # This function performs an alpha shape analysis. It computes the radius of
    # circumcircle/sphere for each simplex.
    if size(p,2) == 3
        return alpha_shape3D(simplices,p)
    elseif size(p,2) == 2
        return alpha_shape2D(simplices,p)
    end
end


function edge_indices_delaunay(simplices,neighbours)
    # Identify points which are edge points.
    edge_index = Array{Int64}(undef, 0)
    dim = size(simplices,2)
    for i = 1:size(simplices,1), j = 1:dim
            if neighbours[i,j]< 1
                push!(edge_index,simplices[i,setdiff(1:dim,j)]...)
            end
    end
    unique!(edge_index)
    return edge_index
end


# Convert vector of vectors to matrix if not matrix already
p_vals(P::Vector) = permutedims(hcat(P...))
p_vals(P::Matrix) = P


function neighbor_relabel(x::Int64,α_keep::Vector{Int64})
    if x == 0
        return x
    else
        return α_keep[x]
    end
end

function Delaunay_find(Positions,α)
    p = p_vals(Positions)
    tri = scipy_qhull.Delaunay(p)

    simplices = tri.simplices.+1 # convert to julia indexing
    neighbours = tri.neighbors.+1
    # simplices contains the triangles/tetrahedrons which form the Delaunay
    # neighbours tells you which of the triangles/tetrahedrons border each other
    # If they border the edge, they border 0


    # Compute the α value of each simplex
    α_val = alpha_shape_eval(simplices,p)
    # Determine α threshold value
    if α == 0
        #α = 2*median(α_val) #For Dom circumradius method
        #println("using default alpha of 2*(median alpha val) = ",α)
        α = median(α_val)+ std(α_val)
        println("using default filter of median(edge length)+stdev(edge length) = ",α)
    else
        println("using custom alpha = ",α)
    end

    if size(p,2) == 3
    # Figure out which simplices we want to keep, and also keep track of their
    # new indices in order to update the neighbour vector
        α_keep = zeros(Int64,size(α_val,1))
        counter = 1
        for i = 1:length(α_keep)
            indices = findall(y -> y < α, α_val[i,:]) #are all 3 edges < maxEdge (i.e. alpha)?
            #if α_val[i] < α
            if length(indices) == 3
                α_keep[i] = counter
                counter = counter + 1
            else
                α_keep[i] = 0
            end
        end
        # Relabel points in neighbours (some simplices are now edge points)
        f = x->neighbor_relabel(x,α_keep)
        neighbours = f.(neighbours)
        # Only keep the small simplices
        simplices = simplices[α_keep.>0,:]
        neighbours = neighbours[α_keep.>0,:]
    end
    # This function keeps track of which points are on the boundary
    edge_index = edge_indices_delaunay(simplices,neighbours)

    return p, simplices, neighbours, edge_index, α_val,α
end



function periodic_extend!(Coords;tol=0.6)
    # Applies the periodic boundary conditions so that no boundary effects apply
    # to any of the original N points. TODO: Currently assumes on a [0,1]^2 grid.
    if length(Coords[1]) == 2
        N = length(Coords)
        for i = 1:N
            for k1 in  -1:1, k2 in -1:1
                if k1^2 + k2^2 != 0
                    pos_mirrored = Coords[i] .+ [k1,k2]
                    if maximum(abs.(pos_mirrored .- 0.5)) < 0.5 + tol
                        push!(Coords,pos_mirrored)

                    end
                end
            end
        end
    elseif length(Coords[1]) == 3
        N = length(Coords)
        for i = 1:N
            for k1 in  -1:1, k2 in -1:1, k3 in -1:1
                if k1^2 + k2^2 + k3^2 != 0
                    pos_mirrored = Coords[i] .+ [k1,k2,k3]
                    if maximum(abs.(pos_mirrored .- 0.5)) < 0.5 + tol
                        push!(Coords,pos_mirrored)

                    end
                end
            end
        end
    end
end

function find_delaunay_network_core(Positions, periodic, α, tol, edge_keep)

    N = length(Positions)
    if periodic; periodic_extend!(Positions,tol=tol); end

    p, simplices, neighbrs, edge_index, α_val, α = Delaunay_find(Positions, α)

    if edge_keep
        not_edge = Vector(1:size(p,1)) # keep all indices
    else
        not_edge = setdiff(1:size(p,1), edge_index)
    end

    dim = size(p,2)
    return TopologicalNetwork(simplices, not_edge, dim, periodic, α,
               edge_keep, tol, N)

end
