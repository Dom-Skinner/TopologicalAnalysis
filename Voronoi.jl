#https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.html
#https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html
#https://graphics.stanford.edu/courses/cs268-11-spring/handouts/AlphaShapes/as_fisher.pdf

using PyCall
using LinearAlgebra: norm, det
using Statistics: median
scipy = pyimport("scipy.spatial.qhull")


function p_val(Positions)
   p = zeros(length(Positions),length(Positions[1]))
   for s = 1:length(Positions)
       p[s,:] = Positions[s]
   end
   return p
end

function edge_find_voronoi(p,index_points,regions)
   edge_index = []
   for s = 1:length(p[:,1])
        nbhd = regions[index_points[s]]
        for s2 in nbhd
            if s2 < 0
               push!(edge_index,s)
            end
         end
   end
   unique!(edge_index)
   return edge_index
end
#p = rand(1500,2)

function circumradius(e1,e2,e3)
    # see http://mathworld.wolfram.com/Circumradius.html
    a = sqrt(e1[1]^2+e1[2]^2)
    b = sqrt(e2[1]^2+e2[2]^2)
    c = sqrt(e3[1]^2+e3[2]^2)
    return a*b*c/sqrt( (a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))
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
    return sqrt(Dx^2+Dy^2+Dz^2-4*a*c)/(2*abs(a))
end

function alpha_shape2D!( α_val,simplices,p)
    for s = 1:length(simplices[:,1])
        e1 = p[simplices[s,1],:] - p[simplices[s,2],:]
        e2 = p[simplices[s,1],:] - p[simplices[s,3],:]
        e3 = p[simplices[s,2],:] - p[simplices[s,3],:]
        α_val[s] = circumradius(e1,e2,e3)
    end
end

function alpha_shape3D!( α_val,simplices,p)
    for s = 1:length(simplices[:,1])
        e1 = p[simplices[s,1],:]
        e2 = p[simplices[s,2],:]
        e3 = p[simplices[s,3],:]
        e4 = p[simplices[s,4],:]
        α_val[s] = circumradius3D(e1,e2,e3,e4)
    end
end

function alpha_shape_eval!(α_val,neighbours,simplices,p; α = 0)

    if size(p,2) == 3
        alpha_shape3D!(α_val,simplices,p)
    elseif size(p,2) == 2
         alpha_shape2D!(α_val,simplices,p)
    end
    if α == 0
        α = 2*median(α_val)         # Was previously α = 2. May need to be adjusted
    end
    #α = minimum([α;63.0]) # For the biofilms
    println("using alpha = ",α) # if there are density variations

    for s in 1:length(neighbours)
        if neighbours[s] > 0
            if  α_val[neighbours[s]] > α
                neighbours[s] = -1.
            end
        end
    end

    for i = 1:length(α_val)
        α_val[i] = α_val[i] < α
    end
end

function edge_indices_delaunay(simplices,neighbours)
    edge_index = Array{Int64}(undef, 0)
    for s = 1:length(simplices[:,1])
        for s2 = 1:length(simplices[1,:])
            if neighbours[s,s2]< 1
                 for s3 = 1:length(simplices[1,:])
                 # including opposite points
                 #push!(edge_index,simplices[s,s3])
                 # not including opposite points
                     if s2 != s3
                        push!(edge_index,simplices[s,s3])
                     end
                 end
            end
        end
    end
    unique!(edge_index)
    return edge_index
end


#=
function Voronoi_find(Positions)
   p = p_val(Positions)
   vor = scipy.Voronoi(p)
   index_points = vor.point_region.+1
   regions = vor.regions
   edge_index = edge_find_voronoi(p,index_points,regions)

   return p, index_points, regions, edge_index
end
=#

function Voronoi_find_shape(Positions)
   p = p_val(Positions)
   vor = scipy.Voronoi(p)
   index_points = vor.point_region.+1
   regions = vor.regions
   vert = vor.vertices
   ridge_vert = vor.ridge_vertices
   ridge_pts = vor.ridge_points.+1

   return index_points, regions, vert, ridge_vert, ridge_pts
end

function Delaunay_find(Positions; α = 0)
    p = p_val(Positions)
    tri = scipy.Delaunay(p)

    simplices = tri.simplices.+1 # convert to julia indexing
    neighbours = tri.neighbors.+1

    α_val = zeros(length(simplices[:,1]))

    alpha_shape_eval!(α_val,neighbours,simplices,p,α=α)

    edge_index = edge_indices_delaunay(simplices,neighbours)

    return p, simplices, neighbours, edge_index,α_val
end



function periodic_extend!(Coords,ref,tol=0.6)
    # Applies the periodic boundary conditions so that no boundary effects apply
    # to any of the original N points.
    if length(Coords[1]) == 2
        N = length(Coords)
        for i = 1:N
            for k1 in  -1:1, k2 in -1:1
                if k1^2 + k2^2 != 0
                    pos_mirrored = Coords[i] .+ [k1,k2]
                    if maximum(abs.(pos_mirrored .- 0.5)) < 0.5 + tol
                        push!(Coords,pos_mirrored)
                        push!(ref,i)

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
                        push!(ref,i)

                    end
                end
            end
        end
    end
end
