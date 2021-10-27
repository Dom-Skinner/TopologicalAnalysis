using Random

include("/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

Data_dir = "/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/"

Random.seed!(1234)
Positions = [randn(3) for k in 1:1000]
p, simplices, neighbours, edge_index = Delaunay_find(Positions)

for k = 1:1000
    k_nbhd = find_nbhd(simplices,k)
    tvec = topological_vec(k_nbhd,k)
end

using BenchmarkTools
using Profile
Profile.clear()

@btime topological_vec(k_nbhd,k)
Profile.print()
Juno.profiler()


k_nbhd = find_nbhd(simplices,5)
tvec = topological_vec(k_nbhd,5)
simplex = zeros(Int64,length(tvec),4)
t_vec_to_simplex!(tvec,simplex)

k_nbhd = [ 1 2 5 6; 1 2 4 6; 5 4 2 1; 5 4 1 3; 1 3 5 6; 1 6 4 3]

tvec = topological_vec(k_nbhd,1)
simplex = zeros(Int64,length(tvec),4)
t_vec_to_simplex!(tvec,simplex)

function perm_simplex!(nbhd)
    P = randperm(length(unique(nbhd)))
    for j=1:4,i=1:size(nbhd,1)
        nbhd[i,j] = P[nbhd[i,j]]
    end
end


using Random,LightGraphs
g = barabasi_albert(1000, 3, seed=123)
is_connected(g)
Random.seed!(123)
ρ1 = rand(nv(g))
ρ2 = rand(nv(g))
ρ1 = ρ1/sum(ρ1)
ρ2 = ρ2/sum(ρ2)

#= === To test iterative method vs Gaussian elimination
rem_self_edges!(g)
L = float.(laplacian_matrix(g))
D = float.(incidence_matrix(g,oriented=true))
W = ρ1 .- ρ2
J = D' * minres(L , W)
iter = sum(abs.(J))
J = D' * (L \ W)
gauss = sum(abs.(J))
=#


approx = distance_mat_lap(g,[ρ1, ρ2])[1,2]
exact = W_dist(g,ρ1 - ρ2)
