using Random

include("/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

Data_dir = "/Users/Dominic/Dropbox (MIT)/3DTopology/LocalCellularStructure/"

Random.seed!(1234)
Positions = [randn(3) for k in 1:1000]
p, simplices, neighbours, edge_index = Delaunay_find(Positions)

for k = 1:1000
    k_nbhd = find_nbhd(simplices,k)
    topological_vec(k_nbhd,k)
end

using BenchmarkTools
using Profile
Profile.clear()

@btime topological_vec(k_nbhd,k)
Profile.print()
Juno.profiler()

function find_nbhd(simplices,k)
    N = size(simplices,1)
    contains_k = falses(N)
    @inbounds for j = 1:4, i=1:N
        if k == simplices[i,j]
            contains_k[i] = true
        end
    end
    simplices[contains_k,:]
end

@btime find_nbhd(simplices,5)
