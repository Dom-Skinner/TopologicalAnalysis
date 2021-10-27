using MPI
include("/home/gridsan/dskinner/Packing/LocalCellularStructure/LocalCellularStructure.jl")
using .LocalCellularStructure

include("Parallel_distance_helper.jl")
using LightGraphs
using SparseArrays
using LinearAlgebra
using IterativeSolvers


MPI.Init()
comm = MPI.COMM_WORLD
my_rank = MPI.Comm_rank(comm)
n_procs = MPI.Comm_size(comm)
root = 0

g,weight = initialize()
n_needed = Int(ceil(0.5*length(weight)*(length(weight)-1)))
my_indices = my_rank+1:n_procs:n_needed

my_d = zeros(length(my_indices))

rem_self_edges!(g)
L = float.(laplacian_matrix(g))
D = float.(incidence_matrix(g,oriented=true))

for k in 1:length(my_indices)
    i,j =  triangle_index(my_indices[k])
    W = weight[i] .- weight[j]
    J = D' * minres(L , W) 
    my_d[k] = sum(abs.(J))
    println("Did calculation $(my_indices[k])")
end

# Gather all calculations onto process 1
if my_rank > 0
    # send to 0
    msg_d = MPI.serialize(my_d)
    msg_i = MPI.serialize(my_indices)
    sreq = MPI.Send(length(msg_d), 0, my_rank+n_procs, comm) # First send the length of the message
    sreq = MPI.Send(length(msg_i), 0, my_rank+(n_procs*2), comm) # First send the length of the message
    sreq = MPI.Send(msg_d, 0, my_rank+(n_procs*3), comm) # Then send the message
    sreq = MPI.Send(msg_i, 0, my_rank+(n_procs*4), comm) # Then send the message
else
    all_d = Array{Array{UInt8,1},1}(undef,n_procs-1)
    all_i = Array{Array{UInt8,1},1}(undef,n_procs-1)
    for i = 1:n_procs-1
        # receive
        global all_d
        msglen_d = Array{Int64,1}(undef,1)
        msglen_i = Array{Int64,1}(undef,1)
        rreq = MPI.Recv!(msglen_d, i,  i+n_procs, comm) # First receive the length of the message
        rreq = MPI.Recv!(msglen_i, i,  i+(n_procs*2), comm) # First receive the length of the message
        all_d[i] = Array{UInt8}(undef,msglen_d[1])
        all_i[i] = Array{UInt8}(undef,msglen_i[1])
        rreq = MPI.Recv!(all_d[i], i,  i+(n_procs*3), comm) # Teh receive the message
        rreq = MPI.Recv!(all_i[i], i,  i+(n_procs*4), comm) # Teh receive the message
    end
    all_d = vcat(MPI.deserialize.(all_d)...,my_d...)
    all_i = vcat(MPI.deserialize.(all_i)...,my_indices...)
    d0 = zeros(length(weight),length(weight))
    for k in 1:length(all_d)
	i,j = triangle_index(all_i[k])
	d0[i,j] = all_d[k]
	d0[j,i] = d0[i,j]
    end
    println("Task completed")
    println(d0)
end
