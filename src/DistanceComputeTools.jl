using Distributed
using IterativeSolvers
using LinearAlgebra
import LightGraphs
const lg = LightGraphs
import StatsBase
using SparseArrays
using Base.Threads

using JuMP, Gurobi


function calculate_distance_matrix(fg::FlipGraph, motif_array;
		optimal_transport= true)
    # This function is a wrapper for all other functions in this file
    # from a flip graph and n motifs in  it returns the
    # n by n distance matrix

	fg = connected_flip_graph(fg)

    weight = [ret_weights(fg,m) for m in motif_array]

    return distance_mat(fg,weight,optimal_transport)


end


function triangle_index(k)
	# Helps us to iterate through upper triangular part of matrix
	T = n ->Int(0.5*n*(n+1))
	kk = Int(ceil( sqrt(2*k + 0.25) - 0.5))
	return k - T(kk-1),kk + 1
end

function abs_value(fg,weight)
	# for caluclating avg. values for velocity and curvatures
	g = fg.g
	rem_self_edges!(g)
	L = float.(laplacian_matrix(g))
	D = float.(incidence_matrix(g,oriented=true))
	f = x -> sum(abs.(D' * minres(L , x)))
	d_flat = pmap(f,weight)
	return d_flat
end

function distance_mat(fg,weight,optimal_transport)
	n_needed = Int(ceil(0.5*length(weight)*(length(weight)-1)))
    d = zeros(length(weight),length(weight))
	W = [weight[triangle_index(i)[1]] .- weight[triangle_index(i)[2]]  for i = 1:n_needed]
	g = fg.g
	rem_self_edges!(g)


	if optimal_transport
		f = x -> min_cost_flow(g,x)
	else
		L = float.(laplacian_matrix(g))
		D = float.(incidence_matrix(g,oriented=true))
		f = x -> sum(abs.(D' * minres(L , x)))
	end

	d_flat = pmap(f,W)
	for i  = 1:n_needed
		d[triangle_index(i)[1],triangle_index(i)[2]] = d_flat[i]
		d[triangle_index(i)[2],triangle_index(i)[1]] = d_flat[i]
	end
	return d
end

function ret_weights(fg::FlipGraph,motif::MotifArray)
	return ret_weights(fg,avg_motif(motif))
end

function ret_weights(fg::FlipGraph,motif::MotifDist)
    # Find the distribution of networks in the connected comp of the flip graph
	weight = motif.map

    weight_arr = zeros(Float64,nv(fg.g))
    key_arr = collect(keys(weight))
    for k in key_arr
        if haskey(fg.motif_code,k)
            weight_arr[fg.motif_code[k]] = weight[k]
        end
    end

    return weight_arr./sum(weight_arr)
end


function sample_dist(probs,N)
	# should this be somewhere else? IDK?
    w = StatsBase.Weights(probs)
    p_emp = zeros(length(probs))
    for i = 1:N
        p_emp[StatsBase.sample( w)] += 1
    end
    return p_emp/sum(p_emp)
end

function rem_self_edges!(g)
    # Removes self loops
    for e in collect(lg.edges(g))
        if lg.src(e) == lg.dst(e)
            lg.rem_edge!(g,e)
        end
    end
end

function min_cost_flow(g,sources)
    # sources to be a vector of length nv(g), with zeros as necessary
    @assert LightGraphs.nv(g) == length(sources)
    # sum of sources to be zero to machine precision
    @assert abs(sum(sources)) < 1e-12

    es = collect(lg.edges(g))
    Es = [lg.src(e) for e in es]
    Ed = [lg.dst(e) for e in es]
    e = hcat(Es,Ed)
    e_rev = hcat(Ed, Es);
    e_tot = [e;e_rev];
    ne = size(e_tot,1);

    i = [e_tot[:,1];e_tot[:,2]];
    j = [(1:ne);(1:ne)];
    v = [ones(ne);-1*ones(ne)];

    model = Model(Gurobi.Optimizer)
    @variable(model, J[1:ne] >= 0)
    @objective(model, Min, sum(J))
    A = sparse(i,j,v)
    @constraint(model, A*J .== sources)

    optimize!(model)
    return objective_value(model), value.(J)
end

function CFTDist(g,p0,p1;k=10)

    # first solve min_cost_flow to identify non-zero edges to solve full minimization over
    v,J = min_cost_flow(g,p0 .- p1)

    # collect only the edges which were used for min_cost_flow
    es = collect(lg.edges(g))
    Es = [lg.src(e) for e in es]
    Ed = [lg.dst(e) for e in es]
    e = hcat(Es,Ed)
    e_rev = hcat(Ed, Es);
    e_tot = [e;e_rev];
    e_reduced = e_tot[J .> 0 ,:]
    nz_verts = unique(e_reduced)

    # collect only vertices which do not change mass between 0 and 1
    p0_v = p0[nz_verts]
    p1_v = p1[nz_verts]


    verts_inv = Dict(nz_verts .=> 1:length(nz_verts))
    num_v = length(nz_verts)
    num_e = size(e_reduced,1)

    # With graph specified, now set up the convex problem
    model = Model(Gurobi.Optimizer)
    @variable(model,q[1:num_v,1:k+1])
    @variable(model,F[1:num_e,1:k])
    @variable(model,s1[1:num_e,1:k])
    @variable(model,s2[1:num_e,1:k])

    # rewrite the objective as linear, with SOC constraints
    @objective(model,Min,sum(s1) + sum(s2))

    @constraint(model,q[:,1] .== p0_v)
    @constraint(model,q[:,end] .== p1_v)

    for i = 1:k , j in nz_verts
    @constraint(model,-sum(F[findall(e_reduced[:,1] .== j),i]) +
        sum(F[findall(e_reduced[:,2] .== j),i]) == q[verts_inv[j],i+1]-q[verts_inv[j],i] )
    end
    for i = 1:k, e = 1:num_e
        # [s;q;F] RotatedSecondOrderCone() is just constraining s*q >= F^2
        @constraint(model, [s1[e,i];q[verts_inv[e_reduced[e,1]],i]; F[e,i]] in RotatedSecondOrderCone() )
        @constraint(model, [s2[e,i];q[verts_inv[e_reduced[e,2]],i+1]; F[e,i]] in RotatedSecondOrderCone())
    end
    optimize!(model)


    q_full = zeros(length(p0),k+1)
    for i = 1:(k+1)
            q_full[nz_verts,i] .= value.(q[:,i])
    end
    return q_full, k*objective_value(model)
end
