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
    #@assert abs(sum(sources)) < 1e-12

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
	set_optimizer_attributes(model, "OptimalityTol" => 1e-9)
	set_optimizer_attributes(model, "FeasibilityTol" => 1e-9)
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
    return q_full, sqrt(k*objective_value(model))
end


function CFTD_perturbation_0_alt(g,p0,p1,r1)
	es = collect(lg.edges(g))
    Es = [lg.src(e) for e in es]
    Ed = [lg.dst(e) for e in es]
    e = hcat(Es,Ed)

	num_v = lg.nv(g)
    num_e = size(e,1)

	I = [e[:,1];e[:,2]]
	J = [1:size(e,1);1:size(e,1)]
	V = [-1*ones(size(e,1));ones(size(e,1))]
	Dt = sparse(I,J,V)
	D = sparse(J,I,V)

	Λ = [1/p0[e[i,1]] + 1/p0[e[i,2]] for i in 1:size(e,1)]

	L = Dt* diagm(1 ./Λ)*D

	if maximum(abs.(r1 .- p1))<1e-10
		λ = zeros(num_e)
	else
		λ = -lsmr(L, r1 .- p1)
	end

	J = - (1 ./Λ) .* (D*λ)

	return e, J, -0.5*sum(λ.*(r1 .- p1))
end

function CFTD_perturbation_0(g,p0,p1,r1)

    # collect edges
    es = collect(lg.edges(g))
    Es = [lg.src(e) for e in es]
    Ed = [lg.dst(e) for e in es]
    e = hcat(Es,Ed)
    e_rev = hcat(Ed, Es);
    e_tot = [e;e_rev];
	#e_tot = e_tot[(p0[e_tot[:,1]] .> 0) .& (p0[e_tot[:,2]] .> 0),:]
	p0[p0.==0] .= 1e-9
	num_v = lg.nv(g)
    num_e = size(e_tot,1)

	I = [e_tot[:,1];e_tot[:,2]]
	J = [1:size(e_tot,1);1:size(e_tot,1)]
	V = [-1*ones(size(e_tot,1));ones(size(e_tot,1))]
	Dt = sparse(I,J,V)
    # With graph specified, now set up the convex problem
    model = Model(Gurobi.Optimizer)
	#set_optimizer_attributes(model, "OptimalityTol" => 1e-7)
	#set_optimizer_attributes(model, "FeasibilityTol" => 1e-7)

    @variable(model,J1[1:num_e] >= 0)

    # rewrite the objective as linear, with SOC constraints
    @objective(model,Min,sum(J1[j]^2 *(1/p0[e_tot[j,1]] + 1/p0[e_tot[j,2]])/2
			for j = 1:num_e ))

	@constraint(model,Dt*J1 .== r1  .- p1)

    #for j in 1:num_v
	 #   @constraint(model,-sum(J1[findall(e_tot[:,1] .== j)]) +
	 #       sum(J1[findall(e_tot[:,2] .== j)]) == r1[j]  - p1[j])
    #end
	println("starting to optimize")
    optimize!(model)

	# After optimizing, only keep the edges that were actually used
	J1_used = zeros(lg.ne(g))
	edge_used = [maximum(value.(J1[i])) > maximum(value.(J1[i+lg.ne(g)])) for i = 1:lg.ne(g)]
	edge_used = vcat(edge_used, .!edge_used)
	J1_used = value.(J1[edge_used])

    return e_tot[edge_used,:],J1_used, objective_value(model)
end

function CFTD_perturbation_1(e,p0,J1,p1,r1,p2,r2)

    # collect edges
	num_v = length(unique(e))
	num_e = size(e,1)

	I = [e[:,1];e[:,2]]
	J = [1:size(e,1);1:size(e,1)]
	V = [-1*ones(size(e,1));ones(size(e,1))]
	Dt = sparse(I,J,V)

	if maximum(abs.(r2 .- p2))<1e-10
		J2 = zeros(num_e)
	else
		J2 = lsmr(Dt, r2 .- p2)
	end
	#F1 = Dt \ (r2 .- p2)

	chi = sum(  J2[i]*J1[i]*( 1/ p0[e[i,1]] + 1/p0[e[i,2]]) for i in 1:num_e)

	lambda = -0.25*sum( J1[i]^2 * ( (p1[e[i,1]] +r1[e[i,1]])/p0[e[i,1]]^2 +
		(p1[e[i,2]] +r1[e[i,2]])/p0[e[i,2]]^2) for i in 1:num_e)
	I1 = chi+lambda
    return I1
end


function CFTD_perturbation_2(e,p0,J1,p1,r1,p2,r2)

	@assert size(e,1) == length(J1)

	num_e = size(e,1)
	verts = unique(e)
	num_v = length(verts)
	verts_inv = Dict(verts .=> 1:num_v)
	e_inv_u = [verts_inv[e[j,1]] for j = 1:num_v]
	e_inv_v = [verts_inv[e[j,2]] for j = 1:num_v]

	model = Model(Gurobi.Optimizer)

	@variable(model,s[1:num_v])
	@variable(model,F[1:num_e])
	@variable(model,G[1:num_e])

	p0_u = p0[e[:,1]]
	p0_v = p0[e[:,2]]

	p1_u = p1[e[:,1]]
	p1_v = p1[e[:,2]]
	r1_u = r1[e[:,1]]
	r1_v = r1[e[:,2]]

	p2_u = p2[e[:,1]]
	p2_v = p2[e[:,2]]
	r2_u = r2[e[:,1]]
	r2_v = r2[e[:,2]]


    for j in verts
    @constraint(model,-sum(F[findall(e[:,1] .== j)]) +
        sum(F[findall(e[:,2] .== j)]) == r2[j] -p2[j] + s[verts_inv[j]])
	@constraint(model,-sum(G[findall(e[:,1] .== j)]) +
        sum(G[findall(e[:,2] .== j)]) == -2*s[verts_inv[j]])
    end

	χ₁ = J1 .* ( 0.5*(p1_u .+ r1_u)./p0_u.^2 .+ 0.5*(p1_v .+ r1_v)./p0_v.^2 )
	χ₂ = J1 .* ( (p1_u/6 .+ r1_u/3)./p0_u.^2 .+ (p1_v/6 .+ r1_v/3)./p0_v.^2 )
	χ₃ = 0.5*(p2_u .+ r2_u)./p0_u.^2 .+ 0.5*(p2_v .+ r2_v)./p0_v.^2

	λ₂ = (1/12)*sum(J1.^2  .*( (p1_u.^2 .+  p1_u.*r1_u .+ r1_u.^2) ./ p0_u.^3 .+
	  		(p1_v.^2 .+  p1_v.*r1_v .+ r1_v.^2) ./p0_v.^3))


	@objective(model,Min,-sum(F.*χ₁ .+ G.*χ₂) +
		0.5*sum( (F.^2 .+ F.*G .+ (1/3)*G.^2) .*( 1 ./p0_u + 1 ./p0_v)) -
		0.5*sum(J1 .^2 .*(χ₃ + (1/6)*s[e_inv_u] ./ p0_u.^2 + (1/6)*s[e_inv_v] ./ p0_v.^2)  ))
	println(model)
	optimize!(model)

    return objective_value(model) +λ₂
end



function CFTD_perturbation_2_alt(e,p0,J1,p1,r1,p2,r2)

	@assert size(e,1) == length(J1)

	num_e = size(e,1)
	verts = unique(e)
	num_v = length(verts)

	p0_u = p0[e[:,1]]
	p0_v = p0[e[:,2]]

	p1_u = p1[e[:,1]]
	p1_v = p1[e[:,2]]
	r1_u = r1[e[:,1]]
	r1_v = r1[e[:,2]]

	p2_u = p2[e[:,1]]
	p2_v = p2[e[:,2]]
	r2_u = r2[e[:,1]]
	r2_v = r2[e[:,2]]

	I = [e[:,1];e[:,2]]
	J = [1:size(e,1);1:size(e,1)]
	V = [-1*ones(size(e,1));ones(size(e,1))]
	Dt = sparse(I,J,V)
	D = sparse(J,I,V)

	Λ = 1 ./p0_u .+ 1 ./p0_v
	χ₁ = J1 .* ( 0.5*(p1_u .+ r1_u)./p0_u.^2 .+ 0.5*(p1_v .+ r1_v)./p0_v.^2 )
	χ₂ = J1 .* ( (p1_u/6 .+ r1_u/3)./p0_u.^2 .+ (p1_v/6 .+ r1_v/3)./p0_v.^2 )
	χ₃ = 0.5*(p2_u .+ r2_u)./p0_u.^2 .+ 0.5*(p2_v .+ r2_v)./p0_v.^2

	λ₂ = (1/12)*sum(J1.^2  .*( (p1_u.^2 .+  p1_u.*r1_u .+ r1_u.^2) ./ p0_u.^3 .+
	  		(p1_v.^2 .+  p1_v.*r1_v .+ r1_v.^2) ./p0_v.^3)) - 0.5*sum(J1.^2 .*χ₃)

	Γ = zeros(num_v)
	for i = 1:num_e
		Γ[e[i,1]] = Γ[e[i,1]] - J1[i]^2/12/p0[e[i,1]]^2
		Γ[e[i,2]] = Γ[e[i,2]] - J1[i]^2/12/p0[e[i,2]]^2
	end

	m = Dt*(χ₁./Λ) .- (r2 .- p2)
	L = Dt* diagm(1 ./Λ)*D

	if maximum(abs.(m))<1e-10
		μ₁ = zeros(num_e)
	else
		μ₁ = lsmr(L, m)
	end



	μ₂ = 0.5*(μ₁ .- Γ)
	G = 6*(1 ./Λ).*(2*χ₂ .-χ₁ .+ D*μ₁ .- 2*D*μ₂)
	F = (1 ./Λ).*(χ₁ .- D*μ₁) .- 0.5*G
	s = -0.5*Dt*G
	I2 = -sum(F.*χ₁) -sum(G.*χ₂) + 0.5*sum(Λ .* F.^2) + 0.5*sum(Λ .* F .*G) +(1/6)*sum(Λ .* G.^2) +λ₂ + sum(Γ.*s)
    return I2
end


function CFTD_curvature(g,p,dp,d2p)

	@assert size(p) == size(dp)
	@assert size(p) == size(d2p)

    z = zeros(size(p))

    e,J1f,I0f = CFTD_perturbation_0_alt(g,p,z,dp)
	#error()
    I1f = CFTD_perturbation_1(e,p,J1f,z,dp,z,d2p/2)
    I2f = CFTD_perturbation_2_alt(e,p,J1f,z,dp,z,d2p/2)

    I2c = CFTD_perturbation_2_alt(e,p,2*J1f,-dp,dp,d2p/2,d2p/2)

    return 4*(I0f)^(-0.75) * sqrt( (I2f-0.25*I2c)/sqrt(I0f) - 0.25*I1f^2/(I0f)^1.5)
end
