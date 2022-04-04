using JuMP, Gurobi
using SparseArrays
using LightGraphs
using Mosek, MosekTools
using Convex

function min_cost_flow(g,sources)
    # sources to be a vector of length nv(g), with zeros as necessary
    @assert LightGraphs.nv(g) == length(sources)
    # sum of sources to be zero to machine precision
    @assert abs(sum(sources)) < 1e-12

    es = collect(edges(g))
    Es = [src(e) for e in es]
    Ed = [dst(e) for e in es]
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
    es = collect(edges(g))
    Es = [src(e) for e in es]
    Ed = [dst(e) for e in es]
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
    #model = Model(optimizer_with_attributes(Mosek.Optimizer))
    #@variable(model,q[1:num_v,1:k+1])
    #@variable(model,F[1:num_e,1:k])
    q = Variable(num_v,k+1)
    F = Variable(num_e,k)

    #@objective(model,Min,sum([(F[e,i]^2/q[verts_inv[e_reduced[e,1]],i]) +
    #        (F[e,i]^2/q[verts_inv[e_reduced[e,2]],i+1])  for i = 1:k for e = 1:num_e]))

    obj_ = sum([quadoverlin(F[e,i], q[verts_inv[e_reduced[e,1]],i]) +
            quadoverlin(F[e,i], q[verts_inv[e_reduced[e,2]],i+1])  for i = 1:k for e = 1:num_e])

    con1 = q >= 0
    con2 = [-sum(F[findall(e_reduced[:,1] .== j),i]) + sum(F[findall(e_reduced[:,2] .== j),i]) == q[verts_inv[j],i+1]-q[verts_inv[j],i] for i = 1:k for j in nz_verts]
    con3 = q[:,1] == p0_v
    con4 = q[:,end] == p1_v
    con5 = F >= 0.0

    problem = minimize(obj_, [con1;con2;con3;con4;con5])
    solve!(problem, () -> Mosek.Optimizer())
    println(problem.status)

    q_full = zeros(length(p0),k+1)
    for i = 1:(k+1)
            q_full[nz_verts,i] .= q.value[:,i]
    end
    return q_full, k*problem.optval
end




function CFTDist_gurobi(g,p0,p1;k=10)

    # first solve min_cost_flow to identify non-zero edges to solve full minimization over
    v,J = min_cost_flow(g,p0 .- p1)

    # collect only the edges which were used for min_cost_flow
    es = collect(edges(g))
    Es = [src(e) for e in es]
    Ed = [dst(e) for e in es]
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
    #model = Model(optimizer_with_attributes(Mosek.Optimizer))
    model = Model(Gurobi.Optimizer)
    @variable(model,q[1:num_v,1:k+1])
    @variable(model,F[1:num_e,1:k])
    @variable(model,s1[1:num_e,1:k])
    @variable(model,s2[1:num_e,1:k])
    #q = Variable(num_v,k+1)
    #F = Variable(num_e,k)
    #s1 = Variable(num_e,k)
    #s2 = Variable(num_e,k)

    @objective(model,Min,sum(s1) + sum(s2))

    #obj_ = sum([quadoverlin(F[e,i], q[verts_inv[e_reduced[e,1]],i]) +
#            quadoverlin(F[e,i], q[verts_inv[e_reduced[e,2]],i+1])  for i = 1:k for e = 1:num_e])

    @constraint(model,q[:,1] .== p0_v)
    @constraint(model,q[:,end] .== p1_v)

    for i = 1:k , j in nz_verts
    @constraint(model,-sum(F[findall(e_reduced[:,1] .== j),i]) +
        sum(F[findall(e_reduced[:,2] .== j),i]) == q[verts_inv[j],i+1]-q[verts_inv[j],i] )
    end
    for i = 1:k, e = 1:num_e
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


N = 2000
g = barabasi_albert(N, 5)
p0 = rand(N)
p0 = p0/sum(p0)
p1 = rand(N)
p1 = p1/sum(p1)

CFTDist_gurobi(g,p0,p1)



model = Model(Gurobi.Optimizer)
@variable(model, x[1:2] >= 2)
@variable(model, 1 <= y[1:2] <= 3)
@variable(model, s[1:2] >= 0)

@constraint(model,  [s[i] == y[i] for i  = 1:2])

@objective(model, Min, sum(s))

optimize!(model)

value(x)
value(y)
value(s)
