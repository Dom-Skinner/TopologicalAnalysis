using JuMP, Gurobi
using SparseArrays
using LightGraphs

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
