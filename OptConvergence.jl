using MathProgBase, Gurobi
using Plots

# Finds the true solution of the minimum cost flow and the approx diffusion
# solution together with successively better approximations
g = barabasi_albert(1000, 3, seed=123)
is_connected(g)
Random.seed!(123)
ρ1 = rand(nv(g))
ρ2 = rand(nv(g))
Δρ = ρ1/sum(ρ1) .- ρ2/sum(ρ2)
L = float.(laplacian_matrix(g))
D = float.(incidence_matrix(g,oriented=true))
J = D' * (L \ Δρ) # The diffusion approx
approx = sum(abs.(J))
exact = W_dist(g,Δρ)



E = collect(edges(g))
shuffle!(E)
Cycles = []
span_tree = SimpleGraph(kruskal_mst(g)) # to find a cycle basis take spanning tree
N = 1:250:2500
sol_plot = []
for n in N
    for edge in E[1:n]
        if !has_edge(span_tree,edge)
            # cycle basis consists of edges not in spanning tree and cycle completed
            # by edges in planning tree. use a_star algorithm to find path
            push!(Cycles,[a_star(span_tree, dst(edge), src(edge)); edge])
        end
    end

    # Construct a sparse array containing what direction the edge is traversed in
    # the graph. This is incredibly slow and should never be used except for small
    # scripts.
    C = spzeros(ne(g),length(Cycles))
    for k in 1:length(Cycles)
        for j in 1:length(Cycles[k])
            for i = 1:length(E)
                if src(E[i]) == src(Cycles[k][j]) && dst(E[i]) == dst(Cycles[k][j])
                    C[i,k] = 1
                elseif dst(E[i]) == src(Cycles[k][j]) && src(E[i]) == dst(Cycles[k][j])
                        C[i,k] = -1
                end
            end
        end
    end

    # Set up the LP problem.
    A = [hcat(sparse(I,ne(g),ne(g)), C); hcat(sparse(I,ne(g),ne(g)), -C) ]
    b = vcat(J,-J)
    lb = vcat(zeros(ne(g)), repeat([-Inf],length(Cycles)))
    f = vcat(ones(ne(g)),zeros(length(Cycles)))
    # pass params as keyword arguments to GurobiSolver
    solution = linprog(f, A, '>', b, lb, Inf, GurobiSolver(Presolve=0))
    push!(sol_plot,solution.objval)
    println(n)
end

plot(N,sol_plot)
