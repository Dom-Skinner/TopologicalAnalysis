using Distributed
using IterativeSolvers
using LinearAlgebra
using LightGraphsFlows
import LightGraphs
const lg = LightGraphs
import StatsBase
using SparseArrays
using Base.Threads
import Clp


function calculate_distance_matrix(fg::FlipGraph,file_out,
		motif_array; optimal_transport= true)
    # This function is a wrapper for all other functions in this file
    # from a flip graph and n motifs in  it returns the
    # n by n distance matrix

	fg = connected_flip_graph(fg)

    weight = [ret_weights(fg,m) for m in motif_array]

    d =  distance_mat(fg,weight,optimal_transport)
	CSV.write(file_out,DataFrame(d,str_arr))

end


function triangle_index(k)
	# Helps us to iterate through upper triangular part of matrix
	T = n ->Int(0.5*n*(n+1))
	kk = Int(ceil( sqrt(2*k + 0.25) - 0.5))
	return k - T(kk-1),kk + 1
end


function distance_mat(fg,weight,optimal_transport)
	n_needed = Int(ceil(0.5*length(weight)*(length(weight)-1)))
    d = zeros(length(weight),length(weight))
	W = [weight[triangle_index(i)[1]] .- weight[triangle_index(i)[2]]  for i = 1:n_needed]
	g = fg.g
	rem_self_edges!(g)


	if optimal_transport
		f = x -> W_dist(g,x)
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
    # Find the distribution of networks in the connected comp of the flip graph
	weight = avg_motif(motif).map

    weight_arr = zeros(Float64,nv(fg.g))
    key_arr = collect(keys(weight))
    for k in key_arr
        if haskey(fg.motif_code,k)
            weight_arr[fg.motif_code[k]] = weight[k]
        end
    end

    return weight_arr./sum(weight_arr)
end




function W_dist(g_undirected,S,ret_flow = false;tol=5e-5)

    # First create a directed version of g_undirected with edges in both
    # directions and capacities/costs equal to 1

    Nv = lg.nv(g_undirected)
    g = lg.DiGraph(Nv+2) # Create a flow-graph
    w = spzeros(Nv+2,Nv+2)
    capacity = spzeros(Nv+2,Nv+2)

    for e in lg.edges(g_undirected)
        lg.add_edge!(g,lg.src(e),lg.dst(e))
        lg.add_edge!(g,lg.dst(e),lg.src(e))
        capacity[lg.src(e),lg.dst(e)] = 1.
        capacity[lg.dst(e),lg.src(e)] = 1.
        w[lg.src(e),lg.dst(e)] = 1.
        w[lg.dst(e),lg.src(e)] = 1.
    end

    # Take the Nv+1th vertex as a source and Nv+2th as a sink.
    demand = spzeros(Nv+2,Nv+2)
    for i in 1:Nv
        if S[i] > 0
            lg.add_edge!(g,Nv+1,i)
            capacity[Nv+1,i] = S[i]*(1+tol)
        elseif S[i] < - 0.1*tol
            lg.add_edge!(g,i,Nv+2)
            demand[i,Nv+2] = -S[i]*(1-tol)
            capacity[i,Nv+2] = 1.

        end
    end
    println("Is connected?: ", lg.is_connected(g))
    println("Demand ", sum(demand))
    println("Sum S ", sum(S))
    println("Out capacity ", sum(capacity[Nv+1,:]))
    println("In capacity ", sum(capacity[:,Nv+2]))
    # call min cost flow
    flow = mincost_flow(g,spzeros(lg.nv(g)), capacity , w, Clp.Optimizer,
                edge_demand=demand, source_nodes=[Nv+1], sink_nodes=[Nv+2])
    if ret_flow
        return flow
    else
        return sum(abs.(flow[1:Nv,1:Nv]))
    end
end

function sample_dist(probs,N)
    w = StatsBase.Weights(probs)
    p_emp = zeros(length(probs))
    for i = 1:N
        p_emp[StatsBase.sample( w)] += 1
    end
    return p_emp/sum(p_emp)
end

function geo_core(w_in,flow,α,N)
    weight_geo = deepcopy(w_in)

    I,J,V = findnz(flow)
    for i = 1:length(I)
        if (I[i] > length(weight_geo)) || (J[i] > length(weight_geo))
            continue
        end
        weight_geo[I[i]] -= α*V[i]
        weight_geo[J[i]] += α*V[i]
    end
    return sample_dist(weight_geo,N)
end



function geodesic(w1,w2,α,network_save_file,N)
    w_vec_in = [w1;w2]
    fg = load(network_save_file)
	fg = connected_flip_graph(fg)
    weight = [ret_weights(fg, m) for m in w_vec_in]
    flow = W_dist(fg.g,weight[1] .- weight[2],true)
    w_geo = [geo_core(weight[1],flow,α0,N) for α0 in α]
    return w_geo
end

function find_reg_geo(p0,p1,g,k)
        flow = W_dist(g,p0 .- p1,true;tol=0.0)
        flow[findall(flow .< 0)] .= 0
        flow = flow[1:length(p0),1:length(p0)]
        I,J,V = findnz(flow)
        W1_cost = sum(abs.(flow))
        verts = unique(vcat(I,J))
        verts_inv = Dict(verts .=> 1:length(verts))
        num_v = length(verts)
        num_e = length(I)
        p0_v = p0[verts]
        p1_v = p1[verts]

        q = Variable(num_v,k+1)
        F = Variable(num_e,k)

        obj_ = sum([quadoverlin(F[e,i], q[verts_inv[I[e]],i]) +
                quadoverlin(F[e,i], q[verts_inv[J[e]],i+1])  for i = 1:k for e = 1:num_e])

        con1 = q >= 0
        con2 = [-sum(F[findall(I .== j),i]) + sum(F[findall(J .== j),i]) == q[verts_inv[j],i+1]-q[verts_inv[j],i] for i = 1:k for j in verts]
        con3 = q[:,1] == p0_v
        con4 = q[:,end] == p1_v
        con5 = F >= 0.0
        con6 = sum(F) <= W1_cost*1.08
        problem = minimize(obj_, [con1;con2;con3;con4;con5;con6])
        solve!(problem, () -> Mosek.Optimizer())
        println(problem.status)

        q_full = zeros(length(p0),k+1)
        for i = 1:(k+1)
                q_full[verts,i] .= q.value[:,i]
        end
        return q_full, k*problem.optval
end

function geodesic_reg(network_save_file,w1,w2,k;α=1.0,β=0.0)
        fg = load(network_save_file)
        w_vec_in = [w1;w2]
        weight = [ret_weights(fg,w_vec_in[i]) for i in 1:length(w_vec_in)]
        q,val = find_reg_geo(α*weight[1]+(1.0-α)*weight[2],β*weight[1]
                                +(1.0-β)*weight[2],g,k)
        return q,val
end

function distance_OT(g,W)

    try
    	d =  W_dist(g,W)
        return d
    catch
    	d = W_dist(g,-W)
        return d
    end
end

function rem_self_edges!(g)
    # Removes self loops
    for e in collect(lg.edges(g))
        if lg.src(e) == lg.dst(e)
            lg.rem_edge!(g,e)
        end
    end
end
