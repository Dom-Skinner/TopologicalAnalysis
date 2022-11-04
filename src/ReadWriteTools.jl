using HDF5
using CSV, DataFrames
using StatsBase: countmap

function save(save_str,top::TopologicalNetwork)
    h5open(save_str, "w") do file
        write(file, "Type", "TopologicalNetwork")
        write(file, "simplices", top.simplices)
        write(file, "not_edge", top.not_edge)
        write(file, "dim", top.dim)
        write(file, "periodic", top.periodic)
        write(file, "alpha", top.alpha)
        write(file, "edge_keep", top.edge_keep)
        write(file, "tolerance", top.tolerance)
        write(file, "original_vertex_number", top.original_vertex_number)
        if ismissing(top.clusters)
            write(file, "clusters", "missing")
        else
            write(file, "clusters", top.clusters)
        end
    end
end

function save(save_str,motif::MotifArray)
    h5open(save_str, "w") do file
        write(file, "Type", "MotifArray")
        write(file, "idx", motif.idx)
        write(file, "tvec", motif_to_matrix(motif.tvec))
        write(file, "dim", motif.dim)
        write(file, "r", motif.r)
        if ismissing(motif.regions)
            write(file, "regions", "missing")
        else
            write(file, "regions", motif.regions)
        end
    end
end
function save(save_str,motif::MotifDist)
    h5open(save_str, "w") do file
        write(file, "Type", "MotifDist")
        write(file, "codes", motif_to_matrix(collect(keys(motif.map))))
        write(file, "values", collect(values(motif.map)))
        write(file, "dim", motif.dim)
        write(file, "r", motif.r)
    end
end

function save(save_str,fg::FlipGraph)
    E = collect(edges(fg.g))
    h5open(save_str, "w") do file
        write(file, "Type", "FlipGraph")
        write(file, "I", [e.src for e in E])
        write(file, "J", [e.dst for e in E])
        write(file, "codes", motif_to_matrix(collect(keys(fg.motif_code))))
        write(file, "values", collect(values(fg.motif_code)))
    end
end

function motif_to_matrix(tvec)

    if length(tvec) ==0
        return Int64[]
    end

    tvec_len = length.(tvec)
    full_arr = zeros(Int64,length(tvec),maximum(tvec_len))
    for i = 1:length(tvec)
        if tvec_len[i] > 0
            full_arr[i,1:tvec_len[i]] .= tvec[i]
        end
    end
    return full_arr
end

function load(save_str)
    input_type = h5read(save_str,"Type")
    if input_type == "TopologicalNetwork"
        return load_topological_network(save_str)
    elseif input_type == "MotifArray"
        return load_motif_array(save_str)
    elseif input_type == "MotifDist"
        return load_motif_dist_array(save_str)
    elseif input_type == "FlipGraph"
        load_flip_graph(save_str)
    end
end

function load_topological_network(save_str)

    simplices = h5read(save_str,"simplices")
    not_edge = h5read(save_str,"not_edge")
    dim = h5read(save_str,"dim")
    periodic = h5read(save_str,"periodic")
    alpha = h5read(save_str,"alpha")
    edge_keep = h5read(save_str,"edge_keep")
    tolerance = h5read(save_str,"tolerance")
    original_vertex_number = h5read(save_str,"original_vertex_number")
    clusters  = h5read(save_str,"clusters")
    if clusters == "missing"
        clusters = missing
    end
    return TopologicalNetwork(simplices, not_edge, dim, periodic, alpha,
        edge_keep, tolerance, original_vertex_number,clusters)
end

function matrix_to_motif(full_arr)
    tvec = [full_arr[i,:] for i in 1:size(full_arr,1)]
    return [r[r.>0] for r in tvec]
end

function load_motif_array(save_str)

    idx = h5read(save_str,"idx")
    tvec = matrix_to_motif(h5read(save_str,"tvec"))
    dim = h5read(save_str,"dim")
    r = h5read(save_str,"r")
    try
        regions  = h5read(save_str,"regions")
    catch
        regions = "missing"
    end
    if regions == "missing"
        regions = missing
    end

    return MotifArray(idx, tvec, dim, r,regions)
end


function load_motif_dist_array(save_str)

    codes = matrix_to_motif(h5read(save_str,"codes"))
    values = h5read(save_str,"values")
    map_ = Dict(codes .=> values)
    dim = h5read(save_str,"dim")
    r = h5read(save_str,"r")
    return MotifDist(map_,dim,r)
end

function load_flip_graph(save_str)

    I = h5read(save_str,"I")
    J = h5read(save_str,"J")
    # TODO: what is the vectorized way to do this? Not so important as this
    # is fast for unweighted graphs
    g = SimpleGraph(maximum(J))
    for i = 1:length(J)
        add_edge!(g,I[i],J[i])
    end

    codes = matrix_to_motif(h5read(save_str,"codes"))
    values = h5read(save_str,"values")
    code_to_idx = Dict(codes .=> values)
    return FlipGraph(g,code_to_idx)
end

