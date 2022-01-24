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
    end
end
function save(save_str,motif::MotifDist)
    h5open(save_str, "w") do file
        write(file, "Type", "MotifDist")
        write(file, "codes", motif_to_matrix(collect(keys(motif.map))))
        write(file, "values", collect(values(motif.map)))
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
    tvec_len = length.(tvec)
    full_arr = zeros(Int64,length(tvec),maximum(tvec_len))
    for i = 1:length(tvec)
        full_arr[i,1:tvec_len[i]] .= tvec[i]
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

    return MotifArray(idx, tvec, dim, r)
end


function load_motif_dist_array(save_str)

    codes = matrix_to_motif(h5read(save_str,"codes"))
    values = h5read(save_str,"values")
    map_ = Dict(codes .=> values)
    return MotifDist(map_)
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






#=

function weight_in(Data_dir,str_arr,m=-1)
    # a wrapper to readin that takes a string array
    if m > 0
        return vcat([readin(Data_dir*s,m) for s in str_arr]...)
    else
        return vcat([readin(Data_dir*s) for s in str_arr]...)
    end
end

function get_weights_in_dir(data_directory)
    # reads a directory and gets the weights of all avg.txt files in that directory
    if isa(data_directory, Array)
        return vcat([get_weights_in_dir(dir) for dir in data_directory]...)
    end

    text_files = filter(x->occursin("avg.txt",x), readdir(data_directory)) # look for avg.txt files
    if length(text_files) == 0
        return
    end
    weight_array = Dict[]
    for file in text_files
        raw_data = CSV.read(data_directory*file)
        push!(weight_array,Dict(raw_data.codes .=> raw_data.freq))
    end
    return weight_array
end

function readin!(w_tot,freq_tot,Data_dir_str,N)
    for n in 1:N
      dat_in = CSV.read(Data_dir_str*"_avg_"*string(n)*".txt")
      push!(w_tot, dat_in.weinberg_found)
      push!(freq_tot,dat_in.freq)
    end
end

function readin(Data_dir_str)
    dat_in = CSV.read(Data_dir_str,DataFrame)
    return Dict(dat_in.codes .=> dat_in.freq)
end

#=
function readin(Data_dir_str,N)
    if N > 0
        W_arr = Array{Dict}(undef,0)
        for i = 1:N
            dat_in = CSV.read(Data_dir_str*string(i)*"_avg.txt")
            push!(W_arr,Dict(dat_in.codes .=> dat_in.freq))
        end
        return W_arr
    else
        dat_in = CSV.read(Data_dir_str,DataFrame)
        return Dict(dat_in.codes .=> dat_in.freq)
    end
end

function readin(Data_dir_str)
    Data_dir =  SubString(Data_dir_str,1,findlast("/", Data_dir_str)[1] )
    str =  SubString(Data_dir_str,(findlast("/", Data_dir_str)[1] +1),length(Data_dir_str))
    DIR = readdir(Data_dir)
    header = [occursin(str,d)  for d in DIR]
    tail = [SubString(d,(length(d)-6):length(d)) for d in DIR]
    avg_files = DIR[header .& (tail.=="avg.txt")]
    Nums = sort!([parse(Int64,chop(s,head=length(str),tail=8)) for s in avg_files])
    W_arr = Array{Dict}(undef,0)
    for i in Nums
        dat_in = CSV.read(Data_dir_str*string(i)*"_avg.txt")
        push!(W_arr,Dict(dat_in.codes .=> dat_in.freq))
    end
    return W_arr
end
=#

function amalg2(w_tot)
    count_tot = Dict{String,Int64}()
    for i = 1:length(w_tot)
        for k in collect(keys(w_tot[i]))
            count_tot[k]= get(count_tot, k,0) + Int(w_tot[i][k])
        end
    end
    return sort(collect(count_tot), by = tuple -> last(tuple), rev=true)
end


function combine_distributions(dict_array)
    count_tot = Dict{String,Int64}()
    for i = 1:length(dict_array)
        for k in collect(keys(dict_array[i]))
            count_tot[k]= get(count_tot, k,0) + Int(dict_array[i][k])
        end
    end
    return Dict(sort(collect(count_tot), by = tuple -> last(tuple), rev=true))
end


function amalg2(w_tot,freq_tot)
    w_col = []
    w_freq = []
    for i = 1:length(w_tot)
            append!(w_col,w_tot[i])
            append!(w_freq,freq_tot[i])
    end
    return accumarray_s(w_col,w_freq)
end

function write_total(idx, tvec,Data_dir_str)
    df = DataFrame(index = idx, topological_vec = tvec)
    CSV.write(Data_dir_str*".txt",  df)
end

function write_avg(Data_dir_str)
    dat = CSV.read(Data_dir_str*".txt")
    code_tot = dat.weinberg
    code_freq = countmap(code_tot)
    df = DataFrame(freq = collect(values(code_freq)), codes = collect(keys(code_freq)))
    CSV.write(Data_dir_str*"_avg.txt",  df)
end
function write_avg(Weinberg,Data_dir_str)
    df = DataFrame(codes = collect(keys(Weinberg)), freq = collect(values(Weinberg)))
    CSV.write(Data_dir_str*"_avg.txt",  df)
end

function get_files_dir(Data_dir; c = 10)
      # Get all relevant file names in a directory
      DIR = readdir(Data_dir)
      str_arr = DIR[[occursin("_avg.txt",d) for d in DIR]]
      str_arr = [String(SubString(d,1:(length(d)-c))) for d in str_arr].*"_"
      return  str_arr
end
=#
