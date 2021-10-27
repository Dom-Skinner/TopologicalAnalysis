using CSV, DataFrames
using StatsBase: countmap

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
function readin(Data_dir_str,N)
    if N > 0
        W_arr = Array{Dict}(undef,0)
        for i = 1:N
            dat_in = CSV.read(Data_dir_str*string(i)*"_avg.txt")
            push!(W_arr,Dict(dat_in.codes .=> dat_in.freq))
        end
        return W_arr
    else
        dat_in = CSV.read(Data_dir_str)
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

function amalg2(w_tot)
    count_tot = Dict{String,Int64}()
    for i = 1:length(w_tot)
        for k in collect(keys(w_tot[i]))
            count_tot[k]= get(count_tot, k,0) + Int(w_tot[i][k])
        end
    end
    return sort(collect(count_tot), by = tuple -> last(tuple), rev=true)
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

function write_total(Positions, Weinberg,S,Data_dir_str)
    df = DataFrame(x = [Positions[k][1] for k in 1:length(Positions)],
    y = [Positions[k][2] for k in 1:length(Positions)], weinberg = Weinberg, S = S)
    CSV.write(Data_dir_str*".txt",  df)
end
function write_total(Positions, tvec,Data_dir_str)
    df = DataFrame(x = [Positions[k][1] for k in 1:length(Positions)],
                   y = [Positions[k][2] for k in 1:length(Positions)],
                   z = [Positions[k][3] for k in 1:length(Positions)],
                   tvec = tvec)
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
