function JS_div(dict1,dict2)
    N1 = sum(values(dict1)); N2 = sum(values(dict2))
    sum_tot = 0.
    for k1 in keys(dict1)
        y = dict1[k1]/N1
        z = 0.5*y +0.5*get(dict2,k1,0)/N2
        sum_tot += y*log(y/z)
    end
    for k2 in keys(dict2)
        x = dict2[k2]/N2
        z = 0.5*x + 0.5*get(dict1,k2,0)/N1
        sum_tot += x*log(x/z)
    end
    return sqrt(sum_tot/(2*log(2)))
end

function weight_in(Data_dir,str_arr,m=-1)
    if m > 0
        return vcat([readin(Data_dir*s,m) for s in str_arr]...)
    else
        return vcat([readin(Data_dir*s) for s in str_arr]...)
    end
end

function fill_JS_distance_mat(weight)
    d = zeros(length(weight),length(weight))
    for i = 1:length(weight), j= (i+1):length(weight)
        d[i,j] = JS_div(weight[i],weight[j])
        d[j,i] = d[i,j]
    end
    return d
end


function fill_W_distance_mat(str)
    d_in  = read(str, String)
    d_in = Meta.parse(d_in).args
    d = zeros(length(d_in),length(d_in))
    for i = 1:length(d_in)
        d[i,:] = d_in[i].args
    end
    return d
end

function subsample_dist(w_in,M)
    key_in = collect(keys(w_in))
    vals = collect(values(w_in))
    c = cumsum(vals)
    D = Dict(i=>findfirst(x->x≥i,c) for i=1:sum(vals))
    I = rand(1:sum(vals),M)
    w_subsampled = Dict()
    for i in I
        j  = D[i]
        w_subsampled[key_in[j]] = 1 + get(w_subsampled,key_in[j],0)
    end
    return w_subsampled
end

function fill_SN_distance_mat(weight)
    # This function calculates a simple L2 distance based on the number of sides
    s_tot = []
    for i = 1:length(weight)
    s_vec = zeros(20)
    w_in = [Int.(Meta.parse(w).args) for w in collect(keys(weight[i]))]
    freq = collect(values(weight[i]))
    for j in 1:length(w_in)
        s_vec[sum(w_in[j] .== 1) - 1] += freq[j]
    end
    push!(s_tot,s_vec./sum(s_vec))
    end

    D = zeros(length(weight),length(weight))
    for i = 1:length(weight),j=(i+1):length(weight)
        D[i,j] = sqrt.(sum(abs2,s_tot[i].-s_tot[j]))
        D[j,i] = D[i,j]
    end
    return D
end
