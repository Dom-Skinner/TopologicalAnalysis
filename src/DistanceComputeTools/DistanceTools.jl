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
