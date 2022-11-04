using StatsBase:countmap
# Tools for investigating how the distribution of motif sizes vary. To be expanded

function tvec_dist(m::MotifArray)
    # find the distribution of tvecs lengths and probabilities
    len = length.(m.tvec)

    len_unique = sort(unique(len))
    freq_len = [length(findall(isequal(l),len)) for l in len_unique]
    p_len = freq_len/sum(freq_len)

    return len_unique,p_len
end

function tvec_dist(m::MotifDist)
    # find the distribution of tvecs lengths and probabilities
    len = length.(collect(keys(m.map)))
    count = collect(values(m.map))

    len_unique = sort(unique(len))
    freq_len = [sum(count[findall(isequal(l),len)]) for l in len_unique]
    p_len = freq_len/sum(freq_len)

    return len_unique,p_len
end

function resample(m::MotifDist,N::Int)
    # sample N points from the distribution m and return back a MotifArray
    m_arr = collect(keys(m.map))
    count = cumsum(collect(values(m.map)))
    probs = count/count[end]

    rnd_pts = rand(N)
    idx = [findfirst(x-> x >= l,probs) for l in rnd_pts]

    return MotifDist(countmap(m_arr[idx]),m.dim,m.r)
end

function moments_find(m,n=2)
    len_unique,p_len = tvec_dist(m)
    
    moments = zeros(n)
    for i = 1:n
        moments[i] = sum((len_unique.^i).* p_len)
    end
    return moments
end
