using Distributions
using ForwardDiff
using RandomMatrices
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

#=
function find_dist_props(weights_arr)
    l_mean = zeros(length(weights_arr))
    l_var = zeros(length(weights_arr))
    l_skew = zeros(length(weights_arr))

    for i = 1:length(weights_arr)
        μ,σ,γ = moments_find(weights_arr[i],3)
        l_mean[i] = μ
        l_var[i]  = σ
        l_skew[i] = γ
    end
    return l_mean,l_var,l_skew
end
=#

#TODO: this should not be part of the package....
TWcdf(x,beta) = cdf(TracyWidom(),x,beta=beta)
TWpdf(s::Number,beta=1) = ForwardDiff.derivative(x->TWcdf(x,beta), s)[1]
TWpdf(s::Vector,beta=1) = [ForwardDiff.derivative(x->TWcdf(x,beta), s1)[1] for s1 in s]
TWpdf(s::StepRangeLen,beta=1) = [ForwardDiff.derivative(x->TWcdf(x,beta), s1)[1] for s1 in s]
Gaussian(x,μ,σ) = 1/sqrt(2*π*σ^2)*exp.(-(x .- μ).^2 / 2 / σ^2)
log_normal(x,μ,σ) = 1 ./x /σ/sqrt(2*π).* exp.( -(log.(x) .- μ).^2 ./2 ./ σ^2)
