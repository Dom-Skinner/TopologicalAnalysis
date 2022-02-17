using Distributions
using ForwardDiff
using RandomMatrices
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

function moments_find(m::MotifArray,n=2)
    len = length.(m.tvec)

    μ = sum(len)/length(len)
    σ2 = sum((len .- μ).^2)/length(len)
    if n == 2
        return μ,σ2
    elseif n== 3
        skew = sum((len .- μ).^3)/length(len)
        return μ,σ2,skew
    else
        error("TODO: Implement moments higher than n=3")
    end
end
function moments_find(x,p,n=2)
    μ = sum(x.*p)
    σ2 = sum((x .- μ).^2 .*p)
    if n == 2
        return μ,σ2
    elseif n== 3
        skew = sum((x .- μ).^3 .*p)
        return μ,σ2,skew
    else
        error("TODO: Implement moments higher than n=3")
    end
end

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

TWcdf(x,beta) = cdf(TracyWidom(),x,beta=beta)
TWpdf(s::Number,beta=1) = ForwardDiff.derivative(x->TWcdf(x,beta), s)[1]
TWpdf(s::Vector,beta=1) = [ForwardDiff.derivative(x->TWcdf(x,beta), s1)[1] for s1 in s]
TWpdf(s::StepRangeLen,beta=1) = [ForwardDiff.derivative(x->TWcdf(x,beta), s1)[1] for s1 in s]
Gaussian(x,μ,σ) = 1/sqrt(2*π*σ^2)*exp.(-(x .- μ).^2 / 2 / σ^2)
log_normal(x,μ,σ) = 1 ./x /σ/sqrt(2*π).* exp.( -(log.(x) .- μ).^2 ./2 ./ σ^2)
