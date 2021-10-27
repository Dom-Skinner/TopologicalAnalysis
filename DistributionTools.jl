# Tools for investigating how the distribution of motif sizes vary. To be expanded

function tvec_dist(weight,dim="2D",N_ret=false)
    # find the distribution of tvecs lengths and probabilities

    # Takes too long to parse as integers, so do everything as strings
    ky = collect(keys(weight))
    len = [length(unique(split(k[2:end-1],", "))) for k in ky]

    p = float(collect(values(weight)))
    N = sum(p)
    len_unique = sort(unique(len))
    p_len = [sum(p[len .== l]) for l in len_unique]/N

    if N_ret
        return len_unique,p_len,N
    else
        return len_unique,p_len
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
        error("Todo")
    end
end

function find_dist_props(weights_arr)
    l_mean = zeros(length(weights_arr))
    l_var = zeros(length(weights_arr))
    l_skew = zeros(length(weights_arr))

    for i = 1:length(weights_arr)
        l,p = tvec_dist(weights_arr[i])
        μ,σ,γ = moments_find(l,p,3)
        l_mean[i] = μ
        l_var[i]  = σ
        l_skew[i] = γ
        println(100.0 *i / length(weights_arr))
    end
    return l_mean,l_var,l_skew
end
