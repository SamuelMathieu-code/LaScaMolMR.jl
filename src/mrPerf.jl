using Distributions
using GLM
using LinearAlgebra
using Random
using Statistics
using LoopVectorization


    ###########################
    #        Utilities        #
    ###########################


# Raw weighted median effect size estimate
function _wm_estimate(β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                     se_β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                     β_X::AbstractVector{<: Union{AbstractFloat, Missing}}
                     )::Float64

θ_v = β_Y ./ β_X
order = sortperm(θ_v)
θ_v, se_βy_ordered, βx_ordered = θ_v[order], se_β_Y[order], β_X[order]
w_v = abs2.(βx_ordered ./ se_βy_ordered)
p_v = accumulate(+, w_v) - (w_v / 2)
p_v /= sum(w_v)
p = findfirst(x -> x ≥ 0.5, p_v)
θ_est = (p != 1) ? θ_v[p-1] + (θ_v[p] - θ_v[p-1])*(0.5 - p_v[p-1])/(p_v[p] - p_v[p-1]) : 0
return θ_est
end



# Numerical esimation of weighted median effect size standard error
function _bootstrap_se_wm(β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
    se_β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
    β_X::AbstractVector{<: Union{AbstractFloat, Missing}},
    se_β_X::AbstractVector{<: Union{AbstractFloat, Missing}}, 
    iterations = 5000,
    seed = 42)::Float64

    Random.seed!(seed)

    dx = MvNormal(Vector{Float32}(β_X), LinearAlgebra.Diagonal(map(abs2, Vector{Float32}(se_β_X))))
    dy = MvNormal(Vector{Float32}(β_Y), LinearAlgebra.Diagonal(map(abs2, Vector{Float32}(se_β_Y))))
    θ_est_v = Vector{Float32}(undef, iterations)
    for i in 1:iterations
        θ_est_v[i] = _wm_estimate(rand(dy), se_β_Y, rand(dx))
    end
    return Statistics.std(θ_est_v)
end


    #########################################
    #                 MrPerf                #
    #########################################

# eventuality : make lm more efficient by not calculating sumstats (only estimate)


"""
Struct encapsulating the outputs of a Mendelian Randomization analysis

fields :

`nivs` : number of ivs included in regression \\
`effect` : effect size estimate ̂γ₀ \\
`se_effect` : standard error of effect size estimate \\
`ci_low` : lower bound of confidence interval for effect size \\
`ci_high` : higher bound of confidence interval for effect size \\
`p` : effect size p-value for H₀ : γ₀ = 0 \\
`intercept` : intercept estimate ̂γ₁ \\
`p_intercept` : p-value for H₀ : γ₁ = 0 \\
`ci_low_intercept` : lower bound of confidence interval for intercept \\
`ci_high_intercept` : higher bound of confidence interval for intercept \\
`heter_stat` : heterogenetity statistic corresponding to Cochran ̂t \\
`heter_p` : heterogeneity p-value for H₀ : t = 0
"""
struct mr_output
    nivs::Int
    effect::Float64
    se_effect::Float64
    ci_low::Float64
    ci_high::Float64
    p::Float64
    intercept::Float64
    p_intercept::Float64
    ci_low_intercept::Float64
    ci_high_intercept::Float64
    heter_stat::Float64
    heter_p::Float64
end


"""
Default values constructor for mr_output with no heterogeneity stats
"""
function mr_output(n::Int, 
    effect::Float64 = NaN,
    se_effect::Float64 = NaN,
    ci_low::Float64 = NaN,
    ci_high::Float64 = NaN,
    p::Float64 = NaN,
    intercept::Float64 = NaN,
    p_intercept::Float64 = NaN,
    ci_low_intercept::Float64 = NaN,
    ci_high_intercept::Float64 = NaN)::mr_output

    mr_output(n, effect, se_effect, ci_low, ci_high, p, intercept, p_intercept, ci_low_intercept, ci_high_intercept, NaN, NaN)
    
end


"""
Wald ratio for Mendelian Randomization with a single instrumental variable
     y is outcome, x is exposure
"""
function mr_wald(β_Y::AbstractFloat, 
                 se_β_Y::AbstractFloat, 
                 β_X::AbstractFloat, 
                 α::AbstractFloat = 0.05)::mr_output # À véerifier que la distribustion normale convient!!!!!!!!
    θ = β_Y / β_X
    se_θ = se_β_Y / abs(β_X)
    dh = Normal(0, se_θ)
    dobs = Normal(θ, se_θ)
    p = 2*cdf(dh, -abs(θ))
    ci_low, ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    
    return mr_output(1, θ, se_θ, ci_low, ci_high, p)
end

function mr_wald(β_y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                 se_β_y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                 β_x::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                 se_β_X::AbstractVector{<: Union{AbstractFloat, Missing}} = Vector{Float32}([]),
                 α::AbstractFloat = 0.05)::mr_output
    return mr_wald(β_y[1], se_β_y[1], β_x[1], α)
end


"""
Inverse variance weighted linear regression with simple weights (se(B_Y)^-2) Mendelian Randomization
    Y is outcome, X is exposure

currently waiting for GLM.jl PR#487 to be merged to use analytical weights instead of doing calculations twice...
"""
function mr_ivw(β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                se_β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                β_X::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                se_β_X::AbstractVector{<: Union{AbstractFloat, Missing}} = Vector{Float32}([]),
                α::AbstractFloat = 0.05)::mr_output # where F <: Union{AbstractFloat, Missing}

    
    m = length(β_X)
    if !(length(β_Y) == length(se_β_Y) == m) throw(ArgumentError("β_x, β_y, se_β_X, se_β_Y Vectors must be of same length")) end
    if m < 2
        return mr_output(m)
    end
    # regression
    w = inv.(abs2.(se_β_Y))
    regressor = lm(@formula(β_Y ~ 0 + β_X), (;β_X, β_Y), wts = w)
    θivw_est = coef(regressor)[1]

    if any(isnan.(stderror(regressor)))
        return mr_output(m)
    end
    
    ϵ = residuals(regressor)
    u = sqrt.(w).*ϵ
    X = sqrt.(w).*β_X
    se_θivw_est = sqrt((u'*u)*(X'*X)^(-1)/(m-1))
    σ = sqrt(sum(abs2.(u))/(m-1))
    se_θivw_est /= min(σ, 1)

    dh = Normal(0, se_θivw_est)
    dobs = Normal(θivw_est, se_θivw_est)
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θivw_est))

    heter_stat = sum((ϵ ./ se_β_Y).^2)
    chisq = Chisq(m - 1)
    heter_p = 1 - cdf(chisq, heter_stat)

    return mr_output(length(β_Y), θivw_est, se_θivw_est, θ_ci_low, θ_ci_high, p, NaN, NaN, NaN, NaN, heter_stat, heter_p) 

end


"""
Egger Mendelian Randomization
    Y is outcome, X is exposure

currently waiting for GLM.jl PR#487 to be merged to use analytical weights instead of doing calculations twice...
"""
function mr_egger(β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                  se_β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                  β_X::AbstractVector{<: Union{AbstractFloat, Missing}}, 
                  se_β_X::AbstractVector{<: Union{AbstractFloat, Missing}} = Vector{Float32}([]),
                  α::AbstractFloat = 0.05)::mr_output
    
    # regression
    m = length(β_X)
    if !(length(β_Y) == length(se_β_Y) == m) throw(ArgumentError("β_x, β_y, se_β_X, se_β_Y Vectors must be of same length")) end
    if m < 3
        return mr_output(m)
    end
    w =  inv.(abs2.(se_β_Y))
    β_Y_abs = sign.(β_X).*β_Y
    β_X_abs = abs.(β_X)
    regressor = lm(@formula(β_Y_abs ~ β_X_abs), (;β_X_abs, β_Y_abs), wts = w)
    θ_est = coef(regressor)
    
    if any(isnan.(stderror(regressor)))
        return mr_output(m)
    end


    ϵ = residuals(regressor)
    X = sqrt.(w) .* [ones(m) β_X_abs]
    u = sqrt.(w) .* ϵ
    se = (u'*u)*(X'*X)^(-1)/(m-2)
    σ = sqrt(sum(abs2.(u))/(m-2))
    se_θ_est = sqrt.(diag(se)) ./ min(σ, 1)

    # index 1 is intercept and index two is effect (se_θ_est, θ_est)
    dh = Normal(0, se_θ_est[2])
    dobs = Normal(θ_est[2], se_θ_est[2])
    dintercept = Normal(0, se_θ_est[1])
    
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θ_est[2]))

    p_intercept = 1 - 2*(cdf(dintercept, -abs(θ_est[1])))
    ci_low_intercept, ci_high_intercept = quantile(dintercept, α/2), quantile(dintercept, 1-α/2)

    # heterogeneity
    heter_stat = sum(abs2.(ϵ ./ se_β_Y))
    chisq = Chisq(m - 2)
    heter_p = 1 - cdf(chisq, heter_stat)

    return mr_output(length(β_Y), θ_est[2], se_θ_est[2], θ_ci_low, θ_ci_high, p, 
                     θ_est[1], p_intercept, ci_low_intercept, ci_high_intercept, 
                     heter_stat, heter_p)
end


"""
Weighted Median Mendelian Randomization (Bowden et al., 2015)
    Y is outcome, X is exposure

arguments :

`β_Y` : vector of outcome effect sizes
`se_β_Y` : vector of 
"""
function mr_wm(β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
               se_β_Y::AbstractVector{<: Union{AbstractFloat, Missing}}, 
               β_X::AbstractVector{<: Union{AbstractFloat, Missing}},
               se_β_X::AbstractVector{<: Union{AbstractFloat, Missing}}, 
               α::AbstractFloat = 0.05;
               iterations::Integer = 5000,
               seed = 42)::mr_output
    
    m = length(β_X)
    if !(length(β_Y) == length(se_β_Y) == m == length(se_β_X)) throw(ArgumentError("β_x, β_y, se_β_X, se_β_Y Vectors must be of same length")) end
    if iterations < 2 throw(ArgumentError("at least two iterations are needed")) end
    if m < 3
        return mr_output(m)
    end
    
    θ_est = _wm_estimate(β_Y, se_β_Y, β_X)
    # bootstrap fo find standard error
    θ_se_est = _bootstrap_se_wm(β_Y, se_β_Y, β_X, se_β_X, iterations, seed)

    dh = Normal(0, θ_se_est)
    dobs = Normal(θ_est, θ_se_est)
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θ_est))

    return mr_output(m, θ_est, θ_se_est, θ_ci_low, θ_ci_high, p)
end

