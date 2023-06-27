using Distributions
using GLM
using LinearAlgebra
using Random
using Statistics


    ###########################
    #        Utilities        #
    ###########################

"""
Raw weighted median effect size estimate
"""
function wm_estimate(β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                     se_β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                     β_X::AbstractVector{F} where F <: Union{Float64, Missing}
                     )::Float64
θ_v = β_Y ./ β_X
order = sortperm(θ_v)
θ_v, se_βy_ordered, βx_ordered = θ_v[order], se_β_Y[order], β_X[order]
w_v = (βx_ordered ./ se_βy_ordered) .^ 2
p_v = accumulate(+, w_v) - (w_v ./ 2)
p_v /= sum(w_v)
p = findfirst(x -> x ≥ 0.5, p_v)
θ_est = (p != 1) ? θ_v[p-1] + (θ_v[p] - θ_v[p-1])*(0.5 - p_v[p-1])/(p_v[p] - p_v[p-1]) : 0
return θ_est
end


"""
Numerical esimation of weighted median effect size standard error
"""
function bootstrap_se_wm(β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
    se_β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
    β_X::AbstractVector{F} where F <: Union{Float64, Missing},
    se_β_X::AbstractVector{F} where F <: Union{Float64, Missing}, 
    iterations = 5000,
    seed = 42)::Float64

    Random.seed!(seed)

    dx = MvNormal(Vector{Float64}(β_X), Vector{Float64}(se_β_X))
    dy = MvNormal(Vector{Float64}(β_Y), Vector{Float64}(se_β_Y))
    θ_est_v = Vector{Float64}(undef, iterations)
    for i in 1:iterations
        θ_est_v[i] = wm_estimate(rand(dy), se_β_Y, rand(dx))
    end
    return Statistics.std(θ_est_v)
end


    #########################################
    #                 MrPerf                #
    #########################################

# eventuality : make lm more efficient by not calculating sumstats (only estimate)


"""
Struct encapsulating the outputs of a Mendelian Randomization analysis
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
Wald ratio for Mendelian Randomization with a single instrumental variable
     y is outcome, x is exposure
"""
function mr_wald(β_Y::F where F <: AbstractFloat, 
                 se_β_Y::F where F <: AbstractFloat, 
                 β_X::F where F <: AbstractFloat, 
                 α::F where F <: AbstractFloat = 0.05)::mr_output # À véerifier que la distribustion normale convient!!!!!!!!
    θ = β_Y / β_X
    se_θ = se_β_Y / abs(β_X)
    dh = Normal(0, se_θ)
    dobs = Normal(θ, se_θ)
    p = 2*cdf(dh, -abs(θ))
    ci_low, ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    
    return mr_output(1, θ, se_θ, ci_low, ci_high, p, NaN, NaN, NaN, NaN, NaN, NaN)
end

function mr_wald(β_y::AbstractVector{F} where F <: Union{AbstractFloat, Missing}, 
                 se_β_y::AbstractVector{F} where F <: Union{AbstractFloat, Missing}, 
                 β_x::AbstractVector{F} where F <: Union{AbstractFloat, Missing}, 
                 se_β_X::AbstractVector{F} where F <: Union{Float64, Missing} = [],
                 α::F where F <: AbstractFloat = 0.05)::mr_output
    return mr_wald(β_y[1], se_β_y[1], β_x[1], α)
end


"""
Inverse variance weighted linear regression with simple weights (se(B_Y)^-2) Mendelian Randomization
    Y is outcome, X is exposure
"""
function mr_ivw(β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                se_β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                β_X::AbstractVector{F} where F <: Union{Float64, Missing}, 
                se_β_X::AbstractVector{F} where F <: Union{Float64, Missing} = [],
                α::F where F <: AbstractFloat = 0.05)::mr_output #where F <: Union{AbstractFloat, Missing}

    
    m = length(β_X)
    if m < 2
        return mr_output(m, NaN, NaN, NaN, NaN, NaN, 
                         NaN, NaN, NaN, NaN, 
                         NaN, NaN)
    end
    # regression
    w = se_β_Y .^ (-2)
    regressor = lm(@formula(β_Y ~ 0 + β_X), (;β_X, β_Y), wts = w)
    θivw_est = coef(regressor)[1]
    
    ϵ = residuals(regressor)
    u = sqrt.(w).*ϵ
    X = sqrt.(w).*β_X
    se_θivw_est = sqrt((u'*u)*(X'*X)^(-1)/(m-1))
    σ = sqrt(sum(u.^2)/(m-1))
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
"""
function mr_egger(β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                  se_β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
                  β_X::AbstractVector{F} where F <: Union{Float64, Missing}, 
                  se_β_X::AbstractVector{F} where F <: Union{Float64, Missing} = [],
                  α::F where F <: AbstractFloat = 0.05)::mr_output
    
    # regression
    m = length(β_X)
    if m < 3
        return mr_output(m, NaN, NaN, NaN, NaN, NaN, 
                        NaN, NaN, NaN, NaN, 
                        NaN, NaN)
    end
    w =  se_β_Y .^ (-2)
    β_Y_abs = sign.(β_X).*β_Y
    β_X_abs = abs.(β_X)
    regressor = lm(@formula(β_Y_abs ~ β_X_abs), (;β_X_abs, β_Y_abs), wts = w)
    θ_est = coef(regressor)

    ϵ = residuals(regressor)
    X = sqrt.(w) .* [ones(m) β_X_abs]
    u = sqrt.(w) .* ϵ
    se = (u'*u)*(X'*X)^(-1)/(m-2)
    σ = sqrt(sum(u.^2)/(m-2))
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
    heter_stat = sum((ϵ ./ se_β_Y).^2)
    chisq = Chisq(m - 2)
    heter_p = 1 - cdf(chisq, heter_stat)

    return mr_output(length(β_Y), θ_est[2], se_θ_est[2], θ_ci_low, θ_ci_high, p, 
                     θ_est[1], p_intercept, ci_low_intercept, ci_high_intercept, 
                     heter_stat, heter_p) 

end


"""
Weighted Median Mendelian Randomization (Bowden et al., 2015)
    Y is outcome, X is exposure
"""
function mr_wm(β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
               se_β_Y::AbstractVector{F} where F <: Union{Float64, Missing}, 
               β_X::AbstractVector{F} where F <: Union{Float64, Missing},
               se_β_X::AbstractVector{F} where F <: Union{Float64, Missing}, 
               α::F where F <: AbstractFloat = 0.05;
               iterations::Integer = 5000,
               seed = 42)::mr_output
    
    m = length(β_X)
    if m < 3
        return mr_output(m, NaN, NaN, NaN, NaN, NaN, 
                        NaN, NaN, NaN, NaN, 
                        NaN, NaN)
    end
    
    θ_est = wm_estimate(β_Y, se_β_Y, β_X)
    # bootstrap fo find standard error
    θ_se_est = bootstrap_se_wm(β_Y, se_β_Y, β_X, se_β_X, iterations, seed)

    dh = Normal(0, θ_se_est)
    dobs = Normal(θ_est, θ_se_est)
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θ_est))

    return mr_output(m, θ_est, θ_se_est, θ_ci_low, θ_ci_high, p, NaN, NaN, NaN, NaN, NaN, NaN)
end

