
using Distributions
using GLM
##################################################
#                     MrPerf                     #
##################################################


"""
Struct encapsulating the outputs of a Mendelian Randomization anaysis
"""
struct mr_output
    nivs::Int
    effect::Float64
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


# """
# Struct encasulating the outputs of Cochan's Q test
# """
# struct cochran_output
#     Q::Float64
#     p::Float64
# end


"""
Wald ratio for Mendelian Randomization with a single instrumentl variable
"""
function mr_wald(β_Y::Float64, 
                 se_β_Y::Float64, 
                 β_X::Float64,  
                 α::Float64 = 0.05)::mr_output  # À véerifier que la distribustion normale convient!!!!!!!!
    
    θ = β_Y / β_X
    se_θ = se_β_Y / abs(β_X)
    dh = Normal(0, se_θ)
    dobs = Normal(θ, se_θ)
    p = 2*cdf(dh, -abs(θ))
    ci_low, ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    
    return mr_output(1, θ, ci_low, ci_high, p, NaN, NaN, NaN, NaN, NaN, NaN)
end


"""
Inverse variance weighted linear regression with simple weights (se(B_Y)^-2) Mendelian Randomization
"""
function mr_ivw(β_Y::Vector{Float64}, 
                se_β_Y::Vector{Float64}, 
                β_X::Vector{Float64}, 
                α::Float64 = 0.05)::mr_output 

    # regression
    regressor = lm(@formula(β_Y ~ 0 + β_X), (;β_X, β_Y), wts = se_β_Y .^ (-2))
    θivw_est = coef(regressor)[1]
    se_θivw_est = stderror(regressor)[1]

    dh = Normal(0, se_θivw_est)
    dobs = Normal(θivw_est, se_θivw_est)
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θivw_est))

    # heterogeneity --> À revoir éventuellement
    if length(β_Y) ≥ 2
        heter_stat = (length(β_X)-1) * (1/(length(β_X) - 2))^2 * sum((β_Y - residuals(regressor)).^2) # (n-1) * rse^2
        chisq = Chisq(length(β_X) - 1) 
        heter_p = 1 - cdf(chisq, heter_stat)
    else
        heter_stat = heter_p = NaN
    end

    return mr_output(length(β_Y), θivw_est, θ_ci_low, θ_ci_high, p, NaN, NaN, NaN, NaN, heter_stat, heter_p) 

end


# """
# Weighted Median linear regression Medndelian Randomization --> To come in the future
# """
# function mr_wm(β_Y::Vector{Float64}, 
#                se_β_Y::Vector{Float64}, 
#                β_X::Vector{Float64}, 
#                α::Float64 = 0.05)::mr_output 

#     return mr_output()
# end


"""
Egger Mendelian Randomization
"""
function mr_egger(β_Y::Vector{Float64}, 
                  se_β_Y::Vector{Float64}, 
                  β_X::Vector{Float64}, 
                  α::Float64 = 0.05)::mr_output
    
    # regression
    β_Y_abs = sign.(β_X).*β_Y
    β_X_abs = abs.(β_X)
    regressor = lm(@formula(β_Y_abs ~ β_X_abs), (;β_X_abs, β_Y_abs), wts = se_β_Y .^ (-2))
    θ_est = coef(regressor)
    se_θ_est = stderror(regressor)

    # index 1 is intercept and index two is effect
    dh = Normal(0, se_θ_est[2])
    dobs = Normal(θ_est[2], se_θ_est[2])
    dintercept = Normal(0, se_θ_est[1])
    θ_ci_low, θ_ci_high = quantile(dobs, α/2), quantile(dobs, 1-α/2)
    p = 2*cdf(dh, -abs(θ_est[2]))
    p_intercept = 1 - 2*(cdf(dintercept, -abs(θ_est[1])))
    ci_low_intercept, ci_high_intercept = quantile(dintercept, α/2), quantile(dintercept, 1-α/2)

    # heterogeneity
    if length(β_Y) ≥ 3
        heter_stat = sum((residuals(regressor) ./ se_β_Y).^2)  # Pourquoi c'est différent de 
        chisq = Chisq(length(β_X) - 2)
        heter_p = 1 - cdf(chisq, heter_stat)
    else
        heter_stat = heter_p = NaN
    end

    return mr_output(length(β_Y), θ_est[2], θ_ci_low, θ_ci_high, p, 
                     θ_est[1], p_intercept, ci_low_intercept, ci_high_intercept, 
                     heter_stat, heter_p) 

end


# """
# Cochran's Q test for heterogeneity
# """
# function cochran(β_Y::Vector{Float64}, 
#                  se_β_Y::Vector{Float64}, 
#                  β_X::Vector{Float64}, 
#                  α::Float64 = 0.05)::cochran_output 

#     return cochran_output()
# end