

PLINK = get(ENV, "PLINK_PATH", Sys.BINDIR) # Get PLINK_PATH as enc variable default being the environement's bin directory.


##################################################
#                    ivSelect                    #
##################################################


"""
Find potential ivs from exposure with cis method
Returns tuple made of array of ivs (chr, pos) and an array of data (Beta, se, pval)
"""
function get_exp_ivs_cis(exposure::GWAS)::Tuple{Array{Int64, 2}, Array{Float64, 2}}  # NOTE : Faire un accesseur pour avoir le (beta, se, pval) rapidement :
                                                                                     #        nouvel argument du struct GWAS et QtlStudy?
    #...
    return [], []
end


"""
Find potential ivs from exposure with trans method
"""
function get_exp_ivs_trans(exposure::GWAS)::Tuple{Array{Int64, 2}, Array{Float64, 2}}
    #...
    return [], []
end


"""
Merge exp iv list with snps available in target + clumping of ivs if corr_treshold is not `nothing` --> PLINK.
"""
function get_corresp_ivs(target::GWAS, pot_ivs::Array{Int64, 2},  corr_threshold::Union{nothing, Float16} = 0.1)::Tuple{Vector{Int64}, Array{Float64, 2}}
    #...
    return [], []
end


"""
Find ivs for MR from `exposure` to `target` using cis method and with maximal correlation treshold between ivs of `corr_threshold`.
If `corr_threshold` is `nothing` : than no clumping is done. 
"""
function ivSelectCis(exposure::GWAS, target::GWAS, corr_threshold::Union{nothing, Float16} = 0.1)::Tuple{Array{Int64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    pot_ivs, data_exp = get_exp_ivs_cis(exposure)
    index, data_target = get_corresp_ivs(target, pot_ivs, corr_threshold, )
    pot_ivs, data_exp = pot_ivs[index], data_exp[index]
    return pot_ivs, data_exp, data_taget
end


"""
Find ivs for MR from `exposure` to `target` using trans method and with maximal correlation treshold between ivs of `corr_threshold`.
If `corr_threshold` is `nothing` : than no clumping is done. 
"""
function ivSelectTrans(exposure::GWAS, target::GWAS, corr_threshold::Union{nothing, Float16} = 0.1)::Tuple{Array{Int64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    pot_ivs, data_exp = get_exp_ivs_trans(exposure)
    index, data_target = get_corresp_ivs(target, pot_ivs, corr_threshold, )
    pot_ivs, data_exp = pot_ivs[index], data_exp[index]
    return pot_ivs, data_exp, data_taget
end
