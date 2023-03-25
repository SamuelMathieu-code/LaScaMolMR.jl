using DataFrames

PLINK = get(ENV, "PLINK_PATH", Sys.BINDIR) # Get PLINK_PATH as enc variable default being the environement's bin directory.


##################################################
#                    ivSelect                    #
##################################################


"""
Find potential ivs from exposure with cis method
Returns tuple made of DataFrame of ivs (chr, pos, effect_allele, other_allele, β, se, pval)
We assume that gwas arguments trait_name, chr, tss are not nothing.
"""
function get_exp_ivs_cis(exposure::GWAS, window::Int64 = 500000, p_threshold::Float64 = 1e-8)::DataFrame  # NOTE : Besoin de vérifier les alleles?
    fp = open(exposure.path)
    iterable = eachline(fp)
    data = DataFrame([[], [], [], [], [], [], []], ["chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval"])
    
    for line in iterable
        trait, variant, effect = exposure.acc_func(line)
        if ((trait==exposure.trait_name || isnothing(trait)) &&
            variant[1]==exposure.chr &&
            variant[2]>=exposure.tss-window &&
            variant[2]<=exposure.tss+window &&
            effect[3] <= p_threshold
        )
            push!(data, [transpose(variant) transpose(effect)])
        end
    end

    close(fp)

    return data
end


"""
Find potential ivs from exposure with trans method
"""
function get_exp_ivs_trans(exposure::GWAS)::DataFrame
    #...
    return [], []
end


"""
Merge exp iv list with snps available in target + clumping of ivs if corr_treshold is not `nothing` --> PLINK.
"""
function get_corresp_ivs(target::GWAS, pot_ivs::Array{Int64, 2},  corr_threshold::Union{Nothing, Float16} = 0.1)
    #...
    return [], []
end


"""
Find ivs for MR from `exposure` to `target` using cis method and with maximal correlation treshold between ivs of `corr_threshold`.
If `corr_threshold` is `nothing` : than no clumping is done. 
"""
function ivSelectCis(exposure::GWAS, target::GWAS, corr_threshold::Union{Nothing, Float16} = 0.1)::Tuple{Array{Int64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    pot_ivs = get_exp_ivs_cis(exposure)
    ivs = get_corresp_ivs(target, pot_ivs, corr_threshold)
    return ivs
end


"""
Find ivs for MR from `exposure` to `target` using trans method and with maximal correlation treshold between ivs of `corr_threshold`.
If `corr_threshold` is `nothing` : than no clumping is done. 
"""
function ivSelectTrans(exposure::GWAS, target::GWAS, corr_threshold::Union{Nothing, Float16} = 0.1)::Tuple{Array{Int64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    pot_ivs = get_exp_ivs_trans(exposure)
    ivs = get_corresp_ivs(target, pot_ivs, corr_threshold)
    return ivs
end
