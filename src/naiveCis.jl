using InMemoryDatasets
using DLMReader
using SnpArrays
using Base.Threads
using ThreadPools
using Statistics


    ########################
    #       Constants      #
    ########################

const mrNamesDict = Dict(mr_egger => "Egger_",
                   mr_ivw => "IVW_",
                   mr_wald => "Wald_ratio_",
                   mr_wm => "Weighted_median_")

const mr_output_fields = fieldnames(mr_output)

const NOut = length(mr_output_fields)

    ########################
    #       NaiveCis       #
    ########################

"""
Implmentation of naive approach for transcriptome wide MR
    Returns a Dataset of results for each exposure (rows) and 
"""
function NaiveCis(data::Union{Dataset, GroupBy}, GenotypesArr::AbstractVector{SnpData}; 
                  one_file_per_chr_plink::Bool = true,
                  r2_tresh::Float64 = 0.1,
                  mr_methodsV::AbstractVector{Function} = [mr_egger, mr_ivw],
                  α::Float64 = 0.05
                  )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        @warn "More than one bimbedfam set of files but not corresponding to chromosomes : expected only 1, only first will be considered"
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # make header
    ncol = NOut*length(mr_methodsV)
    mr_names = Vector{String}(undef, ncol)
    @inbounds for i in 1:ncol
        n = mr_methodsV[div(i, NOut, RoundUp)]
        if haskey(mrNamesDict, n)
            mr_names[i] = mrNamesDict[n]
        else
            mr_names[i] = string(n)
        end
    end

    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methodsV))
    header = ["exposure_name"; mr_names .* fields; "f_min_iv"; "f_max_iv"; "f_med_ivs"]


    if data isa Dataset && !InMemoryDatasets.index(data).grouped[]
        throw(ArgumentError("Expected a grouped Dataset or GroupBy type."))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methodsV))
    f_arr = Array{Float64}(undef, length(eachgroup(data)), 3)
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    for (i, data_group) in enumerate(eachgroup(data))
        
        # prioritize ivs by lowest pvalue
        ivs_d = sort(data_group, :pval_exp)

        # if one plink fileset per chromosome, take file correponfing to exposure chromosome
        if one_file_per_chr_plink
            # clump list of potential ivs to only keep independant
            kept_v_b = clump(GenotypesArr[ivs_d.chr[1]], 
                             Vector{Tuple{Int8, Int32}}(collect(zip(ivs_d.chr, ivs_d.pos))), 
                             r2_tresh = r2_tresh, 
                             formated = true)
        else # else take first
            # clump list of potential ivs to only keep independant
            kept_v_b = clump(GenotypesArr[1], 
                             Vector{Tuple{Int8, Int32}}(collect(zip(ivs_d.chr, ivs_d.pos))), 
                             r2_tresh = r2_tresh,
                             formated = true)
        end

        ivs_d = ivs_d[kept_v_b, :]

        # harmonisation of effect sizes
        @threads for row in axes(ivs_d, 1)
            @inbounds if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
                ivs_d[row, :β_out] = - ivs_d[row, :β_out]
                ivs_d[row, :a_effect_out], ivs_d[row, :a_other_out] = ivs_d[row, :a_other_out], ivs_d[row, :a_effect_out]
            end
        end
        
        β_out, β_exp, se_out, se_exp = ivs_d.β_out.val, ivs_d.β_exp.val, ivs_d.se_out.val, ivs_d.se_exp.val
        f_v = abs2.(β_exp)./abs2.(se_exp)
        if length(f_v) > 0
            f_med = Statistics.median(f_v)
            f_min, f_max = extrema(f_v)
        else
            f_med = f_min = f_max = NaN
        end
        # make mr methods and write ouput in dataset
        @threads for (j, mr_method) in collect(enumerate(mr_methodsV))
            if size(ivs_d, 1) > 0
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields] # badly optimized => think of a rapper of optmized version for users?
        end
        exposureNamesV[i] = data_group.trait[1]
        f_arr[i, :] = [f_min, f_max, f_med] # F stat for ivs (exposure association force)
    end
    
    return Dataset([exposureNamesV outputArr f_arr], header)
end
