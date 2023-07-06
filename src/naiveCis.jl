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
                   mr_wald => "Wald_ratio_")

const mr_output_fields = fieldnames(mr_output)

const NOut = length(mr_output_fields)

    ########################
    #       NaiveCis       #
    ########################

"""
Implmentation of naive approach for transcriptome wide MR
    Returns a Dataset of results for each exposure (rows) and 
"""
function NaiveCis(data::GroupBy, GenotypesArr::AVSD where AVSD <: AbstractVector{SnpData}; 
                  one_file_per_chr_plink = true,
                  r2_tresh::Float64 = 0.1,
                  mr_methodsV::AbstractVector{Function} = [mr_egger, mr_ivw],
                  α::Float64 = 0.05
                  )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        throw(ArgumentError("More than one bimbedfam set of files but not corresponding to chromosomes : expected"))
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methodsV))
    # f_arr = Array{Float64}(undef, length(eachgroup(data)), 3)
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    for (i, data_group) in enumerate(eachgroup(data))
        
        ivs_d = sort(data_group, :pval_exp)
        
        # if one plink fileset per chromosome, take file correponfing to exposure chromosome
        # How could we optimize this if else statement?
        if one_file_per_chr_plink
            kept_v_b = clump(GenotypesArr[ivs_d.chr[1]], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        else # else take first
            kept_v_b = clump(GenotypesArr[1], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        end

        ivs_d = ivs_d[kept_v_b, :]

        # harmonisation -----------> Would it be better to do the harmonisation at the same time as biallelic filtering (in mrStudyCis)? If so : how?
        @threads for row in axes(ivs_d, 1) # could be more efficient? If so : how?
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
                ivs_d[row, :β_out] = - ivs_d[row, :β_out]
                ivs_d[row, :a_effect_out], ivs_d[row, :a_other_out] = ivs_d[row, :a_other_out], ivs_d[row, :a_effect_out]
            end
        end
        
        β_out, β_exp, se_out, se_exp = ivs_d.β_out.val, ivs_d.β_exp.val, ivs_d.se_out.val, ivs_d.se_exp.val
        # f_v = abs2.(β_exp)./abs2.(se_exp)
        # f_med = Statistics.median(f_v)
        # f_min, f_max = extrema(f_v)
        # make mr methods and write ouput in dataset
        @threads for (j, mr_method) in collect(enumerate(mr_methodsV))
            if size(ivs_d, 1) ≥ 1
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields] # badly optimized => think of a rapper of optmized version for users?
        end
        exposureNamesV[i] = data_group.trait[1]
        # f_arr[i, :] = [f_min, f_max, f_med] # F stat for ivs (exposure association force)
    end
    mr_names = [mrNamesDict[mr_methodsV[div(x, NOut, RoundUp)]] for x in 1:(NOut * length(mr_methodsV))]
    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methodsV))
    # header = ["exposure_name"; mr_names .* fields; "f_min_iv"; "f_max_iv"; "f_med_ivs"]
    header = ["exposure_name"; mr_names .* fields]
    
    return Dataset([exposureNamesV outputArr], header)
end

function NaiveCis2(data::GroupBy, GenotypesArr::AVSD where AVSD <: AbstractVector{SnpData}; 
                    one_file_per_chr_plink = true,
                    r2_tresh::Float64 = 0.1,
                    mr_methodsV::AbstractVector{Function} = [mr_egger, mr_ivw],
                    α::Float64 = 0.05
                    )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        throw(ArgumentError("More than one bimbedfam set of files but not corresponding to chromosomes : expected"))
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methodsV))
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    @threads for (i, data_group) in collect(enumerate(eachgroup(data)))
        
        ivs_d = sort(data_group, :pval_exp, threads = false)
        
        # if one plink fileset per chromosome, take file correponfing to exposure chromosome
        # How could we optimize this if else statement?
        if one_file_per_chr_plink
            kept_v_b = clump2(GenotypesArr[ivs_d.chr[1]], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        else # else take first
            kept_v_b = clump2(GenotypesArr[1], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        end

        ivs_d = ivs_d[kept_v_b, :]

        # harmonisation -----------> Would it be better to do the harmonisation at the same time as biallelic filtering (in mrStudyCis)? If so : how?
        for row in axes(ivs_d, 1) # could be more efficient? If so : how?
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
                ivs_d[row, :β_out] = - ivs_d[row, :β_out]
                ivs_d[row, :a_effect_out], ivs_d[row, :a_other_out] = ivs_d[row, :a_other_out], ivs_d[row, :a_effect_out]
            end
        end
        
        β_out, β_exp, se_out, se_exp = ivs_d.β_out.val, ivs_d.β_exp.val, ivs_d.se_out.val, ivs_d.se_exp.val
        # make mr methods and write ouput in dataset
        for (j, mr_method) in enumerate(mr_methodsV)
            if size(ivs_d, 1) ≥ 1
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields] # badly optimized => think of a rapper of optmized version for users?
        end
        exposureNamesV[i] = data_group.trait[1]
    end
    mr_names = [mrNamesDict[mr_methodsV[div(x, NOut, RoundUp)]] for x in 1:(NOut * length(mr_methodsV))]
    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methodsV))
    header = ["exposure_name"; mr_names .* fields]
    
    return Dataset([exposureNamesV outputArr], header)
end


function NaiveCis3(data::GroupBy, GenotypesArr::AVSD where AVSD <: AbstractVector{SnpData}; 
    one_file_per_chr_plink = true,
    r2_tresh::Float64 = 0.1,
    mr_methodsV::AbstractVector{Function} = [mr_egger, mr_ivw],
    α::Float64 = 0.05
    )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        throw(ArgumentError("More than one bimbedfam set of files but not corresponding to chromosomes : expected"))
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methodsV))
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    @qthreads for (i, data_group) in collect(enumerate(eachgroup(data)))

        ivs_d = sort(data_group, :pval_exp, threads = false)

        # if one plink fileset per chromosome, take file correponfing to exposure chromosome
        # How could we optimize this if else statement?
        if one_file_per_chr_plink
            kept_v_b = clump2(GenotypesArr[ivs_d.chr[1]], 
                        collect(zip(ivs_d.chr, ivs_d.pos)), 
                        r2_tresh)
        else # else take first
            kept_v_b = clump2(GenotypesArr[1], 
                        collect(zip(ivs_d.chr, ivs_d.pos)), 
                        r2_tresh)
        end

        ivs_d = ivs_d[kept_v_b, :]

        # harmonisation -----------> Would it be better to do the harmonisation at the same time as biallelic filtering (in mrStudyCis)? If so : how?
        for row in axes(ivs_d, 1) # could be more efficient? If so : how?
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
                ivs_d[row, :β_out] = - ivs_d[row, :β_out]
                ivs_d[row, :a_effect_out], ivs_d[row, :a_other_out] = ivs_d[row, :a_other_out], ivs_d[row, :a_effect_out]
            end
        end

        β_out, β_exp, se_out, se_exp = ivs_d.β_out.val, ivs_d.β_exp.val, ivs_d.se_out.val, ivs_d.se_exp.val
        # make mr methods and write ouput in dataset
        for (j, mr_method) in enumerate(mr_methodsV)
            if size(ivs_d, 1) ≥ 1
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields] # badly optimized => think of a rapper of optmized version for users?
        end
        exposureNamesV[i] = data_group.trait[1]
    end
    mr_names = [mrNamesDict[mr_methodsV[div(x, NOut, RoundUp)]] for x in 1:(NOut * length(mr_methodsV))]
    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methodsV))
    header = ["exposure_name"; mr_names .* fields]

    return Dataset([exposureNamesV outputArr], header)
end