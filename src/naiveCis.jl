using InMemoryDatasets
using DLMReader
using SnpArrays
using Base.Threads
using Iterators


    ########################
    #       Constants      #
    ########################

const mrNamesDict = Dict(mr_egger => "Egger",
                   mr_ivw => "IVW",
                   mr_wald => "Wald_ratio")

macro NOut()
    return :(11)
end

    ########################
    #       NaiveCis       #
    ########################

"""
Implmentation of naive approach for transcriptome wide MR
    Returns a Dataset of results for each exposure (rows) and 
"""
function NaiveCis(data::Dataset, GenotypesArr::AbstractVector{SnpData}, 
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
    outputArr = Array{mr_output}(undef, length(eachgroup(data)), @NOut * length(mr_methodsV))
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    @threads for (i, data_group) in enumerate(eachgroup(data))
        
        ivs_d = sort(data_group, :pval_exp)

        # !!!!!!!!!!!!!!!!!!!!!!!!
        #  ADD HARMONISATION HERE   -> + filter for non biallelic/indels
        # !!!!!!!!!!!!!!!!!!!!!!!!
        
        # if one plink fileset per chromosome, take file correponfing to exposure chromosome
        if one_file_per_chr_plink
            local chr = ivs_d.chr[0]
            kept_v_b = clump(GenotypesArr[chr], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        else # else take first
            kept_v_b = clump(GenotypesArr[1], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        end

        ivs_d = ivs_d[kept_v_b, :]

        
        # make mr methods and write ouput in dataset
        
        for (j, mr_method) in enumerate(mr_methodsV)
            if size(ivs_d, 1) >= 1
                res = mr_method(ivs_d.β_out, ivs_d.se_out, ivs_d.β_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (@NOut * (j - 1) + 1):(@NOut * j)] = [getproperty(res, field) for field in fieldnames(mr_output)] # badly optimized => think of a rapper of optmized version for users?
        end
        exposureNamesV[i] = data_group.trait[1]
    end
    mr_names = [mrNamesDict[mr_methodsV[div(x, @NOut, RoundUp)]] for x in 1:(@NOut * length(mr_methodsV))]
    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methodsV))
    header = ["exposure_name"; mr_names .* fields]
    
    return Dataset([exposureNamesV outputArr], header)
end