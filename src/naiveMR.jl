using InMemoryDatasets
using DLMReader
using SnpArrays
using Base.Threads
using Statistics
using DLMReader

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
    Returns a Dataset of results for each exposure (rows).

**arguments :** 
`data` : grouped dataset (exposure trait) with specific column names : 
    `:trait`, `:chr`, `:pos`, `:β_exp`, `β_out`, `:se_exp`, `:se_out` among others. \\
`GenotypesArr::AbstractVector{SnpData}` : Reference genotypes. (see [`SnpArrays`](https://github.com/OpenMendel/SnpArrays.jl))

**options :**

`one_file_per_chr_plink::Bool` : If true, it is considered that the each index in GenotypesArr corresponds to a chromosome number (default is true) \\
`r2_tresh::AbstractFloat` : The maximal r² correlation coefficient for two snps to be considered independant. (default is 0.1) \\
`mr_methods::AbstractVector{Function}`: Functions to use to estimate effect of exposure on outcome.
    Any Function taking four vectors of same length (βoutcome, se_outcome, βexposure, se_exposure) and a Float (α) 
    and returns a value of type [`mr_output`](@ref) can be used, that includes user defined functions. 
    Functions already implemented in this module include [`mr_ivw`](@ref), [`mr_egger`](@ref), [`mr_wm`](@ref)
    and [`mr_wald`](@ref). default value is `[mr_ivw, mr_egger]` \\
`α:AbstractFloat` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)

## Examples :

```julia
res_d = NaiveCis(d, GenotypesArr, r2_tresh = 0.1)
```
"""
function NaiveCis(data::Union{Dataset, GroupBy}, GenotypesArr::AbstractVector{SnpData}; 
                  one_file_per_chr_plink::Bool = true,
                  r2_tresh::AbstractFloat = 0.1,
                  mr_methods::AbstractVector{Function} = [mr_ivw, mr_egger],
                  α::AbstractFloat = 0.05,
                  write_ivs::Union{AbstractString, Nothing} = nothing
                  )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        @warn "More than one bimbedfam set of files but not corresponding to chromosomes : expected only 1, only first will be considered"
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # make header for final result Dataset
    ncol = NOut*length(mr_methods)
    mr_names = Vector{String}(undef, ncol)
    @inbounds for i in 1:ncol
        n = mr_methods[div(i, NOut, RoundUp)]
        if haskey(mrNamesDict, n)
            mr_names[i] = mrNamesDict[n]
        else
            mr_names[i] = string(n)
        end
    end

    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methods))
    header = ["exposure_name"; mr_names .* fields; "f_min_iv"; "f_max_iv"; "f_med_ivs"]


    if data isa Dataset && !InMemoryDatasets.index(data).grouped[]
        throw(ArgumentError("Expected a grouped Dataset or GroupBy type."))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methods))
    f_arr = Array{Float64}(undef, length(eachgroup(data)), 3)
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    for (i, data_group) in enumerate(eachgroup(data))

        # prioritize ivs by lowest pvalue
        ivs_d = sort(data_group, :pval_exp, view = true)

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
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
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
        @threads for (j, mr_method) in collect(enumerate(mr_methods))
            if size(ivs_d, 1) > 0
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields]
        end

        exposureNamesV[i] = data_group.trait[1]
        f_arr[i, :] = [f_min, f_max, f_med] # F stat for ivs (exposure association force)

        if write_ivs !== nothing
            filewriter(joinpath(write_ivs, data_group.trait[1] * "_ivs.txt"), ivs_d, delimiter = '\t')
        end
    end

    return Dataset([exposureNamesV outputArr f_arr], header)
end


"""
Implmentation of naive approach for Trans Omic-wide MR
    Returns a Dataset of results for each exposure (rows).

**arguments :** 
`data` : grouped dataset (exposure trait) with specific column names : 
    `:trait`, `:chr`, `:pos`, `:β_exp`, `β_out`, `:se_exp`, `:se_out` among others. \\
`GenotypesArr::AbstractVector{SnpData}` : Reference genotypes. (see [`SnpArrays`](https://github.com/OpenMendel/SnpArrays.jl)) formated with [`formatSnpData!`](@ref).

**options :**

`one_file_per_chr_plink::Bool` : If true, it is considered that the each index in GenotypesArr corresponds to a chromosome number (default is true) \\
`r2_tresh::AbstractFloat` : The maximal r² correlation coefficient for two snps to be considered independant. (default is 0.1) \\
`mr_methods::AbstractVector{Function}`: Functions to use to estimate effect of exposure on outcome.
    Any Function taking four vectors of same length (βoutcome, se_outcome, βexposure, se_exposure) and a Float (α) 
    and returns a value of type [`mr_output`](@ref) can be used, that includes user defined functions. 
    Functions already implemented in this module include [`mr_ivw`](@ref), [`mr_egger`](@ref), [`mr_wm`](@ref)
    and [`mr_wald`](@ref). default value is `[mr_ivw, mr_egger]` \\
`α:AbstractFloat` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)

## Examples :

```julia
res_d = NaiveTrans(d, GenotypesArr, r2_tresh = 0.1)
```
"""
function NaiveTrans(data::Union{Dataset, GroupBy}, GenotypesArr::AbstractVector{SnpData}; 
                  one_file_per_chr_plink::Bool = true,
                  r2_tresh::AbstractFloat = 0.1,
                  mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
                  α::AbstractFloat = 0.05,
                  write_ivs::Union{Nothing, AbstractString} = nothing
                  )::Dataset

    # Gestion of bedbimfam file sets
    if !one_file_per_chr_plink && length(GenotypesArr) != 1
        @warn "More than one bimbedfam set of files but not corresponding to chromosomes : expected only 1, only first will be considered"
    elseif one_file_per_chr_plink && (length(GenotypesArr) > 24 || length(GenotypesArr) < 22)
        throw(ArgumentError("One bedbimfam file set per chromosome, expected between 22 and 24 file sets"))
    end

    # make header for final result Dataset
    ncol = NOut*length(mr_methods)
    mr_names = Vector{String}(undef, ncol)
    @inbounds for i in 1:ncol
        n = mr_methods[div(i, NOut, RoundUp)]
        if haskey(mrNamesDict, n)
            mr_names[i] = mrNamesDict[n]
        else
            mr_names[i] = string(n)
        end
    end

    fields = repeat(collect(string.(fieldnames(mr_output))), length(mr_methods))
    header = ["exposure_name"; mr_names .* fields; "f_min_iv"; "f_max_iv"; "f_med_ivs"]


    if data isa Dataset && !InMemoryDatasets.index(data).grouped[]
        throw(ArgumentError("Expected a grouped Dataset or GroupBy type."))
    end

    # for outputs
    outputArr = Array{Float64}(undef, length(eachgroup(data)), NOut * length(mr_methods))
    f_arr = Array{Float64}(undef, length(eachgroup(data)), 3)
    exposureNamesV = Vector{String}(undef, length(eachgroup(data)))

    # MR pipeline for each exposure
    for (i, data_group) in enumerate(eachgroup(data))
        
        # prioritize ivs by lowest pvalue
        ivs_d = sort(data_group, :pval_exp, view = true)

        ### separate into chrs and clump in diff chrs.
        d_each_chr = Dict{Int8, Vector{Int64}}([i => Vector{Int64}([]) for i in 1:22])
        for i in axes(ivs_d, 1)
            push!(d_each_chr[ivs_d[i, :chr]], i)
        end

        ###################################### A verif si ca marche
        all_kept_indx = Vector{Int64}([])
        for i in 1:22

            # if one plink fileset per chromosome, take file correponfing to exposure chromosome
            if one_file_per_chr_plink
                g = GenotypesArr[i]
            else # else take first
                g = GenotypesArr[1]
            end

            # clump for only independant ivs
            ivs_in_chr_idx = d_each_chr[i]
            kept_v_b = clump(g, 
                             Vector{Tuple{Int8, Int32}}(collect(zip(ivs_d[ivs_in_chr_idx, :chr], ivs_d[ivs_in_chr_idx, :pos]))), 
                             r2_tresh = r2_tresh, 
                             formated = true)

            append!(all_kept_indx, ivs_in_chr_idx[kept_v_b])

            
        end

        ivs_d = ivs_d[all_kept_indx, :]
        ################################# Fin A verif si ca marche

        # harmonisation of effect sizes
        @threads for row in axes(ivs_d, 1)
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
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
        @threads for (j, mr_method) in collect(enumerate(mr_methods))
            if size(ivs_d, 1) > 0
                res = mr_method(β_out, se_out, β_exp, se_exp, α)
            else
                res = mr_output(0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
            end
            outputArr[i, (NOut * (j - 1) + 1):(NOut * j)] = [getproperty(res, field) for field in mr_output_fields]
        end

        exposureNamesV[i] = data_group.trait[1]
        f_arr[i, :] = [f_min, f_max, f_med] # F stat for ivs (exposure association force)

        if write_ivs !== nothing
            filewriter(joinpath(write_ivs, data_group.trait[1] * "_ivs.txt"), ivs_d, delimiter = '\t')
        end

    end

    return Dataset([exposureNamesV outputArr f_arr], header)
end

