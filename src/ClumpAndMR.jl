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


"""
Keep independant IVs and perform MR in Omic-wide MR (Trans or Cis)
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
`α:AbstractFloat` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)\\
`write_ivs::AbstractString` : write selected Instrumental variables to specified directory

## Examples :

```julia
res_d = clumpAndMR(d, GenotypesArr, r2_tresh = 0.1)
```
"""
function clumpAndMR(data::Union{Dataset, GroupBy}, GenotypesArr::AbstractVector{SnpData}; 
                  one_file_per_chr_plink::Bool = true,
                  r2_tresh::AbstractFloat = 0.1,
                  mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
                  α::AbstractFloat = 0.05,
                  write_ivs::Union{Nothing, AbstractString} = nothing,
                  min_maf::Real = 0
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
    @inbounds for (i, data_group) in enumerate(eachgroup(data))
        
        ivs_d = sort(data_group, :pval_exp) #prioritize snps by pval
        ivs_d[!, :kept] .= true
        groupby!(ivs_d, :chr, stable = true)

        for ivs_in_chr in eachgroup(ivs_d) #clump different chromosomes separately
            
            if one_file_per_chr_plink
                g = GenotypesArr[ivs_in_chr.chr[1]]
            else # else take first
                g = GenotypesArr[1]
            end

            ivs_in_chr.kept .= clump(g, 
                                 collect(zip(ivs_d.chr.val, ivs_d.pos.val)), 
                                 r2_tresh = r2_tresh, 
                                 formated = true,
                                 min_maf = min_maf)
        end

        filter!(ivs_d, :kept)
        select!(ivs_d, Not(:kept))

        # harmonisation of effect sizes
        @threads for row in axes(ivs_d, 1)
            if ivs_d[row, :a_effect_exp] != ivs_d[row, :a_effect_out]
                ivs_d[row, :β_out] = - ivs_d[row, :β_out]
                ivs_d[row, :a_effect_out], ivs_d[row, :a_other_out] = ivs_d[row, :a_other_out], ivs_d[row, :a_effect_out]
            end
        end
        
        # Approximate F statistic of intrument~exposure association
        β_out, β_exp, se_out, se_exp = ivs_d.β_out.val, ivs_d.β_exp.val, ivs_d.se_out.val, ivs_d.se_exp.val
        f_v = abs2.(β_exp)./abs2.(se_exp)
        # get median and extremas of F stats
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
            filewriter(joinpath(write_ivs, data_group.trait[1] * "_ivs.txt"), 
                       @view(ivs_d[:, [:trait, :chr, :pos, 
                                       :a_effect_exp, :a_other_exp, :β_exp, :se_exp, :pval_exp, 
                                       :a_effect_out, :a_other_out, :β_out, :se_out, :pval_out]]),
                       delimiter = '\t')
        end

    end

    return Dataset([exposureNamesV outputArr f_arr], header)
end

