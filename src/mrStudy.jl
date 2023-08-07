
using DLMReader
using InMemoryDatasets
using Base.Threads
using SnpArrays
using Folds
using Chain
using PooledArrays

    ########################
    #       Constants      #
    ########################

const GenVarInfo_Types = Dict(TRAIT_NAME => String,
                        CHR => Int8,
                        POS => Int32,
                        A_EFFECT => String,
                        A_OTHER => String,
                        BETA => Float64,
                        SE => Float64,
                        PVAL => Float64)

const GenVarInfo_Symbols_exp = Dict(TRAIT_NAME => :trait,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_exp,
                          A_OTHER => :a_other_exp,
                          BETA => :β_exp,
                          SE => :se_exp,
                          PVAL => :pval_exp)

const GenVarInfo_Symbols_out = Dict(TRAIT_NAME => :trait_out_,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_out,
                          A_OTHER => :a_other_out,
                          BETA => :β_out,
                          SE => :se_out,
                          PVAL => :pval_out)

setprecision(BigFloat, 4)


    ########################
    #     utils funcs      #
    ########################

# make a vector of header symbol for QTLStudy/GWAS and a dictionary of types to pass to DLMReader
function make_types_and_headers(file::Union{QTLStudy, GWAS};
                                pval_bigfloat::Bool = false)::Tuple{Dict{Int, DataType}, Vector{Symbol}}
    types::Dict{Int, DataType} = Dict()
    n_cols_file::Int = 0
    GenVarInfo_Types_copy = copy(GenVarInfo_Types)
    if pval_bigfloat
        GenVarInfo_Types_copy[PVAL] = BigFloat
    end

    path_ex = (file isa GWAS) ? file.path : file.path_v[1]

    open(path_ex, "r") do f
        n_cols_file = (file.separator isa AbstractVector) ? count(i->(i in file.separator), readline(path_ex)) + 1 : count(file.separator, readline(path_ex)) + 1
    end

    header = repeat([:x], n_cols_file)

    for i in 1:n_cols_file
        if haskey(file.columns, i)
            types[i] = GenVarInfo_Types_copy[file.columns[i]]
            header[i] = (file isa GWAS) ? GenVarInfo_Symbols_out[file.columns[i]] : GenVarInfo_Symbols_exp[file.columns[i]]
        else
            types[i] = String
        end
    end
    return types, header
end

# verify columns make sense and simplify repeating elements (modifies the columns field of QTLStudy/GWAS)
function verify_and_simplify_columns!(exposure::Union{QTLStudy, GWAS})

    if exposure isa QTLStudy
        trait_each_path_nothing =  nothing in exposure.traits_for_each_path
        if trait_each_path_nothing && !(TRAIT_NAME in values(exposure.columns))
            throw(ArgumentError("exposure misses TRAIT_NAME information"))
        end
        gen_infos = [CHR, POS, BETA, SE, A_EFFECT, A_OTHER, PVAL, TRAIT_NAME]
    else
        trait_each_path_nothing = false
        gen_infos = [CHR, POS, BETA, SE, A_EFFECT, A_OTHER]
    end

    # Treating exposure
    # keep only columns and verify all columns are satisfied 
    new_cols_exposure::Dict{Int, GenVarInfo} = Dict()
    cols_ok = Dict{GenVarInfo, Bool}(CHR => 0,
                   POS => 0,
                   BETA => 0,
                   SE => 0,
                   PVAL => (exposure isa GWAS) ? 1 : 0,
                   A_EFFECT => 0,
                   A_OTHER => 0,
                   TRAIT_NAME => (!trait_each_path_nothing) ? 1 : 0)

    for key in keys(exposure.columns)
        for info in gen_infos
            if exposure.columns[key] == info
                cols_ok[info] = 1
                new_cols_exposure[key] = info
            end
        end
    end

    for key in keys(cols_ok)
        if cols_ok[key] == 0
            throw(ArgumentError("Missing information in columns. GWAS/QTL Should contain at least : CHR, POS, A_EFFECT, A_OTHER, BETA, SE, PVAL."))
        end
    end
    exposure.columns = new_cols_exposure
end


# read all exposure files and filter for in window and pvalue treshold
function read_filter_file(file::GWAS, 
                          add_trait_name_b::Bool, 
                          threads::Bool, 
                          filtered::Bool, 
                          trsf_log_pval::Union{Function, Nothing}, 
                          in_window::Function, 
                          header::Vector{Symbol}, 
                          types::Dict{Int, DataType})::Dataset
    
    d = filereader(file.path, delimiter = file.separator, 
                    header = header, types = types, skipto=2, 
                    makeunique=true, eolwarn=false, threads = threads)[:,collect(keys(file.columns))]
            
    if add_trait_name_b
        d.trait = repeat([file.trait_name], nrow(d))
    end

    modify!(d, :trait => PooledArray, [:a_effect_exp, :a_other_exp] .=> (PooledArray ∘ (x -> lowercase.(x))))
    
    if trsf_log_pval !== nothing
        modify!(d, :pval_exp => trsf_log_pval, threads = threads)
    end

    if !filtered
        d = @chain d begin 
            filter([:chr, :pos, :trait, :pval_exp], 
                type = in_window, threads = threads, missings = false) # dataset filtered for window and significance
            filter(:a_effect_exp, type = x -> length(x) == 1, threads = threads, missings = false) # remove indels
            filter(:a_other_exp, type = x -> length(x) == 1, threads = threads, missings = false)  # remove indels
        end
    end
    
    return d
end


# read all qtl files and concatenate them into a single Dataset
function read_qtl_files_cis(exposure::QTLStudy,
    types::Dict{Int, DataType}, 
    header::Vector{Symbol},
    window::Integer,
    p_tresh::Float64,
    filtered::Bool =false,
    trsf_log_pval::Union{Function, Nothing} = nothing)::Dataset

    # Dicionary of trait -- tss
    ref_dict = Dict(zip(exposure.trait_v, zip(exposure.chr_v, exposure.tss_v)))

    # boolan tells if the variant is significant causal on exposure and if in window arround good tss
    # equivalent to : good chr && in 500kb window && pval lower than threshold
    in_window(s::SubArray) = (
                              s[1] == ref_dict[s[3]][1] && # good chr
                              abs(s[2]-ref_dict[s[3]][2])≤window && # in window
                              s[4]<p_tresh # p_thresh significant
                            )
    
    data_vect = Vector{Dataset}(undef, length(exposure))

    add_trait_name_b = !(TRAIT_NAME in values(exposure.columns))

    for i in 1:lastindex(data_vect)
        file = exposure[i]
        data_vect[i] = read_filter_file(file, add_trait_name_b, true, filtered, trsf_log_pval, in_window, header, types)
    end

    data_filtered = Folds.reduce(vcat, data_vect, init = Dataset())

    return data_filtered
end

function read_qtl_files_trans(exposure::QTLStudy,
                              types::Dict{Int, DataType}, 
                              header::Vector{Symbol},
                              p_tresh::Float64,
                              filtered::Bool = false,
                              trsf_log_pval::Union{Function, Nothing} = nothing)::Dataset

    data_vect = Vector{Dataset}(undef, length(exposure))

    filter_fun(s::SubArray) = s[4] < p_tresh

    add_trait_name_b = !(TRAIT_NAME in values(exposure.columns))

    for i in 1:lastindex(data_vect)
        file = exposure[i]
        data_vect[i] = read_filter_file(file, add_trait_name_b, true, filtered, trsf_log_pval, filter_fun, header, types)
    end

    data_filtered = Folds.reduce(vcat, data_vect, init = Dataset())

    return data_filtered
end


    #########################
    #      MrStudyCis       #
    #########################

# low memory version valid if one file per exposure in QTLStudy
"""
Perform Cis-Mendelian Randomization study from QTL to GWAS and separating exposure in n folds to avoid loading full dataset in memory. 
    (works if separated in multiple files only) This may be slightly slower than mrStudyCis.

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome gas data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\
`n_folds::Integer` : Number f folds to separate QTl into.

see [`mrStudyCis`](@ref) for options and examples.
"""
function mrStudyCisNFolds(exposure::QTLStudy, 
                           outcome::GWAS, 
                           bedbimfam_dirnames::AbstractArray{<:AbstractString};
                           n_folds = 10,
                           approach::String="naive", 
                           p_tresh::Float64 = 1e-3, 
                           window::Int = 500000, 
                           r2_tresh::Float64 = 0.1,
                           exposure_filtered = false,
                           mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
                           α::Float64 = 0.05,
                           trsf_pval_exp::Union{Function, Nothing} = nothing,
                           trsf_pval_out::Union{Function, Nothing} = nothing,
                           pval_bigfloat::Bool = false,
                           write_ivs::Union{AbstractString, Nothing} = nothing,
                           write_filtered_exposure::Union{AbstractString, Nothing} = nothing
                          )::Union{Dataset, GroupBy}
    if approach != "naive" throw(ArgumentError("aproach should not be strict with mrStudyCisNFolds.")) end

    arr_d = Vector{Dataset}([])
    for qtl in nfolds(exposure, n_folds)
        push!(arr_d, mrStudyCis(qtl, 
                                outcome, 
                                bedbimfam_dirnames,
                                approach = approach, 
                                p_tresh = p_tresh, 
                                window = window, 
                                r2_tresh = r2_tresh,
                                exposure_filtered = exposure_filtered,
                                mr_methods = mr_methods,
                                α = α,
                                trsf_pval_exp = trsf_pval_exp,
                                trsf_pval_out = trsf_pval_out,
                                pval_bigfloat = pval_bigfloat,
                                write_ivs = write_ivs,
                                write_filtered_exposure = write_filtered_exposure))
    end

    return Folds.reduce(vcat, arr_d, init = Dataset())
end


"""
Perform a Cis-Mendelian Randomization study with exposure QTL and outcome GWAS

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome gas data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\

**options : ** \\
`approach::String`: name of MR study aproach chosen (either naive, test or strict) (default is "naive")\\
`p_tresh::Float64`: pvalue threshold for a SNP to be considered associated to an exposure (default is 1e-3)\\
`window::Integer`: maximal distance between a potential Instrument Variable and transciption start site of gene exposure (default is 500_000)\\
`r2_tresh::Float64`: maximial corrlation between to SNPs (default is 0.1)\\
`exposure_filtered::Bool` : If true, the exposure files are considered already filtered will not filtered 
    on distance to tss and level of significance (default is false)\\
`mr_methods::AbstractVector{Function}` : Functions to use to estimate effect of exposure on outcome.
    Any Function taking four vectors of same length (βoutcome, se_outcome, βexposure, se_exposure) and a Float (α) 
    and returns a value of type [`mr_output`](@ref) can be used, that includes user defined functions. 
    Functions already implemented in this module include [`mr_ivw`](@ref), [`mr_egger`](@ref), [`mr_wm`](@ref) and [`mr_wald`](@ref). default value is `[mr_ivw, mr_egger]` \\
`α::Float64` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)\\
`trsf_pval_exp::Union{Function, Nothing}` : Transformation to apply to pvalues in exposure dataset\\
`trsf_pval_out::Union{Function, Nothing}` : t = Transormation to apply on pvalues in outcome dataset\\
`low_ram::Bool` : If true, if the exposure files each contain only one exposure trait, [`mrStudyCisNFolds`](@ref) with n_folds of 10 will be used.

## Examples

```julia
results = mrStudyCis(qtl, gwas, genotypes, 10, approach = "naive", window = 250000, trsf_pval_exp = x -> exp10.(x))
```
""" ### Trouver des meilleurs exemples.
function mrStudyCis(exposure::QTLStudy, 
    outcome::GWAS, 
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    approach::String="naive", 
    p_tresh::Float64 = 1e-3, 
    window::Integer = 500000, 
    r2_tresh::Float64 = 0.1,
    exposure_filtered::Bool = false,
    mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
    α::Float64 = 0.05,
    trsf_pval_exp::Union{Function, Nothing} = nothing,
    trsf_pval_out::Union{Function, Nothing} = nothing,
    low_ram::Bool = false, # temporary? if as performant --> set true as default value,
    pval_bigfloat::Bool = false,
    write_ivs::Union{Nothing, AbstractString} = nothing,
    write_filtered_exposure::Union{AbstractString, Nothing} = nothing
    )::Union{Dataset, GroupBy}

    # input validity verification
    l_unique_traits = length(unique(exposure.traits_for_each_path))
    if approach ∉ ["strict", "naive", "test", "test-strict"] throw(ArgumentError("approach must be either : strict, naive, test or test-strict")) end
    if (approach != "naive" || length(exposure.path_v) < 10) && low_ram 
        @warn "low_ram option in mrStudyCis with approach different from \"naive\" or less than 10 files in QTLStudy will not be considered. 
        Pay attention to Memory state." 
    end

    # verify if one qtl file per exposure and low_ram -> make folds to limit ram usage
    if low_ram && l_unique_traits == length(exposure.path_v) && approach == "naive"
        return mrStudyCisNFolds(exposure, 
                                 outcome, 
                                 bedbimfam_dirnames;
                                 approach = approach, 
                                 p_tresh = p_tresh, 
                                 window = window, 
                                 r2_tresh = r2_tresh,
                                 exposure_filtered = exposure_filtered,
                                 mr_methods = mr_methods,
                                 α = α,
                                 trsf_pval_exp = trsf_pval_exp,
                                 trsf_pval_out = trsf_pval_out,
                                 pval_bigfloat = pval_bigfloat)
    end 
    #load and filter qtl data (filter for significan snps to exposure and within specified window)
    verify_and_simplify_columns!(exposure)
    types, header = make_types_and_headers(exposure; pval_bigfloat = pval_bigfloat)
    qtl_d = read_qtl_files_cis(exposure, types, header, window, p_tresh, exposure_filtered, trsf_pval_exp)

    if write_filtered_exposure !== nothing
        filewriter(write_filtered_exposure, qtl_d, delimiter = '\t')
    end
    
    # load gwas data
    verify_and_simplify_columns!(outcome)
    types, header = make_types_and_headers(outcome)
    gwas_d = filereader(outcome.path, 
                    delimiter = outcome.separator, 
                    header = header, types = types, skipto=2, 
                    makeunique=true, eolwarn=false)[:,collect(keys(outcome.columns))]

    # Transform pval if trsf pval is not nothing
    if trsf_pval_out !== nothing
        modify!(gwas_d, :pval => trsf_pval_out)
    end
    # change gwas a_effect_out a_other_out to a Pooled array to save memory
    modify!(gwas_d, [:a_effect_out, :a_other_out] => (PooledArray ∘ x -> lowercase.(x)))

    # keep only biallelic snps
    biallelic(s::SubArray) = (s[1]==s[2] && s[3] == s[4]) || (s[1] == s[4] && s[2] == s[3])
    joined_d = @chain qtl_d begin
        innerjoin(gwas_d, on = [:chr, :pos], makeunique = false)
        filter([:a_effect_exp, :a_effect_out, :a_other_exp, :a_other_out], type = biallelic, missings = false) # filter for "obvious" non biallelic variants
        filter(:, by = !ismissing)
    end

    # if strict remove all redundant snps (associated to more than one exposure)
    if approach == "strict" || approach == "test-strict"
        joined_d.chr_pos = collect(zip(joined_d.chr, joined_d.pos))
        local counts = countmap(joined_d.chr_pos)
        is_unique_iv(s) = counts[s] == 1
        filter!(joined_d, :chr_pos, by = is_unique_iv)
    end

    # load and format reference snp data
    GenotypesArr = Vector{SnpData}(undef, length(bedbimfam_dirnames))
    @threads for i in 1:lastindex(bedbimfam_dirnames)
        GenotypesArr[i] = SnpData(SnpArrays.datadir(bedbimfam_dirnames[i]))
        formatSnpData!(GenotypesArr[i])
    end

    one_file_per_chr_plink = length(bedbimfam_dirnames) > 1

    groupby!(joined_d, :trait, stable = false)

    #### for d in eachgroup(joined_d) -> Plink + MR (implemented in NaiveCis)
    if approach == "naive" || approach == "strict"
        return NaiveCis(joined_d, GenotypesArr, r2_tresh = r2_tresh, one_file_per_chr_plink = one_file_per_chr_plink, mr_methodsV = mr_methods, α = α, write_ivs = write_ivs)
    else
        return joined_d
    end
end




    #########################
    #     MrStudyTrans      #
    #########################


    # low memory version valid if one file per exposure in QTLStudy
"""
Perform Trans-Mendelian Randomization study from QTL to GWAS and separating exposure in n folds to avoid loading full dataset in memory. 
    (works if separated in multiple files only) This may be slightly slower than mrStudyCis.

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome gas data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\
`n_folds::Integer` : Number f folds to separate QTl into.

see [`mrStudyCis`](@ref) for options and examples.
"""
function mrStudyTransNFolds(exposure::QTLStudy, 
                           outcome::GWAS, 
                           bedbimfam_dirnames::AbstractArray{<:AbstractString};
                           n_folds = 10,
                           approach::String="naive", 
                           p_tresh::Float64 = 1e-3, 
                           r2_tresh::Float64 = 0.1,
                           exposure_filtered = false,
                           mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
                           α::Float64 = 0.05,
                           trsf_pval_exp::Union{Function, Nothing} = nothing,
                           trsf_pval_out::Union{Function, Nothing} = nothing,
                           pval_bigfloat::Bool = false,
                           write_ivs::Union{AbstractString, Nothing} = nothing,
                           write_filtered_exposure::Union{AbstractString, Nothing} = nothing,
                           filter_beta_ratio::Real = 0
                          )::Union{Dataset, GroupBy}
    if approach != "naive" throw(ArgumentError("aproach should not be strict with mrStudyCisNFolds.")) end

    arr_d = Vector{Dataset}([])
    for qtl in nfolds(exposure, n_folds)
        push!(arr_d, mrStudyTrans(qtl, 
                                outcome, 
                                bedbimfam_dirnames,
                                approach = approach, 
                                p_tresh = p_tresh,  
                                r2_tresh = r2_tresh,
                                exposure_filtered = exposure_filtered,
                                mr_methods = mr_methods,
                                α = α,
                                trsf_pval_exp = trsf_pval_exp,
                                trsf_pval_out = trsf_pval_out,
                                pval_bigfloat = pval_bigfloat,
                                write_ivs = write_ivs,
                                write_filtered_exposure = write_filtered_exposure,
                                filter_beta_ratio = filter_beta_ratio))
    end

    return Folds.reduce(vcat, arr_d, init = Dataset())
end


"""
Perform a Trans-Mendelian Randomization study with exposure QTL and outcome GWAS

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome gas data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\

**options : ** \\
`approach::String`: name of MR study aproach chosen (either naive, test or strict) (default is "naive")\\
`p_tresh::Float64`: pvalue threshold for a SNP to be considered associated to an exposure (default is 1e-3)\\
`r2_tresh::Float64`: maximial corrlation between to SNPs (default is 0.1)\\
`exposure_filtered::Bool` : If true, the exposure files are considered already filtered will not filtered 
    on distance to tss and level of significance (default is false)\\
`mr_methods::AbstractVector{Function}` : Functions to use to estimate effect of exposure on outcome.
    Any Function taking four vectors of same length (βoutcome, se_outcome, βexposure, se_exposure) and a Float (α) 
    and returns a value of type [`mr_output`](@ref) can be used, that includes user defined functions. 
    Functions already implemented in this module include [`mr_ivw`](@ref), [`mr_egger`](@ref), [`mr_wm`](@ref) and [`mr_wald`](@ref). default value is `[mr_ivw, mr_egger]` \\
`α::Float64` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)\\
`trsf_pval_exp::Union{Function, Nothing}` : Transformation to apply to pvalues in exposure dataset\\
`trsf_pval_out::Union{Function, Nothing}` : t = Transormation to apply on pvalues in outcome dataset\\
`low_ram::Bool` : If true, if the exposure files each contain only one exposure trait, [`mrStudyCisNFolds`](@ref) with n_folds of 10 will be used.

## Examples

```julia
results = mrStudyTrans(mqtl, gwas, genotypes, 10, approach = "naive", trsf_pval_exp = x -> exp10.(x))
```
""" ### Trouver des meilleurs exemples.
function mrStudyTrans(exposure::QTLStudy, 
    outcome::GWAS, 
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    approach::String="naive", 
    p_tresh::Float64 = 1e-3, 
    r2_tresh::Float64 = 0.1,
    exposure_filtered::Bool = false,
    mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
    α::Float64 = 0.05,
    trsf_pval_exp::Union{Function, Nothing} = nothing,
    trsf_pval_out::Union{Function, Nothing} = nothing,
    low_ram::Bool = false, # temporary? if as performant --> set true as default value
    pval_bigfloat::Bool = false,
    write_ivs::Union{AbstractString, Nothing} = nothing,
    write_filtered_exposure::Union{AbstractString, Nothing} = nothing,
    filter_beta_ratio::Real = 0
    )::Union{Dataset, GroupBy}

    # input validity verification
    l_unique_traits = length(unique(exposure.traits_for_each_path))
    if approach ∉ ["strict", "naive", "test", "test-strict"] throw(ArgumentError("approach must be either : strict, naive, test, test-strict")) end
    if (approach != "naive" || length(exposure.path_v) < 10) && low_ram 
        @warn "low_ram option in mrStudyCis with approach different from \"naive\" or less than 10 files in QTLStudy will not be considered. Pay attention to Memory state." 
    end

    # verify if one qtl file per exposure and low_ram -> make folds to limit ram usage
    if low_ram && l_unique_traits == length(exposure.path_v) && approach == "naive"
        return mrStudyTransNFolds(exposure, 
                                 outcome, 
                                 bedbimfam_dirnames;
                                 approach = approach, 
                                 p_tresh = p_tresh, 
                                 r2_tresh = r2_tresh,
                                 exposure_filtered = exposure_filtered,
                                 mr_methods = mr_methods,
                                 α = α,
                                 trsf_pval_exp = trsf_pval_exp,
                                 trsf_pval_out = trsf_pval_out,
                                 pval_bigfloat = pval_bigfloat,
                                 filter_beta_ratio = filter_beta_ratio)
    end 
    #load and filter qtl data (filter for significan snps to exposure and within specified window)
    verify_and_simplify_columns!(exposure)
    types, header = make_types_and_headers(exposure; pval_bigfloat = pval_bigfloat)
    qtl_d = read_qtl_files_trans(exposure, types, header, p_tresh, exposure_filtered, trsf_pval_exp)

    if write_filtered_exposure !== nothing
        filewriter(write_filtered_exposure, qtl_d, delimiter = '\t')
    end
    
    # load gwas data
    verify_and_simplify_columns!(outcome)
    types, header = make_types_and_headers(outcome)
    gwas_d = filereader(outcome.path, 
                        delimiter = outcome.separator, 
                        header = header, types = types, skipto=2, 
                        makeunique=true, eolwarn=false)[:,collect(keys(outcome.columns))]

    # Transform pval if trsf pval is not nothing
    if trsf_pval_out !== nothing
        modify!(gwas_d, :pval => trsf_pval_out)
    end
    # change gwas a_effect_out a_other_out to a Pooled array to save memory
    modify!(gwas_d, [:a_effect_out, :a_other_out] => (PooledArray ∘ x -> lowercase.(x)))

    # keep only biallelic snps
    biallelic(s::SubArray) = (s[1] == s[2] && s[3] == s[4]) || (s[1] == s[4] && s[2] == s[3])
    joined_d = @chain qtl_d begin
        innerjoin(gwas_d, on = [:chr, :pos], makeunique = false)
        filter([:a_effect_exp, :a_effect_out, :a_other_exp, :a_other_out], type = biallelic, missings = false) # filter for "obvious" non biallelic variants
        filter(:, by = !ismissing)
    end

    # if strict remove all redundant snps (associated to more than one exposure)
    if approach == "strict" || approach == "test-strict"
        joined_d.chr_pos = collect(zip(joined_d.chr, joined_d.pos))
        local counts = countmap(joined_d.chr_pos)
        is_unique_iv(s) = counts[s] == 1
        filter!(joined_d, :chr_pos, by = is_unique_iv)
    end

    if filter_beta_ratio > 0
        beta_compare_b(s) = abs(s[2]) ≤ abs(s[1])
        filter!(joined_d, [:β_exp, :β_out], type = beta_compare_b, missings = false)
    end

    # load and format reference snp data
    GenotypesArr = Vector{SnpData}(undef, length(bedbimfam_dirnames))
    @threads for i in 1:lastindex(bedbimfam_dirnames)
        GenotypesArr[i] = SnpData(SnpArrays.datadir(bedbimfam_dirnames[i]))
        formatSnpData!(GenotypesArr[i])
    end

    one_file_per_chr_plink = length(bedbimfam_dirnames) > 1

    groupby!(joined_d, :trait, stable = false)

    #### for d in eachgroup(joined_d) -> Plink + MR (implemented in NaiveCis)
    if approach == "naive" || approach == "strict"
        return NaiveTrans(joined_d, GenotypesArr, r2_tresh = r2_tresh, one_file_per_chr_plink = one_file_per_chr_plink, mr_methodsV = mr_methods, α = α, write_ivs = write_ivs)
    else
        return joined_d
    end
end


function mrStudyTrans(exposure::GWAS, 
                      outcome::GWAS, 
                      bedbimfam_dirnames::AbstractArray{<:AbstractString};
                      approach::String="naive", 
                      p_tresh::Float64 = 1e-3, 
                      r2_tresh::Float64 = 0.1,
                      exposure_filtered::Bool = false,
                      mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
                      α::Float64 = 0.05,
                      trsf_pval_exp::Union{Function, Nothing} = nothing,
                      trsf_pval_out::Union{Function, Nothing} = nothing,
                      low_ram::Bool = false, # temporary? if as performant --> set true as default value
                      pval_bigfloat::Bool = false,
                      write_ivs::Union{AbstractString, Nothing} = nothing,
                      write_filtered_exposure::Union{AbstractString, Nothing} = nothing,
                      filter_beta_ratio::Real = 0
                      )::Union{Dataset, GroupBy}

exposure_name = (exposure.trait_name === nothing) ? "exposure" : exposure.trait_name

qtl_exposure = QTLStudy(exposure.path, [exposure_name], [exposure_name], nothing, nothing, exposure.columns, exposure.separator)

return mrStudyTrans(qtl_exposure, outcome, bedbimfam_dirnames;
                    approach = approach,
                    p_tresh = p_tresh,
                    r2_tresh = r2_tresh,
                    exposure_filtered = exposure_filtered.
                    mr_methods = mr_methods,
                    α = α,
                    trsf_pval_exp = trsf_pval_exp,
                    trsf_pval_out = trsf_pval_out,
                    low_ram = low_ram,
                    pval_bigfloat = pval_bigfloat,
                    write_ivs = write_ivs,
                    write_filtered_exposure = write_filtered_exposure,
                    filter_beta_raio = filter_beta_ratio)

end
