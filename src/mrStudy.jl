using DLMReader
using InMemoryDatasets
using Base.Threads
using SnpArrays
using Folds
using Chain
using PooledArrays
using StatsBase

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

# Set bigfloat to store small pvalues with low precision
setprecision(BigFloat, 4)


########################
#     utils funcs      #
########################

# tell if alleles in exposure data and outcome data correspond to a biallelic situation
function biallelic(s::SubArray)::Bool
    return (s[1] == s[2] && s[3] == s[4]) || (s[1] == s[4] && s[2] == s[3])
end

# make a vector of header symbol for QTLStudy/GWAS and a dictionary of types to pass to DLMReader
function make_types_and_headers(file::Union{QTLStudy,GWAS};
    pval_bigfloat::Bool=false,
    reverse::Bool=false)::Tuple{Dict{Int,DataType},Vector{Symbol}}
    if reverse
        outcome_b = file isa QTLStudy
    else
        outcome_b = file isa GWAS
    end
    types::Dict{Int,DataType} = Dict()
    n_cols_file::Int = 0
    GenVarInfo_Types_copy = copy(GenVarInfo_Types)
    if pval_bigfloat
        GenVarInfo_Types_copy[PVAL] = BigFloat
    end

    path_ex = (file isa GWAS) ? file.path : file.path_v[1]

    open(path_ex, "r") do f
        n_cols_file = (file.separator isa AbstractVector) ? count(i -> (i in file.separator), readline(path_ex)) + 1 : count(file.separator, readline(path_ex)) + 1
    end

    header = repeat([:x], n_cols_file)

    for i in 1:n_cols_file
        if haskey(file.columns, i)
            types[i] = GenVarInfo_Types_copy[file.columns[i]]
            header[i] = (outcome_b) ? GenVarInfo_Symbols_out[file.columns[i]] : GenVarInfo_Symbols_exp[file.columns[i]]
        else
            types[i] = String
        end
    end
    return types, header
end

# verify columns don't miss information
# QTL implementation
function verify_columns(study::QTLStudy, is_exposure::Bool=true)

    if is_exposure
        gen_infos = Set((CHR, POS, BETA, SE, A_EFFECT, A_OTHER, PVAL, TRAIT_NAME))
    else
        gen_infos = Set((CHR, POS, BETA, SE, A_EFFECT, A_OTHER, TRAIT_NAME))
    end

    cols_ok = Dict{GenVarInfo,Bool}(CHR => 0,
        POS => 0,
        BETA => 0,
        SE => 0,
        PVAL => 0,
        A_EFFECT => 0,
        A_OTHER => 0,
        TRAIT_NAME => !(nothing in study.traits_for_each_path))

    for value in values(study.columns)
        cols_ok[value] = 1
    end

    for key in gen_infos
        if cols_ok[key] == 0
            throw(ArgumentError("Missing $(key) in columns. QTL Should contain at least : CHR, POS, A_EFFECT, A_OTHER, BETA, SE, PVAL and TRAIT_NAME."))
        end
    end

    all_infos = values(study.columns)
    if length(all_infos) != length(unique(all_infos))
        throw(ArgumentError("There are repeating elements in columns. Make sure to specify only one column number per GenVarInfo type."))
    end
end

# verify files don't miss information
# GWAS implementation
function verify_columns(study::GWAS, is_exposure::Bool=false)
    if is_exposure
        gen_infos = Set((CHR, POS, BETA, SE, A_EFFECT, A_OTHER, PVAL))
    else
        gen_infos = Set((CHR, POS, BETA, SE, A_EFFECT, A_OTHER))
    end

    cols_ok = Dict{GenVarInfo,Bool}(CHR => 0,
        POS => 0,
        BETA => 0,
        SE => 0,
        PVAL => 0,
        A_EFFECT => 0,
        A_OTHER => 0)

    for value in values(study.columns)
        cols_ok[value] = 1
    end

    for key in gen_infos
        if cols_ok[key] == 0
            throw(ArgumentError("Missing $(key) in columns. GWAS Should contain at least : CHR, POS, A_EFFECT, A_OTHER, BETA, SE, PVAL."))
        end
    end

    all_infos = values(study.columns)
    if length(all_infos) != length(unique(all_infos))
        throw(ArgumentError("There are repeating elements in columns. Make sure to specify only one column number per GenVarInfo type."))
    end
end


# read all exposure files and filter for in window and pvalue treshold
function read_filter_file(file::GWAS,
    add_trait_name_b::Bool,
    threads::Bool,
    filtered::Bool,
    trsf_log_pval::Union{Function,Nothing},
    in_window::Function,
    header::Vector{Symbol},
    types::Dict{Int,DataType};
    reverse::Bool=false)::Dataset

    if reverse
        a_effect = :a_effect_out
        a_other = :a_other_out
        pval = :pval_out
    else
        a_effect = :a_effect_exp
        a_other = :a_other_exp
        pval = :pval_exp
    end

    d = filereader(file.path, delimiter=file.separator,
        header=header, types=types, skipto=2,
        makeunique=true, eolwarn=false, threads=threads)[:, collect(keys(file.columns))]

    if add_trait_name_b
        d.trait = repeat([file.trait_name], nrow(d))
    end

    modify!(d, :trait => PooledArray, [a_effect, a_other] .=> (PooledArray ∘ (x -> lowercase.(x))))

    if trsf_log_pval !== nothing
        modify!(d, pval => trsf_log_pval, threads=threads)
    end

    if !filtered
        d = @chain d begin
            filter([:chr, :pos, :trait, pval],
                type=in_window, threads=threads, missings=false) # dataset filtered for window and significance
            filter(a_effect, type=x -> length(x) == 1, threads=threads, missings=false) # remove indels
            filter(a_other, type=x -> length(x) == 1, threads=threads, missings=false)  # remove indels
        end
    end

    return d
end


# read all qtl files and concatenate them into a single Dataset
function read_qtl_files_cis(exposure::QTLStudy,
    types::Dict{Int,DataType},
    header::Vector{Symbol},
    window::Integer,
    p_tresh::Float64,
    filtered::Bool=false,
    trsf_log_pval::Union{Function,Nothing}=nothing)::Dataset

    # Dicionary of trait -- tss
    ref_dict = Dict(zip(exposure.trait_v, zip(exposure.chr_v, exposure.tss_v)))

    # boolan tells if the variant is significant causal on exposure and if in window arround good tss
    # equivalent to : good chr && in 500kb window && pval lower than threshold
    in_window(s::SubArray) = (
        haskey(ref_dict, s[3]) &&
        s[1] == ref_dict[s[3]][1] && # good chr
        abs(s[2] - ref_dict[s[3]][2]) ≤ window && # in window
        s[4] < p_tresh # p_thresh significant
    )

    data_vect = Vector{Dataset}(undef, length(exposure))

    add_trait_name_b = !(TRAIT_NAME in values(exposure.columns))

    for i in 1:lastindex(data_vect)
        file = exposure[i]
        data_vect[i] = read_filter_file(file, add_trait_name_b, true, filtered, trsf_log_pval, in_window, header, types)
    end

    data_filtered = Folds.reduce(vcat, data_vect, init=Dataset())

    return data_filtered
end

function read_qtl_files_trans(exposure::QTLStudy,
    types::Dict{Int,DataType},
    header::Vector{Symbol},
    p_tresh::Float64,
    filtered::Bool=false,
    trsf_log_pval::Union{Function,Nothing}=nothing;
    reverse::Bool=false)::Dataset

    data_vect = Vector{Dataset}(undef, length(exposure))

    filter_fun(s::SubArray) = s[4] < p_tresh

    add_trait_name_b = !(TRAIT_NAME in values(exposure.columns))

    for i in 1:lastindex(data_vect)
        file = exposure[i]
        data_vect[i] = read_filter_file(file, add_trait_name_b, true, filtered, trsf_log_pval, filter_fun, header, types, reverse=reverse)
    end

    data_filtered = Folds.reduce(vcat, data_vect, init=Dataset())

    return data_filtered
end

# Remove all pleiotropic variants and apply final selection filters
function apply_MiLoP!(joined_d, exposure, window, p_tresh)
    joined_d.chr_pos = collect(zip(joined_d.chr, joined_d.pos))
    counts = countmap(joined_d.chr_pos)
    is_unique_iv(s) = counts[s] == 1
    filter!(joined_d, :chr_pos, by=is_unique_iv)
    filter!(joined_d, :pval_exp, by=(<(p_tresh)))

    ref_dict = Dict(zip(exposure.trait_v, zip(exposure.chr_v, exposure.tss_v)))

    in_window(s::SubArray) = (
        haskey(ref_dict, s[3]) &&
        s[1] == ref_dict[s[3]][1] && # good chr
        abs(s[2] - ref_dict[s[3]][2]) ≤ window # in window
    )
    filter!(joined_d, [:chr, :pos, :trait], type = in_window)
end

#########################
#     mrStudyNFolds     #
#########################

# low memory version valid if one file per exposure in QTLStudy
"""
Perform Mendelian Randomization study from QTL to GWAS and separating exposure in n folds to avoid loading full dataset in memory. 
    (works if separated in multiple files only)

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome GWAS data \\
`type::AbstractString` : specifies if the analysis is a trans or cis MR analysis (is either `\"trans\"` or `\"cis\"`)\\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\


**options :**

mrStudyNFolds contains the same options as [`mrStudy`](@ref)\\
`n_folds::Integer` : Number f folds to separate QTl into.

see [`mrStudy`](@ref) for other options and examples.
"""
function mrStudyNFolds(exposure::QTLStudy,
    outcome::GWAS,
    type::AbstractString,
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    n_folds=10,
    approach::AbstractString="naive",
    p_tresh::Float64=1e-3,
    window::Int=500000,
    r2_tresh::Float64=0.1,
    exposure_filtered=false,
    mr_methods::AbstractVector{Function}=[mr_egger, mr_ivw],
    α::Float64=0.05,
    trsf_pval_exp::Union{Function,Nothing}=nothing,
    trsf_pval_out::Union{Function,Nothing}=nothing,
    pval_bigfloat::Bool=false,
    write_ivs::Union{AbstractString,Nothing}=nothing,
    write_filtered_exposure::Union{AbstractString,Nothing}=nothing,
    filter_beta_ratio::Real=0,
    min_maf::Real=0,
    infos::Bool=true
)::Union{Dataset,GroupBy}
    if approach != "naive"
        throw(ArgumentError("aproach should not be MiLoP with mrStudyNFolds."))
    end

    arr_d = Vector{Dataset}([])
    for qtl in nfolds(exposure, n_folds)
        push!(arr_d, mrStudy(qtl,
            outcome,
            type,
            bedbimfam_dirnames,
            approach=approach,
            p_tresh=p_tresh,
            p_tresh_MiLoP=p_tresh,
            window=window,
            r2_tresh=r2_tresh,
            exposure_filtered=exposure_filtered,
            mr_methods=mr_methods,
            α=α,
            trsf_pval_exp=trsf_pval_exp,
            trsf_pval_out=trsf_pval_out,
            pval_bigfloat=pval_bigfloat,
            write_ivs=write_ivs,
            write_filtered_exposure=write_filtered_exposure,
            filter_beta_ratio = filter_beta_ratio,
            min_maf=min_maf,
            infos=infos))
    end

    return Folds.reduce(vcat, arr_d, init=Dataset())
end

"""
Perform Mendelian Randomization study from GWAS to QTL and separating outcome in n folds to avoid loading full dataset in memory. 
    (works if separated in multiple files only)

**arguments :**

`exposure::GWAS` : exposure GWAS data \\
`outcome::QTLStudy` : outcome QTL data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\


**options :**

mrStudyNFolds contains the same options as [`mrStudy`](@ref)\\
`n_folds::Integer` : Number f folds to separate QTl into.

see [`mrStudy`](@ref) for other options and examples.
"""
function mrStudyNFolds(exposure::GWAS,
    outcome::QTLStudy,
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    n_folds=10,
    approach::String="naive",
    p_tresh::Float64=1e-3,
    r2_tresh::Float64=0.1,
    mr_methods::AbstractVector{Function}=[mr_egger, mr_ivw],
    α::Float64=0.05,
    trsf_pval_exp::Union{Function,Nothing}=nothing,
    trsf_pval_out::Union{Function,Nothing}=nothing,
    pval_bigfloat::Bool=false,
    write_ivs::Union{AbstractString,Nothing}=nothing,
    filter_beta_ratio::Real=0,
    min_maf::Real=0,
    infos::Bool=true
)::Union{Dataset,GroupBy}

    if approach != "naive"
        throw(ArgumentError("aproach should not be MiLoP with mrStudyNFolds."))
    end

    arr_d = Vector{Dataset}([])
    for qtl in nfolds(outcome, n_folds)
        push!(arr_d, mrStudy(exposure,
            qtl,
            bedbimfam_dirnames,
            approach=approach,
            p_tresh=p_tresh,
            r2_tresh=r2_tresh,
            mr_methods=mr_methods,
            α=α,
            trsf_pval_exp=trsf_pval_exp,
            trsf_pval_out=trsf_pval_out,
            pval_bigfloat=pval_bigfloat,
            write_ivs=write_ivs,
            filter_beta_ratio=filter_beta_ratio,
            min_maf=min_maf,
            infos=infos))
    end

    return Folds.reduce(vcat, arr_d, init=Dataset())
end

#########################
#        mrStudy        #
#########################

"""
Perform a Mendelian Randomization study with exposure QTL and outcome GWAS

**arguments :**

`exposure::QTLStudy` : exposure QTL data \\
`outcome::GWAS` : outcome gas data \\
`type::AbstractString` : specifies if the analysis is a trans or cis MR analysis (is either `\"trans\"` or `\"cis\"`)\\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\

**options :**\\

`approach::String`: name of MR study aproach chosen (either naive, test or MiLoP) (default is "naive")\\
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
`low_ram::Bool` : If true, if the exposure files each contain only one exposure trait, [`mrStudyNFolds`](@ref) with `n_folds` of 10 will be used.\\
`write_ivs::AbstractString` : write selected Instrumental variables to specified directory\\
`write_filtered_exposure::AbstractString` : write a filtered version of exposure files to specified file name.
    This file will be tab separated and will only contain columns necessary for further MR Studies.\\
`filter_beta_raio::Real` : Filter IVs for which the exposure effect is `filter_beta_raio` times outcome effect or greater. default is 0.\\
`pval_bigfloat::Bool` : use `true` if pvalues can be under `5e-324`. (default is `false`)\\
`min_maf::Real` : minimal variant maf to be kept as a potential IV. (default is 0)\\
`infos::Bool` : If true, infos about advancement compute are printed to terminal (default is `true`)

## Examples

```julia
results = mrStudy(qtl, gwas, "trans", genotypes, 10, approach = "naive", window = 250000, trsf_pval_exp = x -> exp10.(x))
```
""" 
function mrStudy(exposure::QTLStudy,
    outcome::GWAS,
    type::AbstractString,
    bedbimfam_dirnames::AbstractArray{<:AbstractString};

    # options
    approach::AbstractString="naive",
    p_tresh::Float64=1e-3,
    p_tresh_MiLoP::Float64=p_tresh,
    window::Integer=500000,
    r2_tresh::Float64=0.1,
    exposure_filtered::Bool=false,
    mr_methods::AbstractVector{Function}=[mr_egger, mr_ivw],
    α::Float64=0.05,
    trsf_pval_exp::Union{Function,Nothing}=nothing,
    trsf_pval_out::Union{Function,Nothing}=nothing,
    pval_bigfloat::Bool=false,
    write_ivs::Union{Nothing,AbstractString}=nothing,
    write_filtered_exposure::Union{AbstractString,Nothing}=nothing,
    filter_beta_ratio::Real=0,
    min_maf::Real=0,
    infos::Bool=true
)::Union{Dataset,GroupBy}

    # input validity verification
    if approach ∉ ["MiLoP", "naive", "test", "test-MiLoP"]
        throw(ArgumentError("approach must be either : \"MiLoP\", \"naive\", \"test\" or \"test-MiLoP\""))
    end
    if type ∉ ["trans", "cis"]
        throw(ArgumentError("type must be either \"trans\" or \"cis\""))
    end

    # load and filter qtl data (filter for significan snps to exposure and within specified window)
    if infos
        @info "reading and filtering exposure files..."
    end

    # Verify input format and generate types and names for each column of the exposure dataset
    verify_columns(exposure)
    types, header = make_types_and_headers(exposure; pval_bigfloat=pval_bigfloat)

    # Read exposure according to structure defined by `types` and `header` in trans or cis.
    if type == "cis" && approach ∈ ["naive", "test"]
        qtl_d = read_qtl_files_cis(exposure, types, header, window, p_tresh_MiLoP, exposure_filtered, trsf_pval_exp)
    elseif type == "cis" # If MiLoP, widen window to catch all pleiotropic effeccts
        qtl_d = read_qtl_files_cis(exposure, types, header, 250_000_000, p_tresh_MiLoP, exposure_filtered, trsf_pval_exp)
    else
        qtl_d = read_qtl_files_trans(exposure, types, header, p_tresh_MiLoP, exposure_filtered, trsf_pval_exp)
    end

    # If the user specified to write a filtered version of the exposure to speedup future analysis
    if write_filtered_exposure !== nothing
        filewriter(write_filtered_exposure, qtl_d, delimiter='\t')
    end

    if infos
        @info "reading outcome file..."
    end

    # verify input format and generate types and names for each column of the exposure dataset
    verify_columns(outcome)
    types, header = make_types_and_headers(outcome)

    # Read the outcome data according to `header` and `types`
    gwas_d = filereader(outcome.path,
        delimiter=outcome.separator,
        header=header, types=types, skipto=2,
        makeunique=true, eolwarn=false)[:, collect(keys(outcome.columns))]

    # Transform pval if trsf pval is not nothing
    if trsf_pval_out !== nothing
        modify!(gwas_d, :pval => trsf_pval_out)
    end
    # change gwas a_effect_out a_other_out to a Pooled array to save memory
    modify!(gwas_d, [:a_effect_out, :a_other_out] => (PooledArray ∘ x -> lowercase.(x)))

    if size(gwas_d, 1) == 0 || size(qtl_d, 1) == 0
        @warn "No IVS were found with given parameters."
        return Dataset()
    end

    if infos
        @info "Joining exposure and outcome datasets..."
    end

    # join datasets and filter out non biallelic alleles (alleles from both datasets do not correspond)
    joined_d = @chain qtl_d begin
        innerjoin(gwas_d, on=[:chr, :pos], makeunique=false, check=false)
        # keep only biallelic snps
        filter([:a_effect_exp, :a_effect_out, :a_other_exp, :a_other_out], type=biallelic, missings=false) # filter for "obvious" non biallelic variants
        filter(:, by=!ismissing)
    end

    # if MiLoP remove all redundant snps (associated to more than one exposure)
    # Than apply final filters
    if approach == "MiLoP" || approach == "test-MiLoP"
        apply_MiLoP!(joined_d, exposure, window, p_tresh)
    end

    # filter out variants for which the association to the outcome is has larger effect than association to the exposure if option enabled
    if filter_beta_ratio > 0
        beta_compare_b(s) = abs(s[2]) / abs(s[1]) ≤ filter_beta_ratio
        filter!(joined_d, [:β_exp, :β_out], type=beta_compare_b, missings=false)
    end

    if infos
        @info "loading genotypes..."
    end
    # load and format reference snp data
    GenotypesArr = Vector{SnpData}(undef, length(bedbimfam_dirnames))
    @threads for i in 1:lastindex(bedbimfam_dirnames)
        GenotypesArr[i] = SnpData(SnpArrays.datadir(bedbimfam_dirnames[i]))
        formatSnpData!(GenotypesArr[i])
    end

    one_file_per_chr_plink = length(bedbimfam_dirnames) > 1

    if infos
        @info "grouping traits..."
    end

    groupby!(joined_d, :trait, stable=false)

    if infos
        @info "Performing clumping and MR..."
    end

    # For each exposure, select independant IVs and perform MR analysis 
    if approach == "naive" || approach == "MiLoP"
        return ClumpAndMR(joined_d, GenotypesArr, r2_tresh=r2_tresh,
            one_file_per_chr_plink=one_file_per_chr_plink,
            mr_methods=mr_methods, α=α,
            write_ivs=write_ivs, min_maf=min_maf)
    else
        return joined_d
    end
end


"""
Perform a Trans-Mendelian Randomization study with exposure GWAS and outcome GWAS

**arguments :**

`exposure::GWAS` : exposure QTL data \\
`outcome::GWAS` : outcome GWAS data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\

**options :** \\

`approach::String`: name of MR study aproach chosen (either naive, test or MiLoP) (default is "naive")\\
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
`low_ram::Bool` : If true, if the exposure files each contain only one exposure trait, [`mrStudyNFolds`](@ref) with `n_folds` of 10 will be used.
`write_ivs::AbstractString` : write selected Instrumental variables to specified directory\\
`write_filtered_exposure::AbstractString` : write a filtered version of exposure files to specified file name.
    This file will be tab separated and will only contain columns necessary for further MR Studies.\\
`pval_bigfloat::Bool` : use `true` if pvalues can be under `5e-324`. (default is `false`)\\
`filter_beta_raio::Real` : Filter IVs for which the exposure effect is `filter_beta_raio` times outcome effect or greater. default is 0.\\
`infos::Bool` : If true, infos about advancement compute are printed to terminal (default is `true`)

## Examples

```julia
results = mrStudy(gwas1, gwas2, genotypes, 10, approach = "naive", trsf_pval_exp = x -> exp10.(x))
```
"""
function mrStudy(exposure::GWAS,
    outcome::GWAS,
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    approach::String="naive",
    p_tresh::Float64=1e-3,
    r2_tresh::Float64=0.1,
    exposure_filtered::Bool=false,
    mr_methods::AbstractVector{Function}=[mr_egger, mr_ivw],
    α::Float64=0.05,
    trsf_pval_exp::Union{Function,Nothing}=nothing,
    trsf_pval_out::Union{Function,Nothing}=nothing,
    low_ram::Bool=false, # temporary? if as performant --> set true as default value
    pval_bigfloat::Bool=false,
    write_ivs::Union{AbstractString,Nothing}=nothing,
    write_filtered_exposure::Union{AbstractString,Nothing}=nothing,
    filter_beta_ratio::Real=0,
    min_maf::Real=0,
    infos::Bool=true
)::Union{Dataset,GroupBy}

    exposure_name = (exposure.trait_name === nothing) ? "exposure" : exposure.trait_name

    qtl_exposure = QTLStudy(exposure.path, [exposure_name], [exposure_name], nothing, nothing, exposure.columns, exposure.separator)

    return mrStudy(qtl_exposure, outcome, "trans",
        bedbimfam_dirnames;
        approach=approach,
        p_tresh=p_tresh,
        p_tresh_MiLoP=p_tresh,
        r2_tresh=r2_tresh,
        exposure_filtered=exposure_filtered.
        mr_methods = mr_methods,
        α=α,
        trsf_pval_exp=trsf_pval_exp,
        trsf_pval_out=trsf_pval_out,
        low_ram=low_ram,
        pval_bigfloat=pval_bigfloat,
        write_ivs=write_ivs,
        write_filtered_exposure=write_filtered_exposure,
        filter_beta_raio=filter_beta_ratio,
        min_maf=min_maf,
        infos=infos)

end


"""
Perform a Trans-Mendelian Randomization study with exposure GWAS and outcome QTL

**arguments :**

`exposure::GWAS` : exposure GWAS data \\
`outcome::QTLStudy` : outcome QTL data \\
`bedbimfam_dirnames::AbstractArray{<:AbstractString}` : base names of Plink bedbimfam files for reference genotypes 
    (see [SnpArrays documentation](https://openmendel.github.io/SnpArrays.jl/latest/))\\

**options :** 

`approach::String`: name of MR study aproach chosen (either naive, test or MiLoP) (default is "naive")\\
`p_tresh::Float64`: pvalue threshold for a SNP to be considered associated to an exposure (default is 1e-3)\\
`r2_tresh::Float64`: maximial corrlation between to SNPs (default is 0.1)\\
`mr_methods::AbstractVector{Function}` : Functions to use to estimate effect of exposure on outcome.
    Any Function taking four vectors of same length (βoutcome, se_outcome, βexposure, se_exposure) and a Float (α) 
    and returns a value of type [`mr_output`](@ref) can be used, that includes user defined functions. 
    Functions already implemented in this module include [`mr_ivw`](@ref), [`mr_egger`](@ref), [`mr_wm`](@ref) and [`mr_wald`](@ref). default value is `[mr_ivw, mr_egger]` \\
`α::Float64` : α value for confidance intervals of parameter estimations in MR (e.g. 95% CI is α = 0.05, which is the default value)\\
`trsf_pval_exp::Union{Function, Nothing}` : Transformation to apply to pvalues in exposure dataset\\
`trsf_pval_out::Union{Function, Nothing}` : t = Transormation to apply on pvalues in outcome dataset\\
`low_ram::Bool` : If true, if the exposure files each contain only one exposure trait, [`mrStudyNFolds`](@ref) with `n_folds` of 10 will be used.
`write_ivs::AbstractString` : write selected Instrumental variables to specified directory\\
`pval_bigfloat::Bool` : use `true` if pvalues can be under `5e-324`. (default is `false`)\\
`filter_beta_raio::Real` : Filter IVs for which the exposure effect is `filter_beta_raio` times outcome effect or greater. default is 0.\\
`infos::Bool` : If true, infos about advancement compute are printed to terminal (default is `true`)

"""
function mrStudy(exposure::GWAS,
    outcome::QTLStudy,
    bedbimfam_dirnames::AbstractArray{<:AbstractString};
    approach::String="naive",
    p_tresh::Float64=1e-3,
    r2_tresh::Float64=0.1,
    mr_methods::AbstractVector{Function}=[mr_egger, mr_ivw],
    α::Float64=0.05,
    trsf_pval_exp::Union{Function,Nothing}=nothing,
    trsf_pval_out::Union{Function,Nothing}=nothing,
    low_ram::Bool=false, # temporary? if as performant --> set true as default value
    pval_bigfloat::Bool=false,
    write_ivs::Union{AbstractString,Nothing}=nothing,
    filter_beta_ratio::Real=0,
    min_maf::Real=0,
    infos::Bool=true
)

    # input validity verification
    if approach ∉ ["naive", "test"]
        throw(ArgumentError("approach must be either : naive, test"))
    end

    #load and filter qtl data (filter for significan snps to exposure and within specified window)
    if infos
        @info "reading outcome..."
    end
    verify_columns(outcome, false)
    types, header = make_types_and_headers(outcome; reverse=true)
    qtl_d = read_qtl_files_trans(outcome, types, header, 1.0, true, trsf_pval_out, reverse=true)
    qtl_d = @chain qtl_d begin
        filter(:a_effect_out, type=x -> length(x) == 1, missings=false) # remove indels
        filter(:a_other_out, type=x -> length(x) == 1, missings=false) # remove indels
    end
    # load gwas data
    if infos
        @info "reading and filtering exposure..."
    end
    verify_columns(exposure, true)
    types, header = make_types_and_headers(exposure, pval_bigfloat=pval_bigfloat, reverse=true)
    gwas_d = filereader(exposure.path,
        delimiter=exposure.separator,
        header=header, types=types, skipto=2,
        makeunique=true, eolwarn=false)[:, collect(keys(exposure.columns))]
    # Transform pval if trsf pval is not nothing
    if trsf_pval_exp !== nothing
        modify!(gwas_d, :pval_exp => trsf_pval_exp)
    end
    # change gwas a_effect_out a_other_out to a Pooled array to save memory
    modify!(gwas_d, [:a_effect_exp, :a_other_exp] => (PooledArray ∘ x -> lowercase.(x)))

    gwas_d = @chain gwas_d begin
        filter(:pval_exp, type=x -> x .< p_tresh, missings=false)
        filter(:a_effect_exp, type=x -> length(x) == 1, missings=false) # remove indels
        filter(:a_other_exp, type=x -> length(x) == 1, missings=false)
    end

    ## Do not join if one Dataset is empty
    if size(gwas_d, 1) == 0 || size(qtl_d, 1) == 0
        @warn "No IVS were found with given parameters."
        return Dataset()
    end

    if infos
        @info "Joining exposure and outcome datasets..."
    end

    joined_d = @chain gwas_d begin
        innerjoin(qtl_d, on=[:chr, :pos], makeunique=false, check=false)
        # keep only biallelic snps
        filter([:a_effect_exp, :a_effect_out, :a_other_exp, :a_other_out], type=biallelic, missings=false) # filter for "obvious" non biallelic variants
        filter(:, by=!ismissing)
    end

    if filter_beta_ratio > 0
        beta_compare_b(s) = abs(s[2]) / abs(s[1]) ≤ filter_beta_ratio
        filter!(joined_d, [:β_exp, :β_out], type=beta_compare_b, missings=false)
    end

    if infos
        @info "loading genotypes..."
    end

    # load and format reference snp data
    GenotypesArr = Vector{SnpData}(undef, length(bedbimfam_dirnames))
    @threads for i in 1:lastindex(bedbimfam_dirnames)
        GenotypesArr[i] = SnpData(SnpArrays.datadir(bedbimfam_dirnames[i]))
        formatSnpData!(GenotypesArr[i])
    end

    one_file_per_chr_plink = length(bedbimfam_dirnames) > 1

    if infos
        @info "grouping by trait..."
    end

    groupby!(joined_d, :trait, stable=false)

    if infos
        @info "performing clumping and MR..."
    end
    # for each exposure, select independant IVs and perform MR
    if approach == "naive"
        return ClumpAndMR(joined_d, GenotypesArr,
            r2_tresh=r2_tresh,
            one_file_per_chr_plink=one_file_per_chr_plink,
            mr_methods=mr_methods,
            α=α,
            write_ivs=write_ivs,
            min_maf=min_maf)
    else
        return joined_d
    end

end
