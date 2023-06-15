using DLMReader
using InMemoryDatasets
using Base.Threads
using SnpArrays
using Folds
using Chain # unused

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
                        PVAL => Float64,
                        LOG_PVAL => Float64)

const GenVarInfo_Symbols_exp = Dict(TRAIT_NAME => :trait,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_exp,
                          A_OTHER => :a_other_exp,
                          BETA => :β_exp,
                          SE => :se_exp,
                          PVAL => :pval_exp,
                          LOG_PVAL => :pval_exp)

const GenVarInfo_Symbols_out = Dict(TRAIT_NAME => :trait_out_,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_out,
                          A_OTHER => :a_other_out,
                          BETA => :β_out,
                          SE => :se_out,
                          PVAL => :pval_out,
                          LOG_PVAL => :pval_out)


    ########################
    #     utils funcs      #
    ########################

function make_types_and_headers(file::Union{QTLStudy, GWAS})::Tuple{Dict{Int, DataType}, Vector{Symbol}}
    types::Dict{Int, DataType} = Dict()
    n_cols_file::Int = 0

    path_ex = (file isa GWAS) ? file.path : file.path_v[1]

    open(path_ex, "r") do f
        n_cols_file = (file.separator isa AbstractVector) ? count(i->(i in file.separator), readline(path_ex)) + 1 : count(file.separator, readline(path_ex)) + 1
    end

    header = repeat([:x], n_cols_file)

    for i in 1:n_cols_file
        if haskey(file.columns, i)
            types[i] = GenVarInfo_Types[file.columns[i]]
            header[i] = (file isa GWAS) ? GenVarInfo_Symbols_out[file.columns[i]] : GenVarInfo_Symbols_exp[file.columns[i]]
        else
            types[i] = String
        end
    end
    return types, header
end


function verify_and_simplify_columns(exposure::Union{QTLStudy, GWAS})::Bool

    if exposure isa QTLStudy
        trait_each_path_nothing =  nothing in exposure.traits_for_each_path
        if trait_each_path_nothing && !(TRAIT_NAME in values(exposure.columns))
            throw(ArgumentError("exposure misses TRAIT_NAME information"))
        end
    else
        trait_each_path_nothing = false
    end


    # Treating exposure
    # keep only columns and verify all columns are satisfied       ----> cette section pourrait être optimisée? Encapsuler dans une foncion!
    new_cols_exposure::Dict{Int, GenVarInfo} = Dict()
    col_log_pval = -1
    cols_ok = Dict(CHR => false,
                   POS => false,
                   BETA => false,
                   SE => false,
                   PVAL => false,
                   A_EFFECT => false,
                   A_OTHER => false,
                   TRAIT_NAME => !trait_each_path_nothing)

    for key in keys(exposure.columns)
        for info in [CHR, POS, BETA, SE, A_EFFECT, A_OTHER]
            if exposure.columns[key] == info
                cols_ok[info] = true
                new_cols_exposure[key] = info
            end
        end

        if exposure.columns[key] == PVAL
            cols_ok[PVAL] = true
            new_cols_exposure[key] = PVAL
            col_log_pval = -1
        elseif cols_ok[PVAL] == false && exposure.columns[key] == LOG_PVAL
            col_log_pval = key
        end
    end

    if col_log_pval != -1
        new_cols_exposure[LOG_PVAL] = col_log_pval
        cols_ok[PVAL] = true
    end

    for key in keys(cols_ok)
        if cols_ok[key] != true
            throw(ArgumentError("Missing information in columns. Should contain at least : CHR, POS, A_EFFECT, A_OTHER, BETA, SE, PVAL/MINUS_LOG10_PVAL."))
        end
    end
    exposure.columns = new_cols_exposure
    return col_log_pval != -1

end


function read_qtl_files(exposure::QTLStudy, 
    col_log_pval::Bool, 
    types::Dict{Int, DataType}, 
    header::Vector{Symbol},
    window::Int,
    p_tresh::Float64,
    filtered::Bool =false,
    trsf_log_pval::Function = x -> exp10.(x))::Dataset

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

    if Base.Threads.nthreads() > length(exposure) ### Case when nore threads than files -> read each file in //
        for i in 1:lastindex(data_vect)
            file = exposure[i]
            d = filereader(file.path, delimiter = file.separator, 
                           header = header, types = types, skipto=2, 
                           makeunique=true, eolwarn=false)[:,collect(keys(file.columns))]
            
            if add_trait_name_b
                d.trait = repeat([file.trait_name], nrow(d))
            end
            
            if col_log_pval
                modify!(d, :pval => trsf_log_pval)
            end

            if !filtered
                d = @chain d  begin 
                    filter([:chr, :pos, :trait, :pval_exp], 
                        type = in_window) # dataset filtered for window and significance
                    filter(:a_effect_exp, type = x -> length(x) == 1) # remove indels
                    filter(:a_other_exp, type = x -> length(x) == 1)  # remove indels
                end
            end
            
            modify!(d, [:a_effect_exp, :a_other_exp] => x -> lowercase.(x))

            data_vect[i] = d
        end
    else                                        ### Case when more files than threads -> read multiple files in //
        @threads for i in 1:lastindex(data_vect)
            file = exposure[i]
            d = filereader(file.path, delimiter = file.separator, 
                           header = header, types = types, skipto=2, 
                           makeunique=true, eolwarn=false, 
                           threads = false)[:,collect(keys(file.columns))]
            
            if add_trait_name_b
                d.trait = repeat([file.trait_name], nrow(d))
            end
            
            if col_log_pval
                modify!(d, :pval => trsf_log_pval, threads = false)
            end

            if !filtered
                d = @chain d  begin 
                    filter([:chr, :pos, :trait, :pval_exp], 
                        type = in_window, threads = false) # dataset filtered for window and significance
                    filter(:a_effect_exp, type = x -> length(x) == 1,  # remove indels
                        threads = false)
                    filter(:a_other_exp, type = x -> length(x) == 1,  # remove indels
                        threads = false)
                end
            end

            modify!(d, [:a_effect_exp, :a_other_exp] => x -> lowercase.(x), threads = false)
            
            data_vect[i] = d
        end
    end

    data_filtered = Folds.reduce(vcat, data_vect, init = Dataset())

    return data_filtered
end


    #########################
    #      MrStudyCis       #
    #########################


"""
Perform a Mendelian Randomization study with exposure QTL and outcome GWAS
"""
function mrStudyCis(exposure::QTLStudy, 
    outcome::GWAS, 
    bedbimfam_dirnames::AbstractArray{String};
    approach::String="naive", 
    p_tresh::Float64 = 1e-3, 
    window::Int = 500000, 
    r2_tresh::Float64 = 0.1,
    exposure_filtered = false,
    mr_methods::AbstractVector{Function} = [mr_egger, mr_ivw],
    α::Float64 = 0.05,
    trsf_log_pval_exp::Function = x -> exp10.(x),
    trsf_log_pval_out::Function = x -> exp10.(x)
    )::Union{Dataset, GroupBy}
    
    # plink_files_load_tsk = @async global GenotypesArr = [SnpData(SnpArrays.datadir(file)) for file in bedbimfam_dirnames]

    #qtl_data
    col_log_pval = verify_and_simplify_columns(exposure)
    types, header = make_types_and_headers(exposure)
    qtl_d = read_qtl_files(exposure, col_log_pval, types, header, window, p_tresh, exposure_filtered, trsf_log_pval_exp)
    
    #gwas_data
    col_log_pval = verify_and_simplify_columns(outcome)
    types, header = make_types_and_headers(outcome)
    gwas_d = filereader(outcome.path, 
                        delimiter = outcome.separator, 
                        header = header, types = types, skipto=2, 
                        makeunique=true, eolwarn=false)[:,collect(keys(outcome.columns))]
    if col_log_pval
        modify!(gwas_d, :pval => trsf_log_pval_out)
    end
    modify!(gwas_d, [:a_effect_out, :a_other_out] => x -> lowercase.(x))

    biallelic(s::SubArray) = (s[1]==s[2] && s[3] == s[4]) || (s[1] == s[4] && s[2] == s[3])
    joined_d = @chain qtl_d begin
        innerjoin(gwas_d, on = [:chr, :pos], makeunique = false)
        filter([:a_effect_exp, :a_effect_out, :a_other_exp, :a_other_out], type = biallelic) # filter for "obvious" non biallelic variants
        filter(:, by = !ismissing)
        groupby(:trait, stable = false)
    end

    # wait for loaded plink files and format them
    # wait(plink_files_load_tsk)
    GenotypesArr = [SnpData(SnpArrays.datadir(file)) for file in bedbimfam_dirnames]
    @threads for i in 1:lastindex(GenotypesArr)
        formatSnpData!(GenotypesArr[i])
    end

    one_file_per_chr_plink = length(bedbimfam_dirnames) > 1

    #### for d in eachgroup(joined_d) -> Plink + MR (implement in NaiveCis)
    if approach == "naive"
        return NaiveCis(joined_d, r2_tresh, GenotypesArr, one_file_per_chr_plink, mr_methods, α)
    elseif approach == "test"
        return joined_d
    elseif approach == "strict"
        return Dataset()
    elseif approach == "2ndChance"
        return Dataset()
    else
        throw(ArgumentError("approach should be either \"naive\", \"test\", \"strict\" or \"2ndChance\""))
    end
end