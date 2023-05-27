using DLMReader
using InMemoryDatasets
import Base.Threads.@threads
using Folds
using Chain

    ########################
    #       Constants      #
    ########################

GenVarInfo_Types = Dict(TRAIT_NAME => String,
                        CHR => Int8,
                        POS => Int32,
                        A_EFFECT => String,
                        A_OTHER => String,
                        BETA => Float64,
                        SE => Float64,
                        PVAL => Float64,
                        MINUS_LOG10_PVAL => Float64)

GenVarInfo_Symbols_exp = Dict(TRAIT_NAME => :trait,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_exp,
                          A_OTHER => :a_other_exp,
                          BETA => :β_exp,
                          SE => :se_exp,
                          PVAL => :pval_exp,
                          MINUS_LOG10_PVAL => :pval_exp)

GenVarInfo_Symbols_out = Dict(TRAIT_NAME => :trait,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect_out,
                          A_OTHER => :a_other_out,
                          BETA => :β_out,
                          SE => :se_out,
                          PVAL => :pval_out,
                          MINUS_LOG10_PVAL => :pval_out)


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


function verify_and_simplify_columns(exposure::Union{QTLStudy, GWAS})::Int

    trait_each_path_nothing =  nothing in exposure.traits_for_each_path
    if trait_each_path_nothing && !(TRAIT_NAME in values(exposure.columns))
        throw(ArgumentError("exposure misses TRAIT_NAME information"))
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
        elseif cols_ok[PVAL] == false && exposure.columns[key] == MINUS_LOG10_PVAL
            col_log_pval = key
        end
    end

    if col_log_pval != -1
        new_cols_exposure[MINUS_LOG10_PVAL] = col_log_pval
        cols_ok[PVAL] = true
    end

    for key in keys(cols_ok)
        if cols_ok[key] != true
            throw(ArgumentError("Missing information in columns. Should contain at least : CHR, POS, A_EFFECT, A_OTHER, BETA, SE, PVAL/MINUS_LOG10_PVAL."))
        end
    end
    exposure.columns = new_cols_exposure
    return col_log_pval

end


function read_files(exposure::QTLStudy, 
    col_log_pval::Int, 
    types::Dict{Int, DataType}, 
    header::Vector{Symbol},
    window::Int,
    p_tresh::Float64,
    filtered::Bool =false)::Dataset

    # Dicionary of trait -- tss
    ref_dict = Dict(zip(exposure.trait_v, zip(exposure.chr_v, exposure.tss_v)))

    # boolan tells if the variant is significant causal on exposure and if in window arround good tss
    in_window(s::SubArray) = (s[1] == ref_dict[s[3]][1] && abs(s[2]-ref_dict[s[3]][2])≤window && s[4]<p_tresh)
    
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
            
            if col_log_pval != -1
                modify!(d, :pval => to_p)
            end

            if !filtered
                filter!(d, [:chr, :pos, :trait, :pval], 
                        type = in_window) # dataset filtered for window and significance
            end
            
            data_vect[i] = d
        end
    else                                        ### Case when more files than threads -> read multiple files in //
        @threads for i in 1:lastindex(data_vect)
            file = exposure[i]
            d = filereader(file.path, delimiter = file.separator, 
                           header = header, types = types, skipto=2, 
                           makeunique=true, eolwarn=false, 
                           threads = flase)[:,collect(keys(file.columns))]
            
            if add_trait_name_b
                d.trait = repeat([file.trait_name], nrow(d))
            end
            
            if col_log_pval != -1
                modify!(d, :pval => to_p, threads = false)
            end

            if !filtered
                filter!(d, [:chr, :pos, :trait, :pval], 
                        type = in_window, threads = false) # dataset filtered for window and significance
            end
            
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
    approach::String="naive", 
    p_tresh::Float64 = 5e-3, 
    window::Int = 500000, 
    r2_tresh::Float64 = 0.1,
    exposure_filtered = false)::AbstractDataset
    
    #qtl_data
    col_log_pval = verify_and_simplify_columns(exposure)
    types, header = make_types_and_headers(exposure)
    qtl_d = read_files(exposure, col_log_pval, types, header, window, p_tresh, exposure_filtered)
    
    #gwas_data
    col_log_pval = verify_and_simplify_columns(outcome)
    types, header = make_types_and_headers(outcome)
    gwas_d = filereader(outcome.path, 
                        delimiter = outcome.separator, 
                        header = header, types = types, skipto=2, 
                        makeunique=true, eolwarn=false)[:,collect(keys(outcome.columns))]

    #joined data
    joined_d = innerjoin(gwas_d, qtl_d, on = [:chr, :pos], makeunique = true)
    groupby!(joined_d, :prots, stable = false)
    
    #### for d in eachgroup(joined_d) -> Plink + MR (implement in NaiveCis)
    
    return NaiveCis(joined_d, r2_tresh)

end