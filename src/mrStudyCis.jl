using DLMReader
using InMemoryDatasets

GenVarInfo_Types = Dict(TRAIT_NAME => String,
                        CHR => Int8,
                        POS => Int32,
                        A_EFFECT => String,
                        A_OTHER => String,
                        BETA => Float64,
                        SE => Float64,
                        PVAL => Float64,
                        MINUS_LOG10_PVAL => Float64)
GenVarInfo_Symbols = Dict(TRAIT_NAME => :trait,
                          CHR => :chr,
                          POS => :pos,
                          A_EFFECT => :a_effect,
                          A_OTHER => :a_other,
                          BETA => :β,
                          SE => :se,
                          PVAL => :pval,
                          MINUS_LOG10_PVAL => :pval)

function make_types_and_headers(file::QTLStudy)::Tuple{Dict{Int, DataType}, Vector{Symbols}}
    types::Dict{Int, DataType} = Dict()
    open(file.path_v[1], "r") do f
        n_cols_file = (file.separator isa AbstractVector) ? count(i->(i in file.separator), readline(file.path_v[1])) : count(file.separator, readline(file.path))
    end
    header = repeat([:x], n_cols_file)
    for i in 1:n_cols_file
        if haskey(file.columns, i)
            types[i] = GenVarInfo_Types[file.columns[i]]
            header[i] = GenVarInfo_Symbols[file.columns[i]]
        else
            types[i] = String
        end
    end
    return header, types
end

function verify_and_simplify_columns(exposure::QTLStudy)
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
        elseif cls_ok[PVAL] == false && exposure.columns[key] == MINUS_LOG10_PVAL
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
    exposure.columns = new_cols_exposure;

end

###############################
#         MrStudyCis          #
###############################

"""
Perform a Mendelian Randomization study with exposure QTL and outcome GWAS
"""
function mrStudyCis(exposure::QTLStudy, 
    outcome::GWAS, 
    approach::String="naive", 
    p_thresh::Float = 5e-3, 
    window::Int = 500000, 
    r2_tresh::Float = 0.1)::AbstractDataset
    
    verify_and_simplify_columns(exposure)
    types, header = make_types_and_headers(exposure)
    
    #Dictionary containing tss
    ref_dict = Dict(zip(exposure.trait_v, zip(exposure.chr_v, exposure.tss_v)))

    # boolan tells if the variant is significant causal on exposure and if in window arround good tss
    in_window(s::SubArray) = (s[1] == ref_dict[s[3]][1] && abs(s[2]-ref_dict[s[3]][2])≤window && s[4]<p_tresh)
    
    #first iter
    file = exposure[1]
    qtl_d = filereader(file.path, delimiter = file.separator, header = header, types = types, skipto=2, makeunique=true, eolwarn=false)[:,keys(file.columns)]
    if !(TRAIT_NAME in values(file.columns))
        qtl_d.trait = repeat([file.trait], nrow(qtl_d))
    end
    if col_log_pval != -1
        to_p(x) = exp10.(-x)
        modify!(qtl_d, :pval => to_p)
    end
    filter!(qtl_d, [:chr, :pos, :trait, :pval], type = in_window)
    
    for file in exposure[2:nend]
        d = filereader(file.path, delimiter = file.separator, header = header, types = types, skipto=2, makeunique=true, eolwarn=false)[:,keys(file.columns)]
        if !(TRAIT_NAME in values(file.columns))
            d.trait = repeat([file.trait], nrow(d))
        end
        if col_log_pval != -1
            to_p(x) = exp10.(-x)
            modify!(d, :pval => to_p)
        end
        filter!(d, [:chr, :pos, :trait, :pval], type = in_window) # dataset filtered for window and significance
        qtl_d = [qtl_d; d]
    end

    # Check for biallelic? pour l'instant : faire confiance à PLINK

    return qtl_d


end