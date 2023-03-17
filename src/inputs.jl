
##################################################
#                     Utils                      #
##################################################

inv_logit(x) = exp(x)/(1+exp(x));

##################################################
#                   GenVarInfo                   #
##################################################


@enum GenVarInfo begin
    TRAIT_NAME
    # TRAIT_NAME_ISO
    # TRAIT_ISO       # à voir .....
    CHR
    POS
    RSID
    CHR_POS
    CHR_COLON_POS
    BETA
    SE
    MAF
    ODD_RATIO
    CI_LOW
    CI_HIGH
    PVAL
    # add others ?
    OTHER_INFO
end

QtlPathPattern = Union{GenVarInfo, String}

promote_rule(::Type{S}, ::Type{T}) where {S <: QtlPathPattern, T <: QtlPathPattern} = QtlPathPattern


###################################################
#                      Inputs                     #
###################################################

struct GWAS
    path::String
    columns::Dict{GenVarInfo, Int64}
    separator::Union{String, Char}
    trait_name::Union{Nothing, String}  # When searching pot ivs --> check trait_name id ok with ouput of GWAS.acc_func <=> GWAS.trait_name !== nothing
    chr::Union{Nothing, String}
    acc_func::Function
end


"""
Return function that extracts (trait_name, [chr, pos], [beta, se, pval]) from line of file (string)  #### NOTE : CHANGE COLUMNS from VECTOR to DICT(GenVarInfo -> Int)
    Raises ErrorException if some information is missiong in `columns` argument.
"""
function make_func(columns::Dict{GenVarInfo, Int64}, 
                   sep::Union{String, Char})::Function
    
    if !(haskey(columns, PVAL)); throw(ErrorException("No pvalue info in columns")); end
    # treat cases of effect info
    if haskey(columns, BETA) && haskey(columns, SE)
        get_effect = line -> [line[columns[BETA]], line[columns[SE]], line[columns[PVAL]]]
    elseif haskey(columns, OR) && haskey(columns, CI_LOW) && haskey(columns, CI_HIGH)
        get_effect = line -> [inv_logit(line[columns[ODD_RATIO]]), (inv_logit(line[columns[CI_HIGH]])-inv_logit(line[columns[CI_LOW]]))/3.91992]    #### Can I optimize this better??
    else
        throw(ErrorException("Not either (BETA, SE) or (OR, CI) in columns"))
    end

    # Treat cases of variant info formating
    if haskey(columns, CHR) && haskey(columns, POS)
        get_var = line -> [line[columns[CHR]], line[columns[POS]]]
    elseif haskey(columns, CHR_COLON_POS)
        get_var = line -> split(line[columns[CHR_COLON_POS]], ':')
    elseif haskey(columns, CHR_POS)
        get_var = line -> split(line[columns[CHR_POS]], '_')
    else
        throw(ErrorException("Missing CHR and POS information"))
    end

    if haskey(columns, TRAIT_NAME)                                       #### See how we can treat the case where trait_name is composed of base + iso in different columns
        get_trait = line -> line[columns[TRAIT_NAME]]
    else
        get_trait = line -> nothing
    end

    function f(sline::String)
        line = split(sline, sep)
        return get_trait(line), get_var(line), get_effect(line)
    end

    return f
end


GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{String, Char}) = GWAS(path, columns, separator, nothing, nothing, make_func(columns, separator))

GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{Char, String},
     trait_name::String) = GWAS(path, columns, separator, trait_name, nothing, make_func(columns, separator))

GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{String, Char},
     trait_name::String,
     chr::Int) = GWAS(path, columns, separator, trait_name, chr, make_func(columns, separator))

GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{String, Char},
     chr::Int) = GWAS(path, columns, separator, nothing, chr, make_func(columns, separator))


struct QtlStudy      # Change to see it directly as a collection of gwas??
    paths::Vector{String}
    columns::Dict{GenVarInfo, Int64}
    separator::Union{Char, String}
    acc_func::Function
end

# """
# Return Regex from Vector of Strings and ValTypeGen.GenVarInfo
# """
# function make_regex_path(path::Vector{Any})::Regex
#     str = "";
#     for i in path
#         if (i isa String); str*=i; else; str*="(.)"; end;     # pourra etre augmenté pour prendre en compte les types d'infos et les conditions a respecter.
#     end
#     return Regex(str)
# end

"""
Return all existing files corresponding to the regex defined
"""
function find_all_corresp(path::Vector{QtlPathPattern})::Tuple{Vector{String}, Vector{String}}
    #...
    return [], []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names.
"""
function find_all_corresp(path::Vector{QtlPathPattern}, 
                          genes::Vector{String})::Tuple{Vector{String}, Vector{String}}
    #...
    return [], []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names and chromosomes specified (particularly in a cis iv selection process).
Returns a tuple of arrays corresponding to traits included and corresponding paths.
"""
function find_all_corresp(path::Vector{QtlPathPattern}, 
                          traits::Vector{String},
                          chr::Vector{String})::Tuple{Vector{String}, Vector{String}}
    #...
    return [], []
end


QtlStudy((paths, traits), columns, separator, acc_func) = QtlStudy(paths, columns, separator, acc_func)

QtlStudy(path::Vector{QtlPathPattern}, 
         columns::Dict{GenVarInfo, Int64}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path), columns, separator, make_func(columns, separator))

QtlStudy(path::Vector{QtlPathPattern},
         traits::Vector{String},
         columns::Dict{GenVarInfo, Int64},
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path, traits), columns, separator, make_func)

QtlStudy(path::Vector{QtlPathPattern},
         traits::Vector{String},
         chr::Vector{String},
         columns::Dict{GenVarInfo, Int64}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path, traits, chr), columns, separator)


# Iteration overload for QtlStudy
function Base.iterate(iter::QtlStudy)
    element = GWAS(iter.paths[1], iter.columns, iter.separator);
    return (element, 1)
end


function Base.iterate(iter::QtlStudy, state)    #### Comment faire l'itération pour avoir un GWAS par exposition et pas necéssairement par fichier??
    count = state+1;
    if count>length(iter.paths)
        return nothing
    end
    element = GWAS(iter.paths[count], iter.columns, iter.separator);
    return (element, count)
end

