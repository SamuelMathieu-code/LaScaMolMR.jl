
##################################################
#                   GenVarInfo                   #
##################################################


@enum GenVarInfo begin
    TRAIT_NAME
    TRAIT_NAME_ISO
    TRAIT_ISO
    CHR
    POS
    RSID
    CHR_POS
    CHR_COLONG_POS
    BETA
    SE
    MAF
    ODD_RATIO
    CI
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
    columns::Vector{GenVarInfo}
    separator::Union{String, Char}
    trait_name::Union{nothing, String}
    chr::Union{nothing, Int8} # Int or string?
    trait_iso::Union{nothing, String}
    acc_func::Function # -----------------------> acces infos from line of file
end


"""
Return function that extracts ([chr, pos], [beta, se, pval])
"""
function make_func(path::String, 
                   columns::Vector{GenVarInfo}, 
                   sep::Union{String, Char}, 
                   trait_name::Union{String, nothing} = nothing, 
                   chr::Union{Int8, nothing}=nothing, 
                   trait_iso::Union{nothing, String}=nothing)::Function
    
    # return a function that threats a line (String) and returns a tuple of vectors : ([chr, pos], [beta, se, pval]) --> raise err if not enough info in columns
    # Traiter cas par cas pour ne pas avoir à le faire à chaque MR. Faire une version pour QtlStudy pour accélérer?
    # cas qui foirent -> 
    #                   si un de chr, trait_iso, trait_name != nothing et pas dans columns -> foire.
    #                   si ni (β, SE) ou (OR, CI) dans columns -> foire
    #                   si ni (CHR, POS), CHR_POS, CHR_COLONG_POS dans columns -> foire 
    
    if chr!==nothing && !(CHR in columns) && !(CHR_COLONG_POS in columns) && !(CHR_POS in columns); throw(ErrorException("no argument chr in columns of GWAS")); end
    if trait_name!==nothing && !(TRAIT_NAME in columns) && !(TRAIT_NAME_ISO in columns); throw(ErrorException("no argument TRAIT_NAME in columns of GWAS")); end
    if trait_iso!==nothing && !(TRAIT_ISO in columns) && !(TRAIT_NAME_ISO in columns); throw(ErrorException("no argument trait_iso in columns of GWAS")); end
    if !(CHR_COLONG_POS in cloumns) || !(CHR_POS in columns) || !((CHR in columns) && (POS in columns)); throw(ErrorException("Not enough information about variants in columns")); end
    if !(PVAL in columns); throw(ErrorException("No pvalue info in columns")); end
    
    if BETA in columns && SE in columns
        
        #...
    elseif OR in columns && CI in columns
        
        #...
    else
        throw(ErrorException("Not either (BETA, SE) or (OR, CI) in columns"))
    end

    return (line::String -> ([0, 0], [0, 0, 0]))
end


GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{String, Char}) = GWAS(path, columns, separator, nothing, nothing, nothing, make_func(path, columns, separator))

GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{Char, String},
     trait_name::String) = GWAS(path, columns, separator, trait_name, nothing, nothing, make_func(path, columns, separator, trait_name))

GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{String, Char},
     trait_name::String,
     trait_iso::String) = GWAS(path, columns, separator, trait_name, nothing, trait_iso, make_func(path, columns, trait_name, nothing, trait_iso))

GWAS(path::String,
     columns::Vector{String},
     separator::Union{String, Char},
     trait_name::String,
     chr::Int) = GWAS(path, columns, separator, trait_name, chr, nothing, make_func(path, columns, trait_name, chr))

GWAS(path::String,
     columns::Vector{String},
     separator::Union{String, Char},
     trait_name::String,
     chr::Int,
     trait_iso::String) = GWAS(path, columns, separator, trait_name, chr, trait_iso, make_func(path, columns, trait_name, chr, trait_iso))


struct QtlStudy
    paths::Vector{String}
    columns::Vector{GenVarInfo}
    separator::Union{Char, String}
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
function find_all_corresp(path::Vector{QtlPathPattern})::Vector{String}
    #...
    return []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names.
"""
function find_all_corresp(path2::Vector{QtlPathPattern}, 
                          genes::Vector{String})::Vector{String}
    #...
    return []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names and chromosomes specified (particularly in a cis iv selection process).
"""
function find_all_corresp(path2::Vector{QtlPathPattern}, 
                          genes::Vector{String},
                          chr::Vector{String})::Vector{String}
    #...
    return []
end

QtlStudy(path::Vector{QtlPathPattern}, 
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path), columns, separator)

QtlStudy(path::Vector{QtlPathPattern},
         genes::Vector{String},
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path, genes), columns, separator)

QtlStudy(path::Vector{QtlPathPattern},
         genes::Vector{String},
         chr::Vector{String},
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path, genes, chr), columns, separator)


# Iteration overload for QtlStudy
function Base.iterate(iter::QtlStudy)
    element = GWAS(iter.paths[1], iter.columns, iter.separator);
    return (element, 1)
end


function Base.iterate(iter::QtlStudy, state)
    count = state+1;
    if count>length(iter.paths)
        return nothing
    end
    element = GWAS(iter.paths[count], iter.columns, iter.separator); ##### TO BE MODIF WITH NEW ADDINGS
    return (element, count)
end

