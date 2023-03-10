
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
    β
    SE
    MAF
    ODD_RATIO
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
    chr::Union{nothing, Int, String} # Int or string?
    trait_iso::Union{nothing, String}
end

GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{String, Char}) = GWAS(path, columns, separator, nothing, nothing, nothing)

GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{Char, String}
     trait_name::String) = GWAS(path, columns, separator, trait_name, nothing, nothing)

GWAS(path::String,
     columns::Vector{GenVarInfo},
     separator::Union{String, Char},
     trait_name::String,
     trait_iso::String) = GWAS(path, columns, separator, trait_name, nothing, trait_iso)

GWAS(path::String,
     columns::Vector{String},
     separator::Union{String, Char},
     trait_name::String,
     chr::Int) = GWAS(path, columns, separator, trait_name, chr, nothing)

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

# Iteration overload
function Base.iterate(iter::QtlStudy)
    element = GWAS(iter.paths[1], iter.columns, iter.separator);
    return (element, 1)
end

function Base.iterate(iter::QtlStudy, state)
    count = state+1;
    if count>length(iter.paths)
        return nothing
    end
    element = GWAS(iter.paths[count], iter.columns, iter.separator);
    return (element, count)
end

