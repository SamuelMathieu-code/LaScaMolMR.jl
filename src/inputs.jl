
##################################################
#                   ValTypeGen                   #
##################################################

module ValTypeGen
@enum GenVarInfo begin
    chr
    pos
    rsid
    chr_pos
    chr_colong_pos
    exposition_name
    exposition_isoform
    exposition_name_isoform
    β
    se
    maf
    odd_ratio
    # Add others ???
    other
end
end #ValTypeGen


###################################################
#                      Inputs                     #
###################################################


module Inputs

using ..ValTypeGen

struct GWAS
    path::String
    columns::Vector{ValTypeGen.GenVarInfo}
    separator::Union{String, Char}
end

struct QtlStudy
    paths::Vector{String}
    columns::Vector{ValTypeGen.GenVarInfo}
    separator::Union{Char, String}
end

"""
Return Regex from Vector of Strings and ValTypeGen.GenVarInfo
"""
function make_regex_path(path1::String, path2::Vector{Any})::Regex
    str = "";
    for i in path2
        if (i isa String); str*=i; else; str*="(.)"; end;     # pourra etre augmenté pour prendre en compte les types d'infos et les conditions a respecter.
    end
    str = joinpath(path1, str);
    return Regex(str)
end

"""
Return all existing files corresponding to the regex defined
"""
function find_all_corresp(path1::String, path2::Vector{Any})::Vector{String}
    #...
    return []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names.
"""
function find_all_corresp(path1::String, path2::Vector{Any}, 
                          genes::Vector{String})::Vector{String}
    #...
    return []
end


"""
Find all existing files corresponding to the regex defined and constrait on gene names and chromosomes.
"""
function find_all_corresp(path1::String, path2::Vector{Any}, 
                          genes::Vector{String},
                          chr::Vector{String})::Vector{String}
    #...
    return []
end

QtlStudy(path1::String, 
         path2::Vector{Any}, 
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path1, path2), columns, separator)

QtlStudy(path1::String, 
         path2::Vector{Any},
         genes::Vector{String},
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path1, path2, genes), columns, separator)

QtlStudy(path1::String, 
         path2::Vector{Any},
         genes::Vector{String},
         chr::Vector{String},
         columns::Vector{ValTypeGen.GenVarInfo}, 
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path1, path2, genes, chr), columns, separator)

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

end #Inputs