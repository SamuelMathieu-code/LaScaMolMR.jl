
using Glob

##################################################
#                   GenVarInfo                   #
##################################################

macro exportinstances(enum)
    eval = GlobalRef(Core, :eval)
    return :($eval($__module__, Expr(:export, map(Symbol, instances($enum))...)))
end


@enum GenVarInfo begin
    TRAIT_NAME
    # TRAIT_NAME_ISO
    # TRAIT_ISO       # à voir .....
    CHR
    POS
    RSID
    CHR_POS
    CHR_COLON_POS
    CHR_COLON_POS_ALLELES
    A_EFFECT
    A_OTHER
    A1_COLON_A2
    BETA
    SE
    MAF
    ODD_RATIO
    CI_LOW
    CI_HIGH
    PVAL
    LOG10_PVAL
    MINUS_LOG10_PVAL
    # add others ?
    OTHER_INFO
end

@exportinstances GenVarInfo

QtlPathPattern = Union{GenVarInfo, String}

promote_rule(::Type{S}, ::Type{T}) where {S <: QtlPathPattern, T <: QtlPathPattern} = QtlPathPattern


###################################################
#                      Inputs                     #
###################################################

struct GWAS
    path::String
    columns::Union{Dict{Int64, Union{GenVarInfo, String}}, Dict{Int64, GenVarInfo}}
    separator::Char
    trait_name::Union{Nothing, String}  # When searching pot ivs --> check trait_name id ok with ouput of GWAS.acc_func <=> GWAS.trait_name !== nothing
end

GWAS(path::String,
     columns::Union{Dict{Int64, Union{GenVarInfo, String}}, Dict{Int64, GenVarInfo}},
     separator::Char) = GWAS(path, columns, separator, nothing)


struct QtlStudy      # Change to see it directly as a collection of gwas??
    path_v::Vector{String}
    trait_v::Union{Nothing, Vector{String}}
    chr_v::Vector{Int8}
    tss_v::Vector{Int32}
    columns::Dict{Int64, Union{GenVarInfo, String}}
    separator::Char
end


# Ne fonctionne pas : à recoder plus tard
#                    |
#                    |
#                    V
# """
# Find all existing files corresponding to specified path pattern and needed traits and chromosome specifications. 
# To have all chromosomes for one trait : duplicate trait_name for all chromosomes
# Returns a vector of corresponding paths
# """
# function find_all_corresp(path::Vector{QtlPathPattern},                                         # Note pour le futur : faire la duplication auto si chr ==0 ? ça pourrait etre cool.
#                           trait_v::Vector{String},
#                           chr_v::Vector{Int64})::Vector{String}
#     if (length(chr_v) !== length(trait_v))
#          throw(ErrorException("trait_v, chr_v should be of same length."))
#     end
#     l = length(trait_v)
#     path_temp::String = ""
#     path_v::Vector{String} = []
#     for i = 1:1:l
#         for ob in path
#             if ob == CHR
#                 path_temp = path_temp*string(chr_v[i])
#             elseif ob == TRAIT_NAME
#                 path_temp = path_temp*trait_v[i]
#             elseif ob isa String
#                 path_temp = path_temp*ob
#             else
#                 path_temp = path_temp*"*"
#             end
#         end
#         found_paths = glob(path_temp)

#         append!(path_v, found_paths)
#     end
#     return path_v
# end


# """
# struct defining format of QTL sumstats data and location
# """
# QtlStudy(path_pattern::Vector{QtlPathPattern},
#          trait_v::Vector{String},
#          chr_v::Vector{Int64},
#          tss_v::Vector{Int64},
#          columns::Dict{Int64, Union{GenVarInfo, String}},
#          separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path_pattern, trait_v, chr_v), trait_v, chr_v, tss_v, columns, separator)


# Iteration overload for QtlStudy
function Base.iterate(iter::QtlStudy)
    element = GWAS(iter.path_v[1], iter.columns, iter.separator, iter.trait_v[1])
    return (element, 1)
end


function Base.iterate(iter::QtlStudy, state)
    count = state + 1
    if count > length(iter.path_v)
        return nothing
    end
    element = GWAS(iter.path_v[count], iter.columns, iter.separator)
    return (element, count)
end

