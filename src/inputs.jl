
using Glob

##################################################
#                     Utils                      #
##################################################

inv_logit(x) = exp(x)/(1+exp(x));

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
    BETA
    SE
    MAF
    EFFECT_ALLELE # à rajouter dans l'extraction d'infos du GWAS
    OTHER_ALLELE
    ODD_RATIO
    CI_LOW
    CI_HIGH
    PVAL
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
    columns::Dict{GenVarInfo, Int64}
    separator::Union{String, Char}
    trait_name::Union{Nothing, String}  # When searching pot ivs --> check trait_name id ok with ouput of GWAS.acc_func <=> GWAS.trait_name !== nothing
    chr::Union{Nothing, Int64}
    tss::Union{Nothing, Int64}
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
        get_effect = line -> [parse(Float64, line[columns[BETA]]), parse(Float64, line[columns[SE]]), parse(Float64, line[columns[PVAL]])]
    elseif haskey(columns, ODD_RATIO) && haskey(columns, CI_LOW) && haskey(columns, CI_HIGH)
        get_effect = line -> [inv_logit(parse(Float64, line[columns[ODD_RATIO]])), 
                              (inv_logit(parse(Float64, line[columns[CI_HIGH]]))-inv_logit(parse(Float64, line[columns[CI_LOW]])))/3.91992, 
                              parse(Float64, line[columns[PVAL]])]
    else
        throw(ErrorException("Not either (BETA, SE) or (OR, CI) in columns"))
    end

    # Treat cases of variant info formating
    if haskey(columns, CHR) && haskey(columns, POS)
        get_var = line -> [parse(Int64, line[columns[CHR]]), parse(Int64, line[columns[POS]])]
    elseif haskey(columns, CHR_COLON_POS)
        get_var = line -> [parse(Int64, x) for x in split(line[columns[CHR_COLON_POS]], ':')]
    elseif haskey(columns, CHR_POS)
        get_var = line -> [parse(Int64, x) for x in split(line[columns[CHR_POS]], '_')]
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
     separator::Union{String, Char}) = GWAS(path, columns, separator, nothing, nothing, nothing, make_func(columns, separator))

GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{Char, String},
     trait_name::String) = GWAS(path, columns, separator, trait_name, nothing, nothing, make_func(columns, separator))

GWAS(path::String,
     columns::Dict{GenVarInfo, Int64},
     separator::Union{String, Char},
     trait_name::String,
     chr::Int64,
     tss::Int64) = GWAS(path, columns, separator, trait_name, chr, tss, make_func(columns, separator))


struct QtlStudy      # Change to see it directly as a collection of gwas??
    path_v::Vector{String}
    trait_v::Vector{String}
    chr_v::Vector{Int64}
    tss_v::Vector{Int64}
    columns::Dict{GenVarInfo, Int64}
    separator::Union{Char, String}
    acc_func::Function
end


"""
Find all existing files corresponding to specified path pattern and needed traits and chromosome specifications. 
To have all chromosomes for one trait : duplicate trait_name for all chromosomes
Returns a vector of corresponding paths
"""
function find_all_corresp(path::Vector{QtlPathPattern},                                         # Note pour le futur : faire la duplication auto si chr ==0 ? ça pourrait etre cool.
                          trait_v::Vector{String},
                          chr_v::Vector{Int64},
                          tss_v::Vector{Int64})::Vector{String}
    if (length(tss_v) !== length(chr_v) || length(tss_v) !== length(trait_v) || length(chr_v) !== length(trait_v))
         throw(ErrorException("trait_v, chr_v and tss_v should be of same length."))
    end
    l = length(trait_v)
    path_temp::String = ""
    path_v::Vector{String} = []
    for i = 1:1:l 
        for ob in path
            if ob == CHR
                path_temp = path_temp*string(chr_v[i])
            elseif ob == TRAIT_NAME
                path_temp = path_temp*trait_v[i]
            elseif ob isa String
                path_temp = path_temp*ob
            else
                path_temp = path_temp*"*"
            end
        end
        found_paths = glob(path_temp)
        trait_temp, chr_temp = trait_v[i], chr_v[i]
        ll = length(found_paths)
        if ll > 1
            trait_temp, chr_temp = trait_v[i], chr_v[i]
            @warn "More than one corresponding file was found for trait $trait_temp and chromosome $chr_temp. 
                   (trait, chr, tss) may appear duplicated in QtlStudy"
            for j = 1:1:(ll-1)
                insert!(trait_v, i+1, trait_v[i])
                insert!(chr_v, i+1, chr_v[i])
                insert!(tss_v, i+1, tss_v[i])
            end
        elseif ll == 0
            @warn "No file corresponding to trait $trait_temp and chromosome $chr_temp. Instance will be skipped."
            continue
        end
        append!(path_v, found_paths)
    end
    return path_v
end


QtlStudy(path_pattern::Vector{QtlPathPattern},
         trait_v::Vector{String},
         chr_v::Vector{Int64},
         tss_v::Vector{Int64},
         columns::Dict{GenVarInfo, Int64},
         separator::Union{String, Char}) =  QtlStudy(find_all_corresp(path_pattern, trait_v, chr_v, tss_v), trait_v, chr_v, tss_v, columns, separator, make_func(columns, separator))


QtlStudy(path_v::Vector{String},
         trait_v::Vector{String},
         chr_v::Vector{Int64},
         tss_v::Vector{Int64},
         columns::Dict{GenVarInfo, Int64},
         separator::Union{String, Char}) =  QtlStudy(path_v, trait_v, chr_v, tss_v, columns, separator, make_func(columns, separator))


# Iteration overload for QtlStudy
function Base.iterate(iter::QtlStudy)
    element = GWAS(iter.path_v[1], ietr.columns, iter.separator, iter.trait_v[1], iter.chr_v[1], iter.tss_v[1], iter.acc_func)
    return (element, 1)
end


function Base.iterate(iter::QtlStudy, state)    #### Comment faire l'itération pour avoir un GWAS par exposition et pas necéssairement par fichier??
    count = state + 1
    element = GWAS(iter.path_v[count], ietr.columns, iter.separator, iter.trait_v[count], iter.chr_v[count], iter.tss_v[count], iter.acc_func)
    return (element, count)
end

