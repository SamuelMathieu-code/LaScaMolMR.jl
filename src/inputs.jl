
using Glob
using InMemoryDatasets

##################################################
#                   GenVarInfo                   #
##################################################

macro exportinstances(enum)
    eval = GlobalRef(Core, :eval)
    return :($eval($__module__, Expr(:export, map(Symbol, instances($enum))...)))
end


@enum GenVarInfo begin
    TRAIT_NAME
    CHR
    POS
    A_EFFECT
    A_OTHER
    BETA
    SE
    PVAL
end

@exportinstances GenVarInfo

QtlPathPattern = Union{GenVarInfo, String}

promote_rule(::Type{S}, ::Type{T}) where {S <: QtlPathPattern, T <: QtlPathPattern} = QtlPathPattern


###################################################
#                      Inputs                     #
###################################################

mutable struct GWAS{T <: Integer}
    path::String
    columns::Union{Dict{T, GenVarInfo}}
    separator::Union{Char, Vector{Char}}
    trait_name::Union{Nothing, String} 
end

GWAS(path,
     columns,
     separator) = GWAS(path, columns, separator, nothing)


mutable struct QTLStudy
    path_v::Vector{S1} where S1 <: Union{Missing, AbstractString}
    traits_for_each_path::Vector{SN} where SN <: Union{String, Missing, Nothing}
    trait_v::Vector{S2} where S2 <: Union{Missing, AbstractString}
    chr_v::Vector{I1} where I1 <: Union{Missing, Integer}
    tss_v::Vector{I2} where I2 <: Union{Missing, Integer}
    columns::Union{Dict{I3, GenVarInfo}} where I3 <: Integer
    separator::Union{Char, Vector{Char}}
end

QTLStudy(path_v::Vector,
    trait_v,
    chr_v,
    tss_v,
    columns,
    separator) = QTLStudy(path_v, repeat([nothing], length(path_v)), trait_v, chr_v, tss_v, columns, separator)

QTLStudy(path::String,
    trait_v::Union{AbstractVector{String}, DatasetColumn{Dataset, Vector{Union{Missing, String}}}},
    chr_v,
    tss_v,
    columns::Union{Dict{Int, Any}, Dict{Int, GenVarInfo}},
    separator::Union{Char, AbstractVector{Char}}) = QTLStudy([path], [nothing], trait_v, chr_v, tss_v, columns, separator)


function QTLStudy_from_pattern(folder::String,
    path_pattern::Vector, 
    trait_v, 
    chr_v, 
    tss_v, 
    columns::Union{Dict{Int, Any}, Dict{Int, GenVarInfo}}, 
    separator::Union{Char, AbstractVector{Char}},
    only_corresp_chr::Bool = true)::QTLStudy                 ########### only_corresp_chr not used! (always treated as true)

    if path_pattern[1] isa String && startswith(path_pattern[1], Base.Filesystem.path_separator)
        @warn "The first element of path_pattern starts with the path separator. Paths separators are not necessary in this place. It will be removed"
        path_pattern_=[path_pattern[1][length(Base.Filesystem.path_separator)+1:end]; path_pattern[2:end]]
    else
        path_pattern_ = path_pattern
    end

    if only_corresp_chr == true && !(TRAIT_NAME in path_pattern_)
        only_corresp_chr = false
    end

    new_arr::Vector{String} = map(x -> (x isa String) ? x : "*", path_pattern_)
    
    pattern_str = accumulate(*, new_arr)[end]
    patt = Glob.GlobMatch(pattern_str)
    files = glob(patt, folder)
    
    trait_index = findfirst(x->x==TRAIT_NAME, path_pattern_)
    chr_index = findfirst(x -> x==CHR, path_pattern_)

    #regex version of pattern
    pattern_v = map(x -> (x isa String) ? raw""*x : r"(.*?)", path_pattern_)
    if endswith(folder, Base.Filesystem.path_separator)
        pattern = folder * accumulate(*, pattern_v)[end]
    else
        pattern = folder * Base.Filesystem.path_separator * accumulate(*, pattern_v)[end]
    end
    
    #verif trait in trait_v
    files_traits = nothing
    if !(trait_index isa Nothing) && trait_index != 1
        count = 1
        for i in 1:trait_index-1
            if path_pattern_[i] isa GenVarInfo
                count += 1
            end
        end
        f(x) = match(pattern, x).captures[count]
        files_traits = map(f, files)

    elseif trait_index == 1
        count = 1
        f2(x) = match(pattern, x).captures[count]
        files_traits = map(f2, files)
    else
        files_traits_b = [true for i in 1:length(files)]
    end

    #verify good chr
    files_chr = nothing
    if only_corresp_chr
        if !(chr_index isa Nothing) && chr_index != 1
            count = 1
            for i in 1:chr_index-1
                if path_pattern_[i] isa GenVarInfo
                    count += 1
                end
            end
            g(x) = parse(Int, match(pattern, x).captures[count])
            files_chr = map(g, files)
    
        elseif trait_index == 1
            g2(x) = parse(Int, match(pattern, x).captures[1])
            files_chr = map(g2, files)
        else
            files_chr_b = [true for i in 1:length(files)]
        end
    else
        files_chr_b = [true for i in 1:length(files)]
    end
    if !(files_traits isa Nothing)
        trait_index_files = [findfirst(x->(x==y), trait_v) for y in files_traits] # index of TRAIT in trait_v, and chr_v, tss_v for tss
        files_traits_b = map(x->(!(x isa Nothing)), trait_index_files)           # array of files ok for TRAIT condition
    end
    if !(files_chr isa Nothing) && only_corresp_chr
        files_chr_b = [(trait_index_files[i] isa Nothing) ? false : (files_chr[i] == chr_v[trait_index_files[i]]) for i in 1:length(files)]
    end

    indexes_keep_file = files_traits_b .& files_chr_b

    files = files[indexes_keep_file]    #vector of kept files
    if !(files_traits isa Nothing)
        traits_for_each_file = trait_v[trait_index_files[indexes_keep_file]]   #vectors of traits corresponding to kept files
    else
        traits_for_each_file = repeat([nothing], length(files)) # treat case when files are not trait_specific
    end

    return QTLStudy(files, traits_for_each_file, trait_v, chr_v, tss_v, columns, separator)
end


# Iteration overload for QtlStudy
function Base.iterate(iter::QTLStudy)
    element = GWAS(iter.path_v[1], iter.columns, iter.separator, iter.trait_v[1])
    return (element, 1)
end


function Base.iterate(iter::QTLStudy, state)
    count = state + 1
    if count > length(iter.path_v)
        return nothing
    end
    element = GWAS(iter.path_v[count], iter.columns, iter.separator, iter.traits_for_each_path[count])
    return (element, count)
end

function Base.getindex(iter::QTLStudy, i::Int)
    return GWAS(iter.path_v[i], iter.columns, iter.separator, iter.traits_for_each_path[i])
end

function Base.getindex(iter::QTLStudy, i::Union{AbstractUnitRange, AbstractVector{Int}})                 ########## Or QTLStudy subset?
    return [GWAS(iter.path_v[j], iter.columns, iter.separator, iter.traits_for_each_path[j]) for j in i]
    
end

function Base.lastindex(iter::QTLStudy)
    return lastindex(iter.path_v)
end

function Base.length(iter::QTLStudy)
    return length(iter.path_v)
end
