
using Glob
using InMemoryDatasets

##################################################
#                   GenVarInfo                   #
##################################################

macro exportinstances(enum)
    eval = GlobalRef(Core, :eval)
    return :($eval($__module__, Expr(:export, map(Symbol, instances($enum))...)))
end


"""
Enum of different information present in qtl file path or inside QTL/GWAS text delimited files.

```
TRAIT_NAME  # Name of the trait (e.g. QTL protein name)
CHR         # chromosome
POS         # position in chromosome
A_EFFECT    # effect allele of SNP
A_OTHER     # reference allele of SNP
BETA        # effect size of SNP
SE          # standard error of effect size
PVAL        # Pvalue for H₀ : BETA = 0
```
"""
@enum GenVarInfo begin
    TRAIT_NAME  # Name of the trait (e.g. QTL protein name)
    CHR         # chromosome
    POS         # position in chromosome
    A_EFFECT    # effect allele of SNP
    A_OTHER     # reference allele of SNP
    BETA        # effect size of SNP
    SE          # standard error of effect size
    PVAL        # Pvalue for H₀ : BETA = 0
end

@exportinstances GenVarInfo


###################################################
#                      Inputs                     #
###################################################

"""
Type GWAS which contains informations about GWAS file

```
path        # path to file
columns     # Dictionary of informations contained in columns of file 
separator   # column separator
trait_name  # name of the trait (optional)
```
## Example

```
gwas = GWAS("path/to/file", 
            Dict(1 => PVAL, 2 => CHR, 3 => POS, 8 => BETA, 9 => SE),
            ',')
```
"""
mutable struct GWAS{T <: Integer, S <: AbstractString}
    path::S
    columns::Union{Dict{T, GenVarInfo}}
    separator::Union{Char, Vector{Char}}
    trait_name::Union{Nothing, S} 
end

GWAS(path,
     columns,
     separator) = GWAS(path, columns, separator, nothing)


"""
Type QTLStudy which contains informations about QTL file(s) format and implacement.

The method [`QTLStudy_from_pattern`](@ref) helps building information from patterns in file names and is the prefered methods to construct a QTLStudy.
"""
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


"""
Make QTLStudy which contains informations about QTL file(s) format and implacement from some file pattern

arguments :

`folder` is the main folder cntaining all QTL files. \\
`path_pattern` is a vector of strngs and GenVarInfo characterizing the pattern of every file name in the folder. \\
`trait_v` is an collection of traits that should be incuded in future analysis. \\
`chr_v` is a vector of each chromosome of each QTL trait. \\
`tss_v` is a vector of every Transcription start site. \\
`columns` is a dictionary of all informations contained in columns of QTL files. \\
`separator` is the separator of QTL files. \\
`only_corresp_chr` indicates if only variants on chromosome correspoding to QTL trait should be kept. Default is `true`

TIP : if your file contains a chromosome:position column (e.g. 1:324765), consider setting your separator to `[':', sep]`

## Examples

for a sigle file QTL :

```
QTLStudy_from_pattern("some/folder", ["qtl.txt"], tv, cv, tssv,
                      Dict(1 => CHR, 2 => POS, 8 => PVAL, ...),
                      separator = ' ',
                      only_corresp_chr = true)
```

for this arhitecture (in UNIX) :

```
test/data
├── exposureA
│   ├── exposureA_chr1.txt
│   └── exposureA_chr2.txt
├── exposureB
│   ├── exposureB_chr1.txt
│   └── exposureB_chr2.txt
└── exposureC
    ├── exposureC_chr3.txt
    └── exposureC_chr4.txt
```

we get this code :

```
path_pattern = ["exposure", TRAIT_NAME, "/exposure", TRAIT_NAME, "_chr", CHR, ".txt"]
trait_v = ["A", "B", "C"]
chr_v = [1, 2, 3]
tss_v = [45287, 984276, 485327765]

qtl = QTLStudy_from_pattern("test/data", path_pattern, ttrait_v, chr_v, tss_v,
                            Dict(1 => CHR, 2 => POS, 8 => PVAL, ...),
                            separator = ' ',
                            only_corresp_chr = true)
```
"""
function QTLStudy_from_pattern(folder::AbstractString,
    path_pattern::AbstractVector, 
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

# Divide QTLStudy in N sub qtl data
"""
Partitionate QTLStudy in n folds. Returns a vector of QTLStudy in which each element contains a subsets of file paths and corresponding traits. 

**arguments :""

`x::QTLStudy` : the qtl files and informations (see [`QTLStudy`](@ref)) \\
`m::Integer` : the number of folds
"""
function nfolds(x::QTLStudy, n::Integer)
    if n > length(qtl) throw(ArgumentError("n should be smaller the qtl's nuber of files.")) end
    s = length(x) / n
    [x[round(Int64, (i-1)*s)+1:min(length(x),round(Int64, i*s))] for i=1:n]
end

# Iteration and indexing overload for QtlStudy
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

function Base.getindex(iter::QTLStudy, i::Union{AbstractUnitRange, AbstractVector{Int}})
    return QTLStudy(iter.path_v[i], iter.traits_for_each_path[i], iter.trait_v, iter.chr_v, iter.tss_v, iter.columns, iter.separator)
    
end

function Base.lastindex(iter::QTLStudy)
    return lastindex(iter.path_v)
end

function Base.length(iter::QTLStudy)
    return length(iter.path_v)
end

function Base.eachindex(iter::QTLStudy)
    return eachindex(iter.path_v, iter.traits_for_each_path)
end
