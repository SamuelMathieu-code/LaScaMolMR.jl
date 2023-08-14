using SnpArrays # To read PLINK 1.9 .bim .bed .fam files. 
import IterTools.subsets
using LinearAlgebra
using Base.Threads

########## Genotypes in Plink .bed format #########

### Macros describing the PLINK 1.9 .bed file format are for code readability
macro HOMO_M()  # AA
    return :(0x00)
end
macro MISSING()
    return :(0x01)
end
macro HETER()   # aA / Aa
    return :(0x02)
end
macro HOMO_m()  # aa
    return :(0x03)
end

################## MAF ##################

@inline function maf(s1)::Float64
    m = n = N = 0
    @inbounds for i in eachindex(s1)
        g = s1[i]
        if g == @HOMO_M
            N += 2
            m += 2
        elseif g == @HOMO_m
            n += 2
            m += 2
        elseif g == @HETER
            n += 1
            N +=1
            m += 2
        end
    end

    return min(N/m, n/m)
end

############ LD calculations ############


# LD r² composite of pair of SNPs from two 1 dimension BitArrays
@inline function ld_r²(s1, s2)::Float64
    n = 0
    
    # calculate needed frequencies
    naa = naA = nAA = nbb = nbB = nBB =
		nAABB = naabb = naaBB = nAAbb = 0

    @inbounds @simd for i in eachindex(s1, s2)
        g1, g2 = s1[i], s2[i]
        if  g1 == @HOMO_M
            if g2 == @HOMO_M # AABB
                nBB +=1
                nAABB +=1
                n += 1
                nAA +=1
            elseif g2 == @HOMO_m # AAbb
                nbb +=1
                nAAbb +=1
                n += 1
                nAA +=1
            elseif g2 == @HETER    #HETER -> AAbB
                nbB +=1
                n += 1
                nAA +=1
            end
        elseif g1 == @HOMO_m
            if g2 == @HOMO_M # aaBB
                nBB +=1
                naaBB +=1
                n += 1
                naa +=1
            elseif g2 == @HOMO_m  # aabb
                nbb +=1
                naabb +=1
                n += 1
                naa +=1
            elseif g2 == @HETER    # HETER -> aabB
                nbB +=1
                n += 1
                naa +=1
            end
        elseif g1 == @HETER # @HETER
            if g2 == @HOMO_M # aABB
                nBB +=1
                n += 1
                naA +=1
            elseif g2 == @HOMO_m # aAbb
                nbb+=1
                n += 1
                naA +=1
            elseif g2 == @HETER    # HETER -> aAbB
                nbB +=1
                n += 1
                naA +=1
            end
        end
    end

    # final calculations
    @fastmath begin
        Δ = (nAABB + naabb - naaBB - nAAbb) / (2*n) - 
            (naa-nAA)*(nbb-nBB) / (2*n*n)
        
        pa = (2*naa + naA) / (2*n)
        pb = (2*nbb + nbB) / (2*n)
        pA = 1 - pa
        pB = 1 - pb
        pAA = nAA / n
        pBB = nBB / n
        DA = pAA - pA*pA
        DB = pBB - pB*pB
        t = (pA*pa + DA) * (pB*pb + DB)

        r = Δ^2 / t
    end

    return r
end


"""
threaded implementation of the clumping algorithm prioritising first snps in given Vector and formated Genotypes SnpData (see [`formatSnpData!`](@ref))
    returns a vector of boolean indicating if each given snp is kept

arguments :

ref_genotypes : reference SnpData ([SnpArrays](https://github.com/OpenMendel/SnpArrays.jl)) \\
snps : vector of chromosome position tuple for each variant

options : 

`formated` : indicates if ref SnpData is already formated according to :chr_pos (chr, pos)\\
`r2_tresh` : minimal r² for 2 snps to be considered strongly correlated.

## Examples :

```julia
julia> ref = SnpData(datadir("some/data"));

julia> kept_v_b::Vector{Bool} = clump([(1, 123), (1, 456), (1, 789)], ref)
3-element Vector{Int64}:
 1
 0
 1

julia> formatSnpData!(ref);

julia> kept_v_b::Vector{Bool} = clump([(1, 123), (1, 456), (1, 789)], ref, formated = true)
3-element Vector{Int64}:
 1
 0
 1
```

If formatSnpData has already been called on good snp info type (`:chr_pos` or `:snpid`), `formated = true` option does not verify or modify formating.
See [`formatSnpData!`](@ref).
"""
function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:Tuple{Integer, Integer}}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false,
               min_maf::Real = 0
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1)) # indices 
    indx_v_b = Vector{Bool}(undef, size(snps, 1)) # found or not
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsorted(ref_genotypes.snp_info.chr_pos, chr_pos_sing) # range valid if of length 1 (> 1 => triallelic, < 1 => not found)
        if length(j) == 1
            snps_indx[i] = ref_genotypes.snp_info.idx[j[]]

            # if maf insufficient : discard variant
            if maf(ref_genotypes.snparray[:, snps_indx[i]]) > min_maf
                indx_v_b[i] = true
            else
                indx_v_b[i] = false
            end
        else
            # if variant not biallelic : discard
            indx_v_b[i] = false
        end
    end

    
    for i in 1:(lastindex(snps)-1)
        if @inbounds(indx_v_b[i])
            @threads for j in i+1:lastindex(snps)
                @inbounds if indx_v_b[j]
                    s1 = @view ref_genotypes.snparray[:, snps_indx[i]]
                    s2 = @view ref_genotypes.snparray[:, snps_indx[j]]
                    if ld_r²(s1, s2) > r2_tresh                         # !(abs(snps[i][2] - snps[j][2]) > 1_000_000) && ld_r²(s1, s2) > r2_tresh 
                            indx_v_b[j] = false                         # -> do not calculate if distance is more than 1Mbp
                    end
                end
            end
        end
    end

    return indx_v_b
end


"""
format Genotype information contained in SnpData for optimised snp search based on chromosome position.
Adds a column of tuple (chr::Int8, pos::Int) in snp_info and sorts snp_info accordingly.
    returns nothing

## Examples

```julia
ref = SnpData(datadir("some/data"))

formatSnpData!(ref)
```

Note this function does not take into account the user might want to write SnpData in bed, bim, fam format. 
    The changes done by this function can and will cause problems if writing plink files after calling `formatSnpData!`.
"""
function formatSnpData!(Genotypes::SnpData)
    if !hasproperty(Genotypes.snp_info, :idx)
        Genotypes.snp_info.idx = collect(1:size(Genotypes.snp_info, 1))
    end
    if !hasproperty(Genotypes.snp_info, :chr_pos)
        Genotypes.snp_info.chr_pos = collect(
                zip(parse.(Int8, Genotypes.snp_info.chromosome), 
                    Genotypes.snp_info.position)
            )
    end
    sort!(Genotypes.snp_info, :chr_pos)
end

