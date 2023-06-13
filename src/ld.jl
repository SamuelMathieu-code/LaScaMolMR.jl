using SnpArrays
using IterTools
using LinearAlgebra

########## Genotypes in Plink .bed format #########

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

############ LD calculations ############

"""
LD r² composite of pair of SNPs
"""
function ld_r²(snp1::Vector{UInt8}, snp2::Vector{UInt8})::Float64
    # drop missing values
    idx = (snp1 .!= @MISSING) .& (snp2 .!= @MISSING)
    s1, s2 = snp1[idx], snp2[idx]
    
    # number of not missing samples
    n = length(idx)
    
    # calculate needed frequencies
    naa = naA = nAA = nbb = nbB = nBB =
		nAABB = naabb = naaBB = nAAbb = 0

    for (g1, g2) in zip(s1, s2)
        if  g1 == @HOMO_M
            nAA +=1
            if g2 == @HOMO_M # AABB
                nBB +=1
                nAABB +=1
            elseif g2 == @HOMO_m # AAbb
                nbb +=1
                nAAbb +=1
            else    #HETER -> AAbB
                nbB +=1
            end
        elseif g1 == @HOMO_m
            naa +=1
            if g2 == @HOMO_M # aaBB
                nBB +=1
                naaBB +=1
            elseif g2 == @HOMO_m  # aabb
                nbb +=1
                naabb +=1
            else    # HETER -> aabB
                nbB +=1
            end
        else # @HETER
            naA +=1
            if g2 == @HOMO_M # aABB
                nBB +=1
            elseif g2 == @HOMO_m # aAbb
                nbb+=1
            else    # HETER -> aAbB
                nbB +=1
            end
        end
    end

    # final calculations
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

    return Δ^2 / t

end


"""
LD r² composite matrix of n SNPs
"""
function mat_r²(arr::SnpArray, idx::AbstractVector{Int})::Matrix{Float64}
    M_corr = Matrix{Float64}(I, length(idx), length(idx))
    
    d = Dict(zip(idx, 1:lastindex(idx)))
    for (i, j) in subsets(idx, 2)
        snp1 = arr[:,i]
        snp2 = arr[:,j]
        M_corr[d[i], d[j]] = M_corr[d[j], d[i]] = ld_r²(snp1, snp2) #function implemented following paper pmid :18757931, for r² type r\^2
    end
    return M_corr
end


"""
Get correlation Matrix for specifies snps (tuple of the form (chr, pos) 
    given reference genotype SnpData.
"""
function getLDmat(ref_genotypes::SnpData, 
                  snps::AbstractVector{Tuple{Any, Any}}
                  )::Tuple{Matrix{Float64}, Vector{Bool}}

    snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    for (i, chr_pos_sing) in enumerate(snps)
        local j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        if ref_genotypes.snp_info.chr_pos[j] != chr_pos_sing
            j = -1
        end
        snps_indx[i] = j
    end
    kept_indx = snps_indx[snps_indx .> 0]

    return mat_r²(ref_genotypes.snparray, kept_indx), snps_indx .> 0
end


"""
Implementation of the clumping algoritm prioritising first snps in given Vector
    returns a vector of booean indication if each given snp is kept
"""
function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{Tuple{Any, Any}}, 
               r2_tresh::Float64 = 0.1
               )::Vector{Bool}
    
    r2_mat, indx_v_b = getLDmat(ref_genotypes, snps)
    idx_on_mat = accumulate(+, indx_v_b)
    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            for j in (i+1):lastindex(snps)
                if r2_mat[idx_on_mat[i], idx_on_mat[j]] > r2_tresh
                    indx_v_b[j] = false
                end
            end
        end
    end

    return indx_v_b
end

