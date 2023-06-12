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

