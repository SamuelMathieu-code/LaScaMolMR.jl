using InMemoryDatasets
using DLMReader
using SnpArrays
import Base.Threads.@threads


"""
Get correlation Matrix for specifies snps (tuple of the form (chr, pos) 
    given reference genotype SnpData.
"""
function getLDmat(ref_genotypes::SnpData, 
                  snps::AbstractVector{Tuple{Int8, Int}}
                  )::Tuple{Matrix{Float64}, AbstractVector{Tuple{Int8, Int}}}

    snps_indx = Vector{Int}(undef, size(snps, 1))
    for (i, chr_pos_sing) in enumerate(snps)
        local j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        if ref_genotypes.snp_info.chr_pos[j] != chr_pos_sing
            j = NaN
        end
        snps_indx[i] = j
    end
    kept_indx = snps_indx[.!(isnan.(snps_indx))]

    return mat_r²(ref_genotypes.snparray, kept_indx), .!(isnan.(snps_indx))
end


"""
Implementation of the clumping algoritm prioritising first snps in given Vector
    returns a vector of booean indication if each given snp is kept
"""
function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{Tuple{Int8, Int}}, 
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


"""
Implmentation of naive approach for transcriptome wide MR
"""
function NaiveCis(data::Dataset, GenotypesArr::AbstractVector{SnpData}, 
                  one_file_per_chr_plink = true, 
                  r2_tresh::Float64 = 0.1,
                  mr_methoods::AbstractVector{Function} = [mr_egger, mr_ivw]
                  )::Dataset
    
    # MR pipeline for each exposure
    @threads for data_group in eachgroup(data)
        ivs_d = sort(data_group, :pval)
        if one_file_per_chr_plink
            local chr = ivs_d.chr[0]
            kept_v_b = clump(GenotypesArr[chr], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        else
            kept_v_b = clump(GenotypesArr[1], 
                             collect(zip(ivs_d.chr, ivs_d.pos)), 
                             r2_tresh)
        end
        ivs_d = ivs_d[kept_v_b, :]
        
        # make mr methods and write ouput in dataset
        # ...

    end
    
    return Dataset()
end