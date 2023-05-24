using DLMReader
using InMemoryDatasets

###############################
#         MrStudyCis          #
###############################

function mrStudyCis(exposure::QTLStudy, outcome::GWAS, approach::String="naive", p_thresh::Float = 5e-3, window::Int = 500000, r2_tresh::Float = 0.1)::AbstractDataset
    
end