module MrPainter

# Write your package code here.
export GenVarInfo, GWAS, QTLStudy, QTLStudy_from_pattern, nfolds
export mr_output, mr_egger, mr_ivw, mr_wald, mr_wm
export clump, formatSnpData!
export mrStudyCis, mrStudyCisNFolds, NaiveCis
# ... others to come

include("ld.jl")
include("inputs.jl")
include("mrPerf.jl")
include("utils.jl")
include("naiveCis.jl")
include("mrStudyCis.jl")
# ... others to come
end
