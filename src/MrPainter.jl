module MrPainter

# Write your package code here.
export GenVarInfo, GWAS, QtlPathPattern, QTLStudy, QTLStudy_from_pattern
export ivSelectCis, ivSelectTrans
export mr_output
export mr_egger, mr_ivw, mr_wald
# ... others to come

include("inputs.jl")
include("ivSelect.jl")
include("mrPerf.jl")
# ... others to come
end
