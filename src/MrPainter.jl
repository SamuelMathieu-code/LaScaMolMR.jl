module MrPainter

# Write your package code here.
export GenVarInfo, GWAS, QtlPathPattern, QTLStudy, QTLStudy_from_pattern
export mr_output
export mr_egger, mr_ivw, mr_wald
export ld_r², mat_r²
# ... others to come

include("ld.jl")
include("inputs.jl")
include("mrPerf.jl")
include("utils.jl")
include("naiveCis.jl")
include("mrStudyCis.jl")
# ... others to come
end
