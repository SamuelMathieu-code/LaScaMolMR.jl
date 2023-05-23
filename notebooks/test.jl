include("../src/inputs.jl")

file_pat = ["test/data/exposure", TRAIT_NAME, "/exposure", TRAIT_NAME, "_chr", CHR, ".txt"]
trait_v = ["A", "C"]
chr_v = [1, 3]
tss_v = [398576, 9874365]
columns = Dict(1=>TRAIT_NAME, 2=>"etc...")
sep = ' '

QTLStudy_from_pattern(file_pat, trait_v, chr_v, tss_v, columns, ' ')
