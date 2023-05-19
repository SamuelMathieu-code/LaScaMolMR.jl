
# DATA

path = "some/path"
columns_1 = Dict(1=>TRAIT_NAME, 2=>CHR, 3=>POS, 4=>BETA, 5=>SE, 6=>PVAL)


sep = ' '
trait_name = "trait"
chr = 9
tss = 8764

gwas_1 = GWAS(path, columns_1, sep)
gwas_2 = GWAS(path, columns_1, sep, trait_name)

qtl = QtlStudy([path, path], [trait_name, trait_name], [chr, chr], [tss, tss], columns_1, sep)

# TESTSET

@testset "inputs.jl" begin

    for gwas in qtl
        @test gwas.separator == ' '
        @test gwas.path == path
    end
    
end