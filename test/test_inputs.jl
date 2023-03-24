
path = "some/path"
columns_1 = Dict(TRAIT_NAME => 1, CHR => 2, POS => 3, BETA => 4, SE => 5, PVAL => 6)
columns_fails = Dict(TRAIT_NAME => 1, CHR_COLON_POS => 3, ODD_RATIO => 4, SE => 5, PVAL => 6)
columns_2 = Dict(TRAIT_NAME => 1, CHR_COLON_POS => 3, ODD_RATIO => 4, CI_LOW => 5, CI_HIGH => 6, PVAL => 7)
sep = " "
trait_name = "trait"
chr = 9
tss = 8764

gwas_1 = GWAS(path, columns_1, sep)
gwas_2 = GWAS(path, columns_2, sep, trait_name)
str_1 = "trait 9 8765 0.25 0.09 1e-6"
str_2 = "trait rs8756 9:8765 1.5 1.2 1.8 1.2e-6"
trait_1, variant_1, effect_1 = gwas_1.acc_func(str_1)
trait_2, variant_2, effect_2 = gwas_2.acc_func(str_2)

qtl = QtlStudy([path, path], [trait_name, trait_name], [chr, chr], [tss, tss], columns_1, sep)

@testset "inputs.jl" begin

    @test_throws ErrorException GWAS(path, columns_fails, sep, trait_name, chr, tss)
    @test trait_1 == "trait"
    @test trait_2 == "trait"
    @test variant_1[1] == 9
    @test variant_1[2] == 8765
    @test variant_2[1] == 9 
    @test variant_2[2] == 8765
    @test effect_1[1] == 0.25 
    @test effect_1[2] == 0.09
    @test effect_1[3] == 1e-6
    @test effect_2[3] == 1.2e-6
    @test effect_2[1] ≈ 0.8175 atol=1e-3
    @test effect_2[2] ≈ 0.0228 atol=1e-3

    for gwas in qtl
        @test gwas.trait_name == "trait"
        @test gwas.chr == 9
        @test gwas.tss == 8764
        trait_test, _, _ = gwas.acc_func(str_1)
        @test trait_test == "trait"
        @test gwas.separator == " "
    end
    
end