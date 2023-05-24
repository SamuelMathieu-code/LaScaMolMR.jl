
# DATA

path = "some/path"
columns_1 = Dict(1=>TRAIT_NAME, 2=>CHR, 3=>POS, 4=>BETA, 5=>SE, 6=>PVAL)


sep = ' '
trait_name = "trait"
chr = 9
tss = 8764


@testset "inputs.jl -- iteration qtl" begin
    gwas_1 = GWAS(path, columns_1, sep)
    gwas_2 = GWAS(path, columns_1, sep, trait_name)

    qtl = QTLStudy([path, path*"2"], [trait_name, trait_name*"2"],[trait_name, trait_name*"2"], [chr, chr], [tss, tss], columns_1, sep)
    i = 1
    for gwas in qtl
        @test gwas.separator == ' '
        if i == 1
            @test gwas.path == path
            @test gwas.trait_name == trait_name
        else
            @test gwas.path == path*"2"
            @test gwas.trait_name == trait_name*"2"
        end
        i+=1
    end
end

@testset "inputs.jl -- 1 file qtl" begin
    trait_v = ["A", "C"]
    chr_v = [1, 3]
    tss_v = [398576, 9874365]
    columns = Dict(1=>TRAIT_NAME, 2=>"etc...")
    sep = ' '
    qtl = QTLStudy("my/file", trait_v, chr_v, tss_v, columns, sep)
    @test qtl.traits_for_each_path == [nothing]
end

@testset "inputs.jl -- 1 file per chr qtl" begin
    folder = "../test/data2"
    pattern = ["chr", CHR, ".txt"]
    trait_v = ["A", "C"]
    chr_v = [1, 3]
    tss_v = [398576, 9874365]
    columns = Dict(1=>TRAIT_NAME, 2=>"etc...")
    sep = ' '

    qtl = QTLStudy_from_pattern(folder, pattern, trait_v, chr_v, tss_v, columns, sep, true)
    qtl2 = QTLStudy_from_pattern(folder, pattern, trait_v, chr_v, tss_v, columns, sep, false)

    @test length(qtl.path_v) == 4 && (nothing in qtl.traits_for_each_path) && length(qtl.traits_for_each_path) == 4
    @test length(qtl2.path_v) == 4 && (nothing in qtl2.traits_for_each_path) && length(qtl2.traits_for_each_path) == 4

end

@testset "inputs.jl -- 1 file per exposure qtl" begin
    folder = "../test/data3"
    pattern = ["exposure", TRAIT_NAME, ".txt"]
    trait_v = ["A", "C"]
    chr_v = [1, 3]
    tss_v = [398576, 9874365]
    columns = Dict(1=>TRAIT_NAME, 2=>"etc...")
    sep = ' '

    qtl = QTLStudy_from_pattern(folder, pattern, trait_v, chr_v, tss_v, columns, sep, false)
    qtl2 = QTLStudy_from_pattern(folder, pattern, trait_v, chr_v, tss_v, columns, sep, true)

    @test length(qtl.path_v) == 2 && ("A" in qtl.traits_for_each_path) && !("B" in qtl.traits_for_each_path) && length(qtl.traits_for_each_path) == 2
    @test length(qtl2.path_v) == 2 && ("A" in qtl2.traits_for_each_path) && !("B" in qtl2.traits_for_each_path) && length(qtl2.traits_for_each_path) == 2
end

@testset "inputs -- 1 file per (exposure, chr) qtl" begin
    
    file_pat = ["exposure", TRAIT_NAME, "/exposure", TRAIT_NAME, "_chr", CHR, ".txt"]
    folder = "../test/data/"
    trait_v = ["A", "C"]
    chr_v = [1, 3]
    tss_v = [398576, 9874365]
    columns = Dict(1=>TRAIT_NAME, 2=>"etc...")
    sep = ' '

    qtl = QTLStudy_from_pattern(folder, file_pat, trait_v, chr_v, tss_v, columns, ' ', false)
    qtl2 = QTLStudy_from_pattern(folder, file_pat, trait_v, chr_v, tss_v, columns, ' ')

    @test length(qtl.path_v) == 4 && !("B" in qtl.traits_for_each_path) && ("C" in qtl.traits_for_each_path) && length(qtl.traits_for_each_path) == 4
    @test length(qtl2.path_v) == 2 && !("B" in qtl2.traits_for_each_path) && ("C" in qtl2.traits_for_each_path) && length(qtl2.traits_for_each_path) == 2

end