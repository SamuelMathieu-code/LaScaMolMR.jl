data_dir_plink_files = joinpath(pwd(), "clump_data", "sample_chosen")

chosen_sample = [(1, 87917746),
                 (1, 100046246),
                 (1, 170645774),
                 (1, 201746768),
                 (0, 0)]


@testset "clump" begin
    expected_1e_1 = [1, 1, 1, 1, 0]
    expected_2e_3 = [1, 0, 1, 1, 0]
    expected_1e_3 = [1, 0, 0, 0, 0]
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))

    @test expected_1e_1 == clump(data, chosen_sample)
    @test expected_1e_3 == clump(data, chosen_sample, r2_tresh = 0.001)
    @test expected_2e_3 == clump(data, chosen_sample, r2_tresh = 0.002)
    
end
