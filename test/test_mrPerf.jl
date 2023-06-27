# DATA
data = (;
    βx = [-0.25, 0.35, 0.42, 0.87, 0.7],
    βy = [-0.1, 0.2, 0.22, 0.5, 0.8],
    se_βx = [0.65, 1.65, 0.8, 1.1, 0.3],
    se_βy = [0.02, 0.1, 0.4, 0.3, 0.09])

data1 = (;
    βx = -0.25,
    βy = -0.1,
    se_βx = 1.65,
    se_βy = 0.1)

data2 = (;
    βx = [-0.25, 0.35],
    βy = [-0.1, 0.2],
    se_βx = [0.65, 1.65],
    se_βy = [0.02, 0.1])

data3 = (;
    βx = [-0.25, 0.35, 0.42],
    βy = [-0.1, 0.2, 0.22],
    se_βx = [0.65, 1.65, 0.8],
    se_βy = [0.02, 0.1, 0.4])

data_error = (;
    βx = [-0.25, 0.35],
    βy = [-0.1, 0.2, 0.22],
    se_βx = [0.65, 1.65, 0.8],
    se_βy = [0.02, 0.1, 0.4])

data_perfect = (;
    βx = [0.1, 0.2, 0.3],
    βy = 2*[0.1, 0.2, 0.3],
    se_βx = [0.01, 0.01, 0.01],
    se_βy = [0.01, 0.01, 0.01])

@testset "mrPerf.jl" begin
    # All regressions
    res0_ivw = mr_ivw(data.βy, data.se_βy, data.βx, data.se_βx, 0.025)
    res0_egger = mr_egger(data.βy, data.se_βy, data.βx)
    res1 = mr_wald(data1.βy, data1.se_βy, data1.βx, data.se_βx, 0.05)
    res0_wm = mr_wm(data.βy, data.se_βy, data.βx, data.se_βx, 0.025)
    res2_ivw = mr_ivw(data2.βy, data2.se_βy, data2.βx, data2.se_βx, 0.05)
    res2_egger = mr_egger(data2.βy, data2.se_βy, data2.βx, data2.se_βx, 0.025)
    res2_wm = mr_wm(data2.βy, data2.se_βy, data2.βx, data2.se_βx, 0.025)
    res3_ivw = mr_ivw(data3.βy, data3.se_βy, data3.βx)
    res3_egger = mr_egger(data3.βy, data3.se_βy, data3.βx)
    res3_wm = mr_wm(data3.βy, data3.se_βy, data3.βx, data3.se_βx)
    res_perfect_ivw = mr_ivw(data_perfect.βy, data_perfect.se_βy, data_perfect.βx, data_perfect.se_βx)
    res_perfect_egger = mr_egger(data_perfect.βy, data_perfect.se_βy, data_perfect.βx, data_perfect.se_βx)
    res_perfect_wm = mr_wm(data_perfect.βy, data_perfect.se_βy, data_perfect.βx, data_perfect.se_βx)

    @test_throws Exception mr_egger(data_error.βy, data_error.se_βy, data_error.βx)
    @test_throws Exception mr_ivw(data_error.βy, data_error.se_βy, data_error.βx)
    @test_throws Exception mr_wm(data_error.βy, data_error.se_βy, data_error.βx, data_error.se_βx)

    @test res0_ivw.nivs == 5 && isnan(res0_ivw.intercept) && !isinf(res0_ivw.effect)

    @test res0_wm.nivs == 5 && isnan(res0_wm.intercept) && !isinf(res0_wm.effect)
    
    @test res0_egger.nivs == 5 && 
        all([!isnan(getfield(res0_egger, n)) && !isinf(getfield(res0_egger, n)) for n in fieldnames(mr_output)])
       
    @test res3_ivw.nivs == 3 && isnan(res3_ivw.intercept) && !isinf(res3_ivw.effect)
    
    @test res3_egger.nivs == 3 && 
        all([!isnan(getfield(res3_egger, n)) && !isinf(getfield(res3_egger, n)) for n in fieldnames(mr_output)])
    
    @test res2_ivw.nivs == 2 && isnan(res2_ivw.intercept) && !isinf(res2_ivw.effect)
    
    @test res2_wm.nivs == 2 && isnan(res2_wm.intercept) && !isinf(res2_wm.effect)

    @test isnan(res2_egger.heter_p) && isnan(res2_egger.heter_stat)
    
    @test !isnan(res1.p) && !iszero(res1.p) && !isinf(res1.p) && 
        !isnan(res1.ci_high) && !isinf(res1.ci_high) &&
        !isnan(res1.ci_low) && !isinf(res1.ci_low) &&
        !isnan(res1.effect) && !isinf(res1.effect)
    
    @test isnan(res1.heter_p) && isnan(res1.heter_stat) && 
        isnan(res1.intercept) && isnan(res1.p_intercept)

    @test isapprox(res_perfect_egger.effect, 2, atol = 1e-10)
    @test isapprox(res_perfect_ivw.effect, 2, atol = 1e-10)
    @test isapprox(res_perfect_wm.effect, 2, atol = 1e-10)

end
