@testset "Logitnormal numerical moments" begin
    DN = LogitNormal(1.3)
    x = rand(DN, 1_000_000);
    m = mean(DN)
    @test abs(m - mean(x))/mean(x) <= 1e-3
    s2 = var(DN)
    @test abs(s2 - var(x))/var(x) <= 1e-2 #1e-4 too strong for random numbers
    s2m = var(DN; mean=m) # specify so that not need to recompute
    @test abs(s2m - var(x))/var(x) <= 1e-2 #1e-4 too strong for random numbers
    sd = std(DN) # specify so that not need to recompute
    @test abs(sd - std(x))/std(x) <= 1e-2 #1e-4 too strong for random numbers
    # sd2 = std(DN, mean=m) # specify so that not need to recompute
    # @test abs(sd2 - std(x))/std(x) <= 1e-2 #1e-4 too strong for random numbers
end

function test_lognormal_mode(g)
    mo = mode(g)
    @test pdf(g, mo-1e-6) < pdf(g,mo) && pdf(g, mo+1e-6) < pdf(g,mo)
end

###### numerical estimation of moments
@testset "Logitnormal seek mode" begin
    #plot(g); vline!([mode(g)])
    # median
    g = LogitNormal(); test_lognormal_mode(g)
    g = LogitNormal(0,1.4); test_lognormal_mode(g) 
    # larger one of two modes, infer the lower one
    g = LogitNormal(0,1.6); test_lognormal_mode(g)
    # two modes, larger is on the right
    g = LogitNormal(0.1,1.6); test_lognormal_mode(g)
    # two modes, larger is on the left
    g = LogitNormal(-0.1,1.6); test_lognormal_mode(g)
end

@testset "fit by mode and upper quantile" begin
    #plot(g); vline!([mode(g)])
    # median
    #g = fit_mode_quantile(LogitNormal, 0.68, @qp(0.71,0.99))
    g = fit(LogitNormal, 0.68, @qp_u(0.9), Val(:mode))
    @test mode(g) ≈ 0.68 
    @test quantile(g, 0.95) ≈ 0.9
    #plot(g)
end

#@testset "fit by single mode and flat" begin
    #plot(g); vline!([mode(g)])
    # median
    #g = fit_mode_quantile(LogitNormal, 0.68, @qp(0.71,0.99))
    xmode = 0.2
    xmode = 0.9
    #g = fit(LogitNormal, xmode, @qp(0.999,1e-4), Val(:mode))
    g = fit_mode_flat(LogitNormal, xmode)
    @test isapprox( mode(g), xmode, atol=1e-4)
    
    plot(g)
end