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

function is_logit_slope_monotone(d, upper=0.5, lower=0.0, decreasing = false) 
  x = range(lower,upper,41)[2:40]	# plotting grid
  dx = pdf.(d, x)	#density function
  if decreasing
    all(diff(dx) .<= 0)
  else
    all(diff(dx) .>= 0)
  end
end

@testset "fit by single mode and flat" begin
    d9 = d = fit_mode_flat(LogitNormal, 0.9)
    @test isapprox( mode(d), 0.9, atol=1e-4)
    @test is_logit_slope_monotone(d, 0.9)
    @test is_logit_slope_monotone(d, 0.0, 0.9, true)
    d1 = d = fit_mode_flat(LogitNormal, 0.1)
    @test isapprox( mode(d), 0.1, atol=1e-4)
    @test is_logit_slope_monotone(d, 0.1)
    @test is_logit_slope_monotone(d, 0.0, 0.1, true)
    @test d1.σ == d9.σ
    @test d1.μ == -d9.μ
    d5 = d = fit_mode_flat(LogitNormal, 0.5)
    @test isapprox( mode(d), 0.5, atol=1e-4)
    @test is_logit_slope_monotone(d, 0.5)
    @test is_logit_slope_monotone(d, 0.0, 0.5, true)
    #plot(d1)
end