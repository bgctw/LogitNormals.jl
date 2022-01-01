function fit_mode_quantile(::Type{LogitNormal}, mode::Real, qp::QuantilePoint)
    matchModeUpper(mode, qp.q, Val(40); perc = qp.p)
end

function matchModeUpper(
    mode::T, upper::T
    , ::Val{nTry}
    ; perc::Real = 0.99
) where {nTry, T<:Real}
    mode == 0.5 && return matchMedianUpper(LogitNormal, 0.5, upper; perc=perc)
    # for given mu we can compute sigma by mode and upper quantile
    # hence univariate search for mu
    # we now that mu is in (\code{logit(mode)},0) for \code{mode < 0.5} and in 
    # \code{(0,logit(mode))} for \code{mode > 0.5} for unimodal distribution
    # there might be a maximum in the middle and optimized misses the low part
    # hence, first get near the global minimum by a evaluating the cost at a 
    # grid that is spaced narrower at the edge
    #
    logitMode = logit(mode)
    logitUpper = logit(upper)
    upperMu = abs(logitMode) - eps() 
    muTry = SVector{nTry}(
        sign(logitMode) .* log.(range(1,stop=exp(upperMu),length=nTry)))
    oF(mu) = ofLogitNormalModeUpper(mu, mode, logitMode, logitUpper, perc)
    ofMuTry = oF.(muTry)
    iMin = argmin(ofMuTry)
    # on postive side muTry are increasing, on negative side muTry decreasing
    # neet to have the lower value at the beginning of the interval
    interval = (logitMode >= 0) ? 
        (muTry[max(1,iMin-1)], muTry[min(nTry,iMin+1)]) :
        (muTry[max(1,iMin+1)], muTry[min(nTry,max(1,iMin-1))])
    resOpt = optimize(oF, interval...)
    @show ofMuTry, iMin, interval, resOpt
    Optim.converged(resOpt) || error("could not find minimum")
    μ = Optim.minimizer(resOpt)
    σ = sqrt((logitMode - μ)/(2*mode - 1))
    LogitNormal{T}(μ,σ)
end

function matchModeUpper(
    mode::S, upper::T, ::Val{nTry}; kwargs...
) where {nTry, T<:Real, S<:Real} 
    p = promote(mode, upper)
    matchModeUpper(p[1], p[2], Val(nTry); kwargs...)
end


"objective function used by `matchModeUpper(LogitNormal,...)`"
function ofLogitNormalModeUpper(mu, mode, logitMode, logitUpper, perc)
  # given mu and mode, we can calculate sigma, 
  # predict the percentile of logitUpper
  # and return the squared difference as cost to be minimized
  sigma2 = (logitMode - mu)/(2.0*mode - 1.0)
  sigma2 > 0.0 || return prevfloat(Inf)
  predp = normcdf(mu, sqrt(sigma2), logitUpper)
  diff = predp - perc
  diff*diff 
end

function fit_mode_flat(::Type{LogitNormal}, mode::T, ::Val{nTry} = Val(40)) where {T<:Real,nTry}
    # broader distributions are most flat at the edges, i.e. most close to zero
    perc = 1.0 - 1e-6
    upper = 0.0
    mode == 0.5 && return matchMedianUpper(LogitNormal, 0.5, upper; perc=perc)
    # handle case with mode left of 0.5 by afterwards providing negative mu
    is_mode_left = mode < 0.5 
    moder = is_mode_left ? (1-mode) : mode 
    logitMode = logit(moder)
    upperMu = abs(logitMode) - eps() 
    muTry = SVector{nTry}(
        sign(logitMode) .* log.(range(1,stop=exp(upperMu),length=nTry)))
    oF(mu) = ofLogitNormalModeUpperFlat(mu, moder, logitMode)
    ofMuTry = oF.(muTry)
    iMin = argmin(ofMuTry)
    interval = (muTry[max(1,iMin-1)], muTry[min(nTry,iMin+1)]) 
    resOpt = optimize(oF, interval...)
    @show ofMuTry, iMin, interval, resOpt
    Optim.converged(resOpt) || error("could not find minimum")
    μr = Optim.minimizer(resOpt)
    σ = sqrt((logitMode - μr)/(2*moder - 1))
    μ = is_mode_left ? -μr : μr
    LogitNormal{T}(μ,σ)
end


function ofLogitNormalModeUpperFlat(mu, xmode, logitMode)
    # given mu and mode, we can calculate sigma, 
    # predict the percentile of logitUpper
    # and return the squared difference as cost to be minimized
    #
    # for flat we optimize
    logitUpper = 0.0
    perc = 1.0 - 1e-6
    sigma2 = (logitMode - mu)/(2.0*xmode - 1.0)
    sigma2 > 0.0 || return prevfloat(Inf)
    predp = normcdf(mu, sqrt(sigma2), logitUpper)
    diff_pupper = predp - perc
    # And we add a penalty for negative slopes left of the mode
    # to prevent bimodal distributions
    # Hence it only works for the case of mode > 0.5
    # For the other case negate mu: (1-X) ~ logitnormal(-mu, sigma)
    # logitx = logit.([0.05,0.1,0.2,0.3,0.4,0.5])
    logitx = @SVector [-3.0, -2.2, -1.4, -0.8, -0.4, 0.0]
    px = normcdf.(mu, sqrt(sigma2), logitx)
    dpx = Base.diff(px)
    penalty_negslope = 1e10*min(0,minimum(dpx))^2
    diff_pupper*diff_pupper + penalty_negslope
end
  

