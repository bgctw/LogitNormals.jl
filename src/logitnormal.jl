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
    mode == 0.5 && return(LogitNormal{T}(0.0, sqrt(2)))
    is_right = mode > 0.5
    mode_r = is_right ? mode : 1.0 - mode
    res_opt = optimize(x -> of_mode_flat(x, mode_r, logit(mode_r)), 0.0, 0.5)
    Optim.converged(res_opt) || error("could not find minimum")
    xt = res_opt.minimizer
    σ2 = (1/xt + 1/(1-xt))/2
    μr = logit(mode_r) - σ2*(2.0*mode_r - 1.0)
    μ = is_right ? μr : -μr
    LogitNormal{T}(μ,sqrt(σ2))
end

function of_mode_flat(x, m, logitm = logit(m))
  # lhs = 1/x + 1/(1-x)
  # rhs = (logitm - logit(x))/(m-x)
  # lhs = (m-x)/x + (m-x)/(1-x)
  # rhs = (logitm - logit(x))
  lhs = m/x + (m-1)/(1-x)
  rhs = logitm - logit(x)
  d = lhs - rhs
  d*d
end



  

