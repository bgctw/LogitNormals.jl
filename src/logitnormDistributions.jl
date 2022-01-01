# removed from Distribution.jl to here to avoid dependencies

mean(d::LogitNormal; kwargs...) = estimateMean(d,kwargs...) # moved to NormalTransforms

function mode(d::LogitNormal)
    (μ, σ) = params(d)
    # if mu<0 then maximum is left of median, if mu>0 right of median
    # if mu=0 the maximum is either at mu or there are two maxima of the same height
    # here, report the left mode
    interval = μ <= 0.0 ? (0.0,logistic(μ)) : (logistic(μ),1.0)
    resOpt = optimize(x -> -pdf(d, x), interval[1], interval[2])
    Optim.minimizer(resOpt)
end

function var(d::LogitNormal; mean=missing, kwargs...) 
    # estimateVariance does not pass kwargs to mean, need to do it here
    m = ismissing(mean) ? estimateMean(d,kwargs...) : mean
    estimateVariance(d; kwargs..., mean=m)
end

### Helpers

# opjective function minimized in mode(d)
# function objectiveModeLogitNormal(mode::Real, mu::Real, sd2::Real)
#     mode <= 0.0 || mode >= 1.0 && return prevfloat(T(Inf))
#     diff =  mu - sd2*(1 - 2*mode) - logit(mode)
#     diff^2
# end
  

