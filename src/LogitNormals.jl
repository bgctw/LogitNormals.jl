module LogitNormals

using Distributions, DistributionFits
using StaticArrays
using Optim
import StatsFuns: logit, logistic, normcdf

# for extension
import Distributions: mean, var, mode
import StatsBase: fit
import DistributionFits: fit_mode_quantile, fit_mean_quantile

export fit_mode_flat

include("logitnormal.jl")
include("util.jl")
include("logitnormDistributions.jl")

end
