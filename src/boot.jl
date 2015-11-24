# This file is part of Kpax3. License is MIT.

#################
# Load packages #
#################
using StatsBase
using Distributions

import FileIO: File, @format_str
import JLD: load, save

##############
# Load Types #
##############
include("types/abstracts.jl")
include("types/exceptions.jl")
include("types/aacolprior.jl")
include("types/ewenspitman.jl")
include("types/data.jl")
include("types/mcmc.jl")

########################
# Load basic functions #
########################
include("initpartition.jl")

#################
# Load the rest #
#################
include("dataprocessing/dataprocessing.jl")
include("likelihoods/likelihoods.jl")
include("priors/priors.jl")
include("posteriors/posteriors.jl")
include("mcmc/mcmc.jl")
