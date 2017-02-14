# This file is part of Kpax3. License is MIT.

#################
# Load packages #
#################
import Clustering
import Distances
import FileIO
import JLD
import Gadfly
import Measures
import StatsBase

##############
# Load Types #
##############
include("types/types.jl")

########################
# Load basic functions #
########################
include("misc/misc.jl")

#################
# Load the rest #
#################
include("data_processing/data_processing.jl")
include("distances/distances.jl")
include("model/model.jl")
include("optimizer/optimizer.jl")
include("mcmc/mcmc.jl")
include("estimate/estimate.jl")
include("plots/plots.jl")
