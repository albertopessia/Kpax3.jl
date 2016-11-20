# This file is part of Kpax3. License is MIT.

#################
# Load packages #
#################
import StatsBase: sample, sample!, WeightVec, values
import Clustering: kmedoids
import Distances: hamming
import FileIO: File, @format_str
import JLD: load, save

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
