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
include("types/types.jl")

########################
# Load basic functions #
########################
include("initpartition.jl")

#################
# Load the rest #
#################
include("dataprocessing/dataprocessing.jl")
include("model/model.jl")
include("mcmc/mcmc.jl")
include("kpax3.jl")
