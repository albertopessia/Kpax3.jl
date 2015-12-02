# This file is part of Kpax3. License is MIT.

include("joint/densities.jl")

include("likelihoods/marglik.jl")

include("partitioncols/aminoacids/densities.jl")
include("partitioncols/aminoacids/simulate.jl")

include("partitionrows/ewenspitman/densities.jl")
include("partitionrows/ewenspitman/split.jl")
