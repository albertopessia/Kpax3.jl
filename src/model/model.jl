# This file is part of Kpax3. License is MIT.

include("joint/densities.jl")

include("likelihoods/marglik.jl")
include("likelihoods/merge.jl")
include("likelihoods/split.jl")
include("likelihoods/biasedrw.jl")

include("partitioncols/aminoacids/densities.jl")
include("partitioncols/aminoacids/simulate.jl")

include("partitionrows/ewenspitman/densities.jl")
include("partitionrows/ewenspitman/split.jl")
include("partitionrows/ewenspitman/merge.jl")
include("partitionrows/ewenspitman/biasedrw.jl")
