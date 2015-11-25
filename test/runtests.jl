# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

tests = ["dataprocessing",
         "datatypes",
         "likelihoods",
         "posterior",
         "priorcol",
         "priorrow"]

for t in tests
  include(string(t, ".jl"))
end
