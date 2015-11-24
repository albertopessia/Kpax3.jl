# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

tests = ["dataprocessing",
         "priorpartition"]

for t in tests
  include(t * ".jl")
end
