# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

ε = 1.0e-13

tests = ["dataprocessing";
         "typessettings";
         "typessupport";
         "typesdata";
         "typespartitioncols";
         "typespartitionrows";
         "typesmcmc";
         "likelihoods";
         "partitioncols";
         "partitionrows";
         "mcmcfunctions";
         "mcmcpartitionratio";
         "mcmcloglik";
         "mcmcbrw";
         "posterior"]

for t in tests
  include(string(t, ".jl"))
end
