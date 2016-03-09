# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

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
         "mcmcloglik"
         "mcmcbrw"]

for t in tests
  include(string(t, ".jl"))
end
