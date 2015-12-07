# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

tests = ["dataprocessing";
         "typescluster";
         "typessettings";
         "typessupport";
         "typesdata";
         "typespartitioncols";
         "typespartitionrows";
         "typesmcmc";
         "likelihoods";
         "partitioncols";
         "partitionrows";
         "posterior"]

for t in tests
  include(string(t, ".jl"))
end
