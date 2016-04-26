# This file is part of Kpax3. License is MIT.

using Base.Test

include("../src/boot.jl")

ε = 1.0e-13
srand(1427328000)

tests = ["data_processing/data_processing";
         "types/settings";
         "types/support";
         "types/data";
         "types/partition_cols";
         "types/partition_rows";
         "types/state";
         "model/likelihoods";
         "model/partition_cols";
         "model/partition_rows";
         "model/loss_binder";
         "mcmc/functions";
         "mcmc/partition_ratios";
         "mcmc/log_likelihoods";
         "mcmc/biased_random_walk";
         "mcmc/posterior";
         "estimate/write"]

for t in tests
  begin
    include(string(t, ".jl"))
    nothing
  end
end
