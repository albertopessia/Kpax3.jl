# This file is part of Kpax3. License is MIT.

# TODO: refactor test files because it is a mess

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
    f = string(t, ".jl")
    @printf("Going through tests in '%s'... ", f)
    include(f)
    @printf("PASSED!\n")
    nothing
  end
end
