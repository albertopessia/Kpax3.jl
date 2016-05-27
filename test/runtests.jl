# This file is part of Kpax3. License is MIT.

# TODO: refactor test files because it is a mess

using Base.Test

include("../src/boot.jl")

ε = 1.0e-13
srand(1427371200)

tests = ["data_processing/data_processing";
         "misc/basic_functions";
         "misc/partition_functions";
         "types/settings";
         "types/data";
         "types/prior_col";
         "types/prior_row";
         "types/state";
         "types/state_list";
         "types/support";
         "model/likelihoods";
         "model/partition_cols";
         "model/partition_rows";
         "model/loss_binder";
         "optimizer/local_mode";
         "optimizer/selection";
         "optimizer/crossover";
         "optimizer/mutation";
         "mcmc/functions";
         "mcmc/weight";
         "mcmc/partition_ratios";
         "mcmc/log_likelihoods";
         "mcmc/biased_random_walk";
         "mcmc/posterior";
         "estimate/write"
]

for t in tests
  f = string(t, ".jl")
  @printf("Going through tests in '%s'... ", f)
  include(f)
  @printf("PASSED!\n")
end
