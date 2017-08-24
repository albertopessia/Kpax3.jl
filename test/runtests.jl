# This file is part of Kpax3. License is MIT.

importall Base.Test
import Kpax3

cd(dirname(@__FILE__))

ε = 1.0e-13
srand(1427371200)

function runtests()
  # this module is used by "model/partition_rows" tests
  include("data/partitions.jl")

  tests = ["data_processing/fasta_data_processing";
           "data_processing/csv_data_processing";
           "distances/simovici_jaroszewicz";
           "distances/jaccard";
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
           "model/loss_binder"
           "optimizer/local_mode";
           "optimizer/selection";
           "optimizer/crossover";
           "optimizer/mutation";
           "mcmc/partition_ratios";
           "mcmc/log_likelihoods";
           "mcmc/weight";
           "mcmc/merge";
           "mcmc/split";
           "mcmc/gibbs";
           "mcmc/biased_random_walk";
           "mcmc/posterior";
           "mcmc/diagnostics";
           "estimate/write"
  ]

  for t in tests
    f = string(t, ".jl")
    @printf("Going through tests in '%s'... ", f)
    include(f)
    @printf("PASSED!\n")
  end

  nothing
end

runtests()
