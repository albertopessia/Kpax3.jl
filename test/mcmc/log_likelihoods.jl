# This file is part of Kpax3. License is MIT.

function test_mcmc_logmarglikmerge()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile)

  data = UInt8[0 0 0 0 0 1;
               1 1 1 1 1 0;
               0 0 1 0 1 1;
               1 1 0 1 0 0;
               1 1 0 0 0 0;
               0 0 0 1 1 0;
               1 1 1 0 0 0;
               0 0 0 1 1 1;
               0 0 1 0 0 0;
               1 0 0 1 0 1;
               0 1 0 0 1 0;
               0 0 0 0 0 1;
               1 1 1 0 0 0;
               0 0 0 1 1 0;
               1 1 0 0 1 1;
               0 0 1 1 0 0;
               1 1 0 1 0 0;
               0 0 1 0 1 1]

  (m, n) = size(data)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

  # [1; 1; 1; 2; 3; 3] => [1; 1; 1; 3; 3; 3]
  state1 = AminoAcidState(data, [1; 1; 1; 2; 3; 3], priorR, priorC, settings)
  support1 = MCMCSupport(state1, priorC)

  ni = vec(sum(float(data[:, [4; 5; 6]]), 2))
  vi = 3

  updatelogmargliki!(ni, vi, priorC, support1)

  logmarglikmerge!(state1.cl, state1.k, 3, 2, priorC, support1)

  state2 = AminoAcidState(data, [1; 1; 1; 2; 2; 2], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  @test_approx_eq_eps support1.logmlikcandidate support2.logmlik ε

  nothing
end

test_mcmc_logmarglikmerge()

function test_mcmc_logmargliksplit()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile)

  data = UInt8[0 0 0 0 0 1;
               1 1 1 1 1 0;
               0 0 1 0 1 1;
               1 1 0 1 0 0;
               1 1 0 0 0 0;
               0 0 0 1 1 0;
               1 1 1 0 0 0;
               0 0 0 1 1 1;
               0 0 1 0 0 0;
               1 0 0 1 0 1;
               0 1 0 0 1 0;
               0 0 0 0 0 1;
               1 1 1 0 0 0;
               0 0 0 1 1 0;
               1 1 0 0 1 1;
               0 0 1 1 0 0;
               1 1 0 1 0 0;
               0 0 1 0 1 1]

  (m, n) = size(data)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

  # [1; 1; 1; 2; 2; 2] => [1; 1; 1; 2; 3; 3]
  state1 = AminoAcidState(data, [1; 1; 1; 2; 2; 2], priorR, priorC, settings)
  support1 = MCMCSupport(state1, priorC)

  support1.ni = vec(sum(float(data[:, 4]), 2))
  support1.vi = 1

  support1.nj = vec(sum(float(data[:, [5; 6]]), 2))
  support1.vj = 2

  updatelogmargliki!(priorC, support1)
  updatelogmarglikj!(priorC, support1)

  logmargliksplit!(state1.cl, state1.k, 2, priorC, support1)

  state2 = AminoAcidState(data, [1; 1; 1; 2; 3; 3], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  @test_approx_eq_eps support1.logmlikcandidate support2.logmlik ε

  nothing
end

test_mcmc_logmargliksplit()
