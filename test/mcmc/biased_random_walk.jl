# This file is part of Kpax3. License is MIT.

function test_mcmc_brw()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test.bin"

  settings = KSettings(ifile, ofile, maxclust=1, maxunit=1, op=[0.0; 1.0])

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

  # [2; 2; 2; 1; 1; 3] => [2; 2; 2; 1; 1; 1]
  state1 = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)
  support1 = MCMCSupport(state1, priorC)

  state2 = AminoAcidState(data, [2; 2; 2; 1; 1; 1], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  i = 6
  hi = state1.R[i]
  hj = state1.R[5]

  initsupportbrwmerge!(i, hj, data, support1, state1)

  support1.lograR = logratiopriorrowbrwmerge(state2.k, state1.v[hj], priorR)

  updatelogmarglikj!(priorC, support1)

  logmarglikbrwmerge!(state1.cl, state1.k, hi, hj, priorC, support1)

  performbrwmerge!(i, hi, hj, priorC, settings, support1, state1)

  @test state1.R == state2.R
  @test state1.k == state2.k

  @test !state1.emptycluster[1]
  @test state1.cl[1] == state2.cl[1]
  @test state1.v[1] == state2.v[1]
  @test state1.n1s[1, :] == state2.n1s[1, :]
  @test state1.unit[1][1:state1.v[1]] == [4; 5; 6]

  @test !state1.emptycluster[2]
  @test state1.cl[2] == state2.cl[2]
  @test state1.v[2] == state2.v[2]
  @test state1.n1s[2, :] == state2.n1s[2, :]
  @test state1.unit[2][1:state1.v[2]] == [1; 2; 3]

  @test state1.emptycluster[3]

  @test_approx_eq_eps state1.logpR state2.logpR ε

  @test_approx_eq_eps support1.logmlik support2.logmlik ε

  # [2; 2; 2; 1; 1; 3] => [2; 2; 3; 1; 1; 3]
  state1 = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)
  support1 = MCMCSupport(state1, priorC)

  state2 = AminoAcidState(data, [2; 2; 3; 1; 1; 3], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  i = 3
  hi = state1.R[i]
  hj = state1.R[6]

  initsupportbrwmove!(i, hi, hj, data, support1, state1)

  support1.lograR = logratiopriorrowbrwmove(state1.v[hi], state1.v[hj], priorR)

  updatelogmargliki!(priorC, support1)
  updatelogmarglikj!(priorC, support1)

  logmarglikbrwmove!(state1.cl, state1.k, hi, hj, priorC, support1)

  performbrwmove!(i, hi, hj, support1, state1)

  @test state1.R == state2.R
  @test state1.k == state2.k

  @test !state1.emptycluster[1]
  @test state1.cl[1] == state2.cl[1]
  @test state1.v[1] == state2.v[1]
  @test state1.n1s[1, :] == state2.n1s[1, :]
  @test state1.unit[1][1:state1.v[1]] == [4; 5]

  @test !state1.emptycluster[2]
  @test state1.cl[2] == state2.cl[2]
  @test state1.v[2] == state2.v[2]
  @test state1.n1s[2, :] == state2.n1s[2, :]
  @test state1.unit[2][1:state1.v[2]] == [1; 2]

  @test !state1.emptycluster[3]
  @test state1.cl[3] == state2.cl[3]
  @test state1.v[3] == state2.v[3]
  @test state1.n1s[3, :] == state2.n1s[3, :]
  @test state1.unit[3][1:state1.v[3]] == [6; 3]

  @test_approx_eq_eps state1.logpR state2.logpR ε

  @test_approx_eq_eps support1.logmlik support2.logmlik ε

  # [3; 3; 3; 2; 2; 4] => [3; 3; 1; 2; 2; 4]
  state1 = AminoAcidState(data, [3; 3; 3; 2; 2; 4], priorR, priorC, settings)
  resizestate!(state1, 4, settings)

  state1.R = [3; 3; 3; 2; 2; 4]

  state1.emptycluster[1] = true
  state1.emptycluster[2:4] = false
  state1.cl[1:3] = [2; 3; 4]
  state1.k = 3

  state1.v[2:4] = copy(state1.v[1:3])
  state1.n1s[2:4, :] = copy(state1.n1s[1:3, :])
  state1.unit[4] = copy(state1.unit[3])
  state1.unit[3] = copy(state1.unit[2])
  state1.unit[2] = copy(state1.unit[1])

  support1 = MCMCSupport(state1, priorC)

  state2 = AminoAcidState(data, [3; 3; 1; 2; 2; 4], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  k = state2.k
  i = 3
  hi = state1.R[i]

  initsupportbrwsplit!(k, i, hi, data, priorC, settings, support1, state1)

  support1.lograR = logratiopriorrowbrwsplit(k, state1.v[hi], priorR)

  updatelogmargliki!(priorC, support1)
  updatelogmarglikj!(priorC, support1)

  logmarglikbrwsplit!(state1.cl, state1.k, hi, priorC, support1)

  performbrwsplit!(i, hi, settings, support1, state1)

  @test state1.R == state2.R
  @test state1.k == state2.k

  @test !state1.emptycluster[1]
  @test state1.cl[1] == state2.cl[1]
  @test state1.v[1] == state2.v[1]
  @test state1.n1s[1, :] == state2.n1s[1, :]
  @test state1.unit[1][1:state1.v[1]] == [i]

  @test !state1.emptycluster[2]
  @test state1.cl[2] == state2.cl[2]
  @test state1.v[2] == state2.v[2]
  @test state1.n1s[2, :] == state2.n1s[2, :]
  @test state1.unit[2][1:state1.v[2]] == [4; 5]

  @test !state1.emptycluster[3]
  @test state1.cl[3] == state2.cl[3]
  @test state1.v[3] == state2.v[3]
  @test state1.n1s[3, :] == state2.n1s[3, :]
  @test state1.unit[3][1:state1.v[3]] == [1; 2]

  @test !state1.emptycluster[4]
  @test state1.cl[4] == state2.cl[4]
  @test state1.v[4] == state2.v[4]
  @test state1.n1s[4, :] == state2.n1s[4, :]
  @test state1.unit[4][1:state1.v[4]] == [6]

  @test_approx_eq_eps state1.logpR state2.logpR ε

  @test_approx_eq_eps support1.logmlik support2.logmlik ε

  # [2; 2; 2; 1; 1; 3] => [2; 2; 4; 1; 1; 3]
  # allocate new resources
  state1 = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)
  support1 = MCMCSupport(state1, priorC)

  state2 = AminoAcidState(data, [2; 2; 4; 1; 1; 3], priorR, priorC, settings)
  support2 = MCMCSupport(state2, priorC)

  k = state2.k
  i = 3
  hi = state1.R[i]

  initsupportbrwsplit!(k, i, hi, data, priorC, settings, support1, state1)

  support1.lograR = logratiopriorrowbrwsplit(k, state1.v[hi], priorR)

  updatelogmargliki!(priorC, support1)
  updatelogmarglikj!(priorC, support1)

  logmarglikbrwsplit!(state1.cl, state1.k, hi, priorC, support1)

  performbrwsplit!(i, hi, settings, support1, state1)

  @test state1.R == state2.R
  @test state1.k == state2.k

  @test !state1.emptycluster[1]
  @test state1.cl[1] == state2.cl[1]
  @test state1.v[1] == state2.v[1]
  @test state1.n1s[1, :] == state2.n1s[1, :]
  @test state1.unit[1][1:state1.v[1]] == [4; 5]

  @test !state1.emptycluster[2]
  @test state1.cl[2] == state2.cl[2]
  @test state1.v[2] == state2.v[2]
  @test state1.n1s[2, :] == state2.n1s[2, :]
  @test state1.unit[2][1:state1.v[2]] == [1; 2]

  @test !state1.emptycluster[3]
  @test state1.cl[3] == state2.cl[3]
  @test state1.v[3] == state2.v[3]
  @test state1.n1s[3, :] == state2.n1s[3, :]
  @test state1.unit[3][1:state1.v[3]] == [6]

  @test !state1.emptycluster[4]
  @test state1.cl[4] == state2.cl[4]
  @test state1.v[4] == state2.v[4]
  @test state1.n1s[4, :] == state2.n1s[4, :]
  @test state1.unit[4][1:state1.v[4]] == [i]

  @test_approx_eq_eps state1.logpR state2.logpR ε

  @test_approx_eq_eps support1.logmlik support2.logmlik ε

  nothing
end

test_mcmc_brw()
