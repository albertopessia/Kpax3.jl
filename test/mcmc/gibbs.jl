# This file is part of Kpax3. License is MIT.

function test_mcmc_gibbs_init()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=3, maxunit=1)

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

  # remove unit 5 from cluster 2 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 5
  hi = 2

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportgibbs!(hi, priorR, state, support)

  @test support.k == 2
  @test support.cl[1:support.k] == [1; 3]
  @test support.t[1:support.k] == [clusterweight(state.v[1], priorR);
                                   clusterweight(state.v[3], priorR)]

  nothing
end

test_mcmc_gibbs_init()

function test_mcmc_gibbs_move()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=4, maxunit=1)

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

  # remove unit 5 from cluster 2 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 5
  hi = 2
  li = 3
  lj = 4

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  support.ni = vec(float(data[:, 4]))
  support.vi = 1

  initsupportgibbs!(hi, priorR, state, support)

  support.t[li] = clusterweight(support.vi, priorR)
  support.t[lj] = clusterweight(1, state.k, priorR)

  g = 0
  l = 0
  y = 0.0
  v = 0
  z = zeros(Float64, 2)
  for b in 1:m
    support.lpi[1, b] = priorC.logγ[1]
    support.lpi[2, b] = priorC.logγ[2]
    support.lpi[3, b] = priorC.logγ[3]
    support.lpi[4, b] = priorC.logγ[3]

    g = 1
    l = 1
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    @test support.lq[1, l, b] == logmarglik(y, v, priorC.A[1, b],
                                            priorC.B[1, b])
    @test support.lq[2, l, b] == logmarglik(y, v, priorC.A[2, b],
                                            priorC.B[2, b])
    @test support.lq[3, l, b] == logmarglik(y, v, priorC.A[3, b],
                                            priorC.B[3, b])
    @test support.lq[4, l, b] == logmarglik(y, v, priorC.A[4, b],
                                            priorC.B[4, b])

    g = 3
    l = 2
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    @test support.lq[1, l, b] == logmarglik(y, v, priorC.A[1, b],
                                            priorC.B[1, b])
    @test support.lq[2, l, b] == logmarglik(y, v, priorC.A[2, b],
                                            priorC.B[2, b])
    @test support.lq[3, l, b] == logmarglik(y, v, priorC.A[3, b],
                                            priorC.B[3, b])
    @test support.lq[4, l, b] == logmarglik(y, v, priorC.A[4, b],
                                            priorC.B[4, b])

    @test_approx_eq_eps support.lpi[1, b] (priorC.logγ[1] +
                                           support.lp[1, 1, b] +
                                           support.lp[1, 3, b]) ε

    @test_approx_eq_eps support.lpi[2, b] (priorC.logγ[2] +
                                           support.lp[2, 1, b] +
                                           support.lp[2, 3, b]) ε

    z[1] = log(exp(priorC.logω[state.k][1] + support.lp[3, 1, b]) +
               exp(priorC.logω[state.k][2] + support.lp[4, 1, b]))
    z[2] = log(exp(priorC.logω[state.k][1] + support.lp[3, 3, b]) +
               exp(priorC.logω[state.k][2] + support.lp[4, 3, b]))

    @test_approx_eq_eps support.lpi[3, b] (priorC.logγ[3] + z[1] + z[2]) ε

    z[1] = support.lq[1, 1, b] + support.lq[1, li, b] - support.lp[1, 1, b]
    z[2] = support.lq[1, 2, b] + support.lq[1, li, b] - support.lp[1, 3, b]
    @test_approx_eq_eps support.lr[1, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[1, 2, b] z[2] ε

    z[1] = support.lq[2, 1, b] + support.lq[2, li, b] - support.lp[2, 1, b]
    z[2] = support.lq[2, 2, b] + support.lq[2, li, b] - support.lp[2, 3, b]
    @test_approx_eq_eps support.lr[2, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[2, 2, b] z[2] ε

    z[1] = log(exp(priorC.logω[state.k][1] + support.lq[3, 1, b]) +
               exp(priorC.logω[state.k][2] + support.lq[4, 1, b])) -
           log(exp(priorC.logω[state.k][1] + support.lp[3, 1, b]) +
               exp(priorC.logω[state.k][2] + support.lp[4, 1, b]))

    z[2] = log(exp(priorC.logω[state.k][1] + support.lq[3, 2, b]) +
               exp(priorC.logω[state.k][2] + support.lq[4, 2, b])) -
           log(exp(priorC.logω[state.k][1] + support.lp[3, 3, b]) +
               exp(priorC.logω[state.k][2] + support.lp[4, 3, b]))

    @test_approx_eq_eps support.lr[3, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[3, 2, b] z[2] ε

    gibbsupdateclusteri!(state.k, b, li, priorC, support)

    y = support.ni[b]
    v = support.vi

    @test support.lq[1, li, b] == logmarglik(y, v, priorC.A[1, b],
                                             priorC.B[1, b])
    @test support.lq[2, li, b] == logmarglik(y, v, priorC.A[2, b],
                                             priorC.B[2, b])
    @test support.lq[3, li, b] == logmarglik(y, v, priorC.A[3, b],
                                             priorC.B[3, b])
    @test support.lq[4, li, b] == logmarglik(y, v, priorC.A[4, b],
                                             priorC.B[4, b])

    z[1] = priorC.logω[state.k][1] + support.lq[3, li, b]
    z[2] = priorC.logω[state.k][2] + support.lq[4, li, b]

    @test_approx_eq_eps support.lpj[1, b] log(exp(z[1]) + exp(z[2])) ε

    gibbsupdateclusterj!(data[b, i], state.k + 1, b, li, lj, priorC,
                             support)

    y = support.ni[b]
    v = 1

    z[1] = priorC.logω[state.k + 1][1] + logmarglik(y, v, priorC.A[3, b],
                                                    priorC.B[3, b])
    z[2] = priorC.logω[state.k + 1][2] + logmarglik(y, v, priorC.A[4, b],
                                                    priorC.B[4, b])

    @test_approx_eq_eps support.lpj[2, b] log(exp(z[1]) + exp(z[2])) ε

    @test support.lq[1, lj, b] == logmarglik(data[b, i], 1, priorC.A[1, b],
                                             priorC.B[1, b])
    @test support.lq[2, lj, b] == logmarglik(data[b, i], 1, priorC.A[2, b],
                                             priorC.B[2, b])
    @test support.lq[3, lj, b] == logmarglik(data[b, i], 1, priorC.A[3, b],
                                             priorC.B[3, b])
    @test support.lq[4, lj, b] == logmarglik(data[b, i], 1, priorC.A[4, b],
                                             priorC.B[4, b])

    z[1] = priorC.logω[state.k + 1][1] + support.lq[3, lj, b]
    z[2] = priorC.logω[state.k + 1][2] + support.lq[4, lj, b]

    @test_approx_eq_eps support.lpj[3, b] log(exp(z[1]) + exp(z[2])) ε
  end

  gibbscomputeprobi!(li, support)

  for b in 1:support.m
    gibbscomputeprobg!(b, 1, li, support)
    gibbscomputeprobg!(b, 2, li, support)
    gibbscomputeprobj!(b, li, lj, support)
  end

  st = AminoAcidState(data, [1; 1; 1; 2; 1; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(state.v[1], priorR) + su.logmlik

  @test_approx_eq_eps support.t[1] t ε

  st = AminoAcidState(data, [1; 1; 1; 2; 3; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(state.v[3], priorR) + su.logmlik

  @test_approx_eq_eps support.t[2] t ε

  st = AminoAcidState(data, [1; 1; 1; 2; 2; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(state.v[2] - 1, priorR) + su.logmlik

  @test_approx_eq_eps support.t[li] t ε

  st = AminoAcidState(data, [1; 1; 1; 2; 4; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(1, state.k, priorR) + su.logmlik

  @test_approx_eq_eps support.t[lj] t ε

  logc = gibbscomputenormconst(lj, support)

  @test_approx_eq_eps logc log(sum(exp(support.t[1:4]))) ε

  nothing
end

test_mcmc_gibbs_move()

function test_mcmc_gibbs_merge()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=3, maxunit=1)

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

  # remove unit 6 from cluster 3 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 6
  hi = 3
  li = 3
  lj = 3

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportgibbs!(hi, priorR, state, support)

  support.t[lj] = clusterweight(1, state.k - 1, priorR)

  g = 0
  l = 0
  y = 0.0
  v = 0
  z = zeros(Float64, 2)
  for b in 1:m
    support.lpi[1, b] = priorC.logγ[1]
    support.lpi[2, b] = priorC.logγ[2]
    support.lpi[3, b] = priorC.logγ[3]

    g = 1
    l = 1
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmerge!(y, v, b, l, g, support.k, priorC, support)

    @test support.lq[1, l, b] == logmarglik(y, v, priorC.A[1, b],
                                            priorC.B[1, b])
    @test support.lq[2, l, b] == logmarglik(y, v, priorC.A[2, b],
                                            priorC.B[2, b])
    @test support.lq[3, l, b] == logmarglik(y, v, priorC.A[3, b],
                                            priorC.B[3, b])
    @test support.lq[4, l, b] == logmarglik(y, v, priorC.A[4, b],
                                            priorC.B[4, b])

    g = 2
    l = 2
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmerge!(y, v, b, l, g, support.k, priorC, support)

    @test support.lq[1, l, b] == logmarglik(y, v, priorC.A[1, b],
                                            priorC.B[1, b])
    @test support.lq[2, l, b] == logmarglik(y, v, priorC.A[2, b],
                                            priorC.B[2, b])
    @test support.lq[3, l, b] == logmarglik(y, v, priorC.A[3, b],
                                            priorC.B[3, b])
    @test support.lq[4, l, b] == logmarglik(y, v, priorC.A[4, b],
                                            priorC.B[4, b])

    @test_approx_eq_eps support.lpi[1, b] (priorC.logγ[1] +
                                           support.lp[1, 1, b] +
                                           support.lp[1, 2, b]) ε

    @test_approx_eq_eps support.lpi[2, b] (priorC.logγ[2] +
                                           support.lp[2, 1, b] +
                                           support.lp[2, 2, b]) ε

    z[1] = log(exp(priorC.logω[support.k][1] + support.lp[3, 1, b]) +
               exp(priorC.logω[support.k][2] + support.lp[4, 1, b]))
    z[2] = log(exp(priorC.logω[support.k][1] + support.lp[3, 2, b]) +
               exp(priorC.logω[support.k][2] + support.lp[4, 2, b]))

    @test_approx_eq_eps support.lpi[3, b] (priorC.logγ[3] + z[1] + z[2]) ε

    z[1] = support.lq[1, 1, b] + support.lq[1, li, b] - support.lp[1, 1, b]
    z[2] = support.lq[1, 2, b] + support.lq[1, li, b] - support.lp[1, 2, b]
    @test_approx_eq_eps support.lr[1, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[1, 2, b] z[2] ε

    z[1] = support.lq[2, 1, b] + support.lq[2, li, b] - support.lp[2, 1, b]
    z[2] = support.lq[2, 2, b] + support.lq[2, li, b] - support.lp[2, 2, b]
    @test_approx_eq_eps support.lr[2, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[2, 2, b] z[2] ε

    z[1] = log(exp(priorC.logω[support.k][1] + support.lq[3, 1, b]) +
               exp(priorC.logω[support.k][2] + support.lq[4, 1, b])) +
           support.lpj[1, b] -
           log(exp(priorC.logω[support.k][1] + support.lp[3, 1, b]) +
               exp(priorC.logω[support.k][2] + support.lp[4, 1, b]))
    z[2] = log(exp(priorC.logω[support.k][1] + support.lq[3, 2, b]) +
               exp(priorC.logω[support.k][2] + support.lq[4, 2, b])) +
           support.lpj[1, b] -
           log(exp(priorC.logω[support.k][1] + support.lp[3, 2, b]) +
               exp(priorC.logω[support.k][2] + support.lp[4, 2, b]))

    @test_approx_eq_eps support.lr[3, 1, b] z[1] ε
    @test_approx_eq_eps support.lr[3, 2, b] z[2] ε
  end

  gibbscomputeprobi!(lj, support)

  for b in 1:support.m
    gibbscomputeprobg!(b, 1, support)
    gibbscomputeprobg!(b, 2, support)
  end

  st = AminoAcidState(data, [1; 1; 1; 2; 2; 1], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(state.v[1], priorR) + su.logmlik

  @test_approx_eq_eps support.t[1] t ε

  st = AminoAcidState(data, [1; 1; 1; 2; 2; 2], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(state.v[2], priorR) + su.logmlik

  @test_approx_eq_eps support.t[2] t ε

  st = AminoAcidState(data, [1; 1; 1; 2; 2; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)
  t = clusterweight(1, state.k - 1, priorR) + su.logmlik

  @test_approx_eq_eps support.t[lj] t ε

  logc = gibbscomputenormconst(lj, support)

  @test_approx_eq_eps logc log(sum(exp(support.t[1:3]))) ε

  nothing
end

test_mcmc_gibbs_merge()

function test_mcmc_gibbs_sample_cluster()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=4, maxunit=1)

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

  # remove unit 5 from cluster 2 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  support.t[1] = log(0.001) + log(0.4)
  support.t[2] = log(0.001) + log(0.3)
  support.t[3] = log(0.001) + log(0.2)
  support.t[4] = log(0.001) + log(0.1)

  logc = log(0.001)

  N = 1000000
  g = 0
  counter = zeros(Float64, 6)
  for t in 1:N
    g = gibbssamplecluster(logc, 4, support)
    counter[g] += 1
  end

  counter /= N

  @test maximum(abs(counter - [0.4; 0.3; 0.2; 0.1; 0.0; 0.0])) <= 0.005

  nothing
end

test_mcmc_gibbs_sample_cluster()

function test_mcmc_gibbs_perform_move()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=4, maxunit=1)

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

  # remove unit 5 from cluster 2 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 5
  hi = 2
  li = 3
  lj = 4

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  support.ni = vec(float(data[:, 4]))
  support.vi = 1

  initsupportgibbs!(hi, priorR, state, support)

  support.t[li] = clusterweight(support.vi, priorR)
  support.t[lj] = clusterweight(1, state.k, priorR)

  g = 0
  l = 0
  y = 0.0
  v = 0
  z = zeros(Float64, 2)
  for b in 1:m
    support.lpi[1, b] = priorC.logγ[1]
    support.lpi[2, b] = priorC.logγ[2]
    support.lpi[3, b] = priorC.logγ[3]
    support.lpi[4, b] = priorC.logγ[3]

    g = 1
    l = 1
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    g = 3
    l = 2
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    gibbsupdateclusteri!(state.k, b, li, priorC, support)
    gibbsupdateclusterj!(data[b, i], state.k + 1, b, li, lj, priorC,
                             support)
  end

  gibbscomputeprobi!(li, support)

  for b in 1:support.m
    gibbscomputeprobg!(b, 1, li, support)
    gibbscomputeprobg!(b, 2, li, support)
    gibbscomputeprobj!(b, li, lj, support)
  end

  # R = [1; 1; 1; 2; 3; 3]
  l = 2

  gibbsmove!(i, hi, support.cl[l], li, l, data, priorR, support, state)

  st = AminoAcidState(data, [1; 1; 1; 2; 3; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)

  @test state.R == st.R
  @test state.k == st.k

  @test !state.emptycluster[1]
  @test state.cl[1] == st.cl[1]
  @test state.v[1] == st.v[1]
  @test state.n1s[1, :] == st.n1s[1, :]
  @test state.unit[1][1:state.v[1]] == [1; 2; 3]

  @test !state.emptycluster[2]
  @test state.cl[2] == st.cl[2]
  @test state.v[2] == st.v[2]
  @test state.n1s[2, :] == st.n1s[2, :]
  @test state.unit[2][1:state.v[2]] == [4]

  @test !state.emptycluster[3]
  @test state.cl[3] == st.cl[3]
  @test state.v[3] == st.v[3]
  @test state.n1s[3, :] == st.n1s[3, :]
  @test state.unit[3][1:state.v[3]] == [6; 5]

  @test_approx_eq_eps state.logpR st.logpR ε

  @test_approx_eq_eps support.logmlik su.logmlik ε

  @test support.lp == su.lp

  nothing
end

test_mcmc_gibbs_perform_move()

function test_mcmc_gibbs_perform_merge()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile, maxclust=3, maxunit=1)

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

  # remove unit 6 from cluster 3 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 6
  hi = 3
  li = 3
  lj = 3

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportgibbs!(hi, priorR, state, support)

  support.t[lj] = clusterweight(1, state.k - 1, priorR)

  g = 0
  l = 0
  y = 0.0
  v = 0
  z = zeros(Float64, 2)
  for b in 1:m
    support.lpi[1, b] = priorC.logγ[1]
    support.lpi[2, b] = priorC.logγ[2]
    support.lpi[3, b] = priorC.logγ[3]

    g = 1
    l = 1
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmerge!(y, v, b, l, g, support.k, priorC, support)

    g = 2
    l = 2
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmerge!(y, v, b, l, g, support.k, priorC, support)
  end

  gibbscomputeprobi!(lj, support)

  for b in 1:support.m
    gibbscomputeprobg!(b, 1, support)
    gibbscomputeprobg!(b, 2, support)
  end

  # R = [1; 1; 1; 2; 2; 1]
  l = 1

  gibbsmerge!(i, hi, support.cl[l], l, data, priorR, support, state)

  st = AminoAcidState(data, [1; 1; 1; 2; 2; 1], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)

  @test state.R == st.R
  @test state.k == st.k

  @test !state.emptycluster[1]
  @test state.cl[1] == st.cl[1]
  @test state.v[1] == st.v[1]
  @test state.n1s[1, :] == st.n1s[1, :]
  @test state.unit[1][1:state.v[1]] == [1; 2; 3; 6]

  @test !state.emptycluster[2]
  @test state.cl[2] == st.cl[2]
  @test state.v[2] == st.v[2]
  @test state.n1s[2, :] == st.n1s[2, :]
  @test state.unit[2][1:state.v[2]] == [4; 5]

  @test state.emptycluster[3]

  @test_approx_eq_eps state.logpR st.logpR ε

  @test_approx_eq_eps support.logmlik su.logmlik ε

  @test support.lp[:, 1:2, :] == su.lp[:, 1:2, :]

  nothing
end

test_mcmc_gibbs_perform_merge()

function test_mcmc_gibbs_perform_split()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  # state should be resized
  settings = KSettings(ifile, ofile, maxclust=3, maxunit=1)

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

  # remove unit 5 from cluster 2 and put it into another cluster
  R = [1; 1; 1; 2; 2; 3]

  i = 5
  hi = 2
  li = 3
  lj = 4

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  resizesupport!(support, state.k + 1)

  support.ni = vec(float(data[:, 4]))
  support.vi = 1

  initsupportgibbs!(hi, priorR, state, support)

  support.t[li] = clusterweight(support.vi, priorR)
  support.t[lj] = clusterweight(1, state.k, priorR)

  g = 0
  l = 0
  y = 0.0
  v = 0
  z = zeros(Float64, 2)
  for b in 1:m
    support.lpi[1, b] = priorC.logγ[1]
    support.lpi[2, b] = priorC.logγ[2]
    support.lpi[3, b] = priorC.logγ[3]
    support.lpi[4, b] = priorC.logγ[3]

    g = 1
    l = 1
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    g = 3
    l = 2
    y = state.n1s[g, b] + float(data[b, i])
    v = state.v[g] + 1

    gibbsupdateclustergmove!(y, v, b, l, g, state.k, priorC, support)

    gibbsupdateclusteri!(state.k, b, li, priorC, support)
    gibbsupdateclusterj!(data[b, i], state.k + 1, b, li, lj, priorC,
                             support)
  end

  gibbscomputeprobi!(li, support)

  for b in 1:support.m
    gibbscomputeprobg!(b, 1, li, support)
    gibbscomputeprobg!(b, 2, li, support)
    gibbscomputeprobj!(b, li, lj, support)
  end

  # R = [1; 1; 1; 2; 4; 3]
  l = 4

  gibbssplit!(i, hi, li, lj, data, priorR, settings, support, state)

  st = AminoAcidState(data, [1; 1; 1; 2; 4; 3], priorR, priorC, settings)
  su = MCMCSupport(st, priorC)

  @test state.R == st.R
  @test state.k == st.k

  @test !state.emptycluster[1]
  @test state.cl[1] == st.cl[1]
  @test state.v[1] == st.v[1]
  @test state.n1s[1, :] == st.n1s[1, :]
  @test state.unit[1][1:state.v[1]] == [1; 2; 3]

  @test !state.emptycluster[2]
  @test state.cl[2] == st.cl[2]
  @test state.v[2] == st.v[2]
  @test state.n1s[2, :] == st.n1s[2, :]
  @test state.unit[2][1:state.v[2]] == [4]

  @test !state.emptycluster[3]
  @test state.cl[3] == st.cl[3]
  @test state.v[3] == st.v[3]
  @test state.n1s[3, :] == st.n1s[3, :]
  @test state.unit[3][1:state.v[3]] == [6]

  @test !state.emptycluster[4]
  @test state.cl[4] == st.cl[4]
  @test state.v[4] == st.v[4]
  @test state.n1s[4, :] == st.n1s[4, :]
  @test state.unit[4][1:state.v[4]] == [i]

  @test_approx_eq_eps state.logpR st.logpR ε

  @test_approx_eq_eps support.logmlik su.logmlik ε

  @test support.lp[:, 1:4, :] == su.lp

  nothing
end

test_mcmc_gibbs_perform_split()
