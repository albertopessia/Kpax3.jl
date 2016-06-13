# This file is part of Kpax3. License is MIT.

function test_mcmc_merge_init()
  # merge cluster 2 and cluster 3
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test.bin"

  settings = KSettings(ifile, ofile, maxclust=2, maxunit=1)

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

  R = [13; 13; 13; 42; 42; 76]
  k = length(unique(R)) - 1

  ij = [4; 6]
  S = 1
  u = 5

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportmerge!(ij, k, data, priorC, support)

  wi = zeros(Float64, 4, m)
  for col in 1:m
    wi[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                            priorC.B[1, col])

    wi[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                            priorC.B[2, col])

    wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                            priorC.B[3, col])

    wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                            priorC.B[4, col])
  end

  ci = Float64[log(sum(exp(support.wi.w[:, b]))) for b in 1:m]

  wj = zeros(Float64, 4, m)
  for col in 1:m
    wj[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                            priorC.B[1, col])

    wj[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                            priorC.B[2, col])

    wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                            priorC.B[3, col])

    wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                            priorC.B[4, col])
  end

  cj = Float64[log(sum(exp(support.wj.w[:, b]))) for b in 1:m]

  @test support.vi == 1
  @test support.ni == float(data[:, ij[1]])
  @test support.ui == [ij[1]; 0; 0; 0; 0; 0]

  @test maximum(abs(support.wi.w - wi)) <= ε
  @test maximum(abs(support.wi.c - ci)) <= ε
  @test support.wi.z == zeros(Float64, 4, m)

  @test support.vj == 1
  @test support.nj == float(data[:, ij[2]])
  @test support.uj == [ij[2]; 0; 0; 0; 0; 0]

  @test maximum(abs(support.wj.w - wj)) <= ε
  @test maximum(abs(support.wj.c - cj)) <= ε
  @test support.wj.z == zeros(Float64, 4, m)

  nothing
end

test_mcmc_merge_init()

function test_mcmc_merge_updatei()
  # move the first unit (u) to cluster 2 (test inverse split operator)
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test.bin"

  settings = KSettings(ifile, ofile, maxclust=2, maxunit=1)

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

  R = [13; 13; 13; 42; 42; 76]
  k = length(unique(R)) - 1

  ij = [4; 6]
  S = 1
  u = 5

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportmerge!(ij, k, data, priorC, support)

  lcp = zeros(Float64, 2)
  for b in 1:m
    lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, support)
    lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, support)
  end
  updateclusteri!(u, data, support)

  wi = zeros(Float64, 4, m)
  for col in 1:m
    wi[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                            priorC.B[1, col]) +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[1, col], priorC.B[1, col])

    wi[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                            priorC.B[2, col]) +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[2, col], priorC.B[2, col])

    wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                            priorC.B[3, col]) +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[3, col], priorC.B[3, col])

    wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                            priorC.B[4, col]) +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[4, col], priorC.B[4, col])
  end

  zi = copy(wi)

  ci = Float64[log(sum(exp(support.wi.w[:, b]))) for b in 1:m]

  wj = zeros(Float64, 4, m)
  zj = zeros(Float64, 4, m)
  for col in 1:m
    wj[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                            priorC.B[1, col])

    wj[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                            priorC.B[2, col])

    wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                            priorC.B[3, col])

    wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                            priorC.B[4, col])

    zj[1, col] = wj[1, col] +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[1, col], priorC.B[1, col])

    zj[2, col] = wj[2, col] +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[2, col], priorC.B[2, col])

    zj[3, col] = wj[3, col] +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[3, col], priorC.B[3, col])

    zj[4, col] = wj[4, col] +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[4, col], priorC.B[4, col])

  end

  cj = Float64[log(sum(exp(support.wj.w[:, b]))) for b in 1:m]

  @test support.vi == 2
  @test support.ni == float(data[:, ij[1]]) + float(data[:, u])
  @test support.ui == [ij[1]; u; 0; 0; 0; 0]

  @test maximum(abs(support.wi.w - wi)) <= ε
  @test maximum(abs(support.wi.c - ci)) <= ε
  @test maximum(abs(support.wi.z - zi)) <= ε

  @test support.vj == 1
  @test support.nj == float(data[:, ij[2]])
  @test support.uj == [ij[2]; 0; 0; 0; 0; 0]

  @test maximum(abs(support.wj.w - wj)) <= ε
  @test maximum(abs(support.wj.c - cj)) <= ε
  @test maximum(abs(support.wj.z - zj)) <= ε

  nothing
end

test_mcmc_merge_updatei()

function test_mcmc_merge_updatej()
  # move the first unit (u) to cluster 2 (test inverse split operator)
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test.bin"

  settings = KSettings(ifile, ofile, maxclust=2, maxunit=1)

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

  R = [13; 13; 13; 42; 42; 76]
  k = length(unique(R)) - 1

  ij = [4; 6]
  S = 1
  u = 5

  state = AminoAcidState(data, R, priorR, priorC, settings)
  support = MCMCSupport(state, priorC)

  initsupportmerge!(ij, k, data, priorC, support)

  lcp = zeros(Float64, 2)
  for b in 1:m
    lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, support)
    lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, support)
  end
  updateclusterj!(u, data, support)

  wi = zeros(Float64, 4, m)
  zi = zeros(Float64, 4, m)
  for col in 1:m
    wi[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                            priorC.B[1, col])

    wi[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                            priorC.B[2, col])

    wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                            priorC.B[3, col])

    wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                            priorC.B[4, col])

    zi[1, col] = wi[1, col] +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[1, col], priorC.B[1, col])

    zi[2, col] = wi[2, col] +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[2, col], priorC.B[2, col])

    zi[3, col] = wi[3, col] +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[3, col], priorC.B[3, col])

    zi[4, col] = wi[4, col] +
                 logcondmarglik(data[col, u], data[col, ij[1]], 1,
                                priorC.A[4, col], priorC.B[4, col])
  end

  ci = Float64[log(sum(exp(support.wi.w[:, b]))) for b in 1:m]

  wj = zeros(Float64, 4, m)
  for col in 1:m
    wj[1, col] = priorC.logγ[1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                            priorC.B[1, col]) +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[1, col], priorC.B[1, col])

    wj[2, col] = priorC.logγ[2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                            priorC.B[2, col]) +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[2, col], priorC.B[2, col])

    wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                            priorC.B[3, col]) +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[3, col], priorC.B[3, col])

    wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
                 logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                            priorC.B[4, col]) +
                 logcondmarglik(data[col, u], data[col, ij[2]], 1,
                                priorC.A[4, col], priorC.B[4, col])
  end
  zj = copy(wj)

  cj = Float64[log(sum(exp(support.wj.w[:, b]))) for b in 1:m]

  @test support.vi == 1
  @test support.ni == float(data[:, ij[1]])
  @test support.ui == [ij[1]; 0; 0; 0; 0; 0]

  @test maximum(abs(support.wi.w - wi)) <= ε
  @test maximum(abs(support.wi.c - ci)) <= ε
  @test maximum(abs(support.wi.z - zi)) <= ε

  @test support.vj == 2
  @test support.nj == float(data[:, ij[2]]) + float(data[:, u])
  @test support.uj == [ij[2]; u; 0; 0; 0; 0]

  @test maximum(abs(support.wj.w - wj)) <= ε
  @test maximum(abs(support.wj.c - cj)) <= ε
  @test maximum(abs(support.wj.z - zj)) <= ε

  nothing
end

test_mcmc_merge_updatej()
