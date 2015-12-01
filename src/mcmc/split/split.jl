# This file is part of Kpax3. License is MIT.

function addunit!(l::Int,
                  g::Int,
                  x::Array{UInt8, 1},
                  support::Kpax3Support)
  support.smR[l] = g
  support.smV[g] += 1.0

  for j in 1:length(x)
    support.smCO[g, j] += x[j]
  end

  nothing
end

function initializeCO!(support::Kpax3Support,
                       xi::Array{UInt8, 1},
                       xj::Array{UInt8, 1})
  for j in 1:length(xi)
    support.smCO[1, j] = xi[j]
    support.smCO[2, j] = xj[j]
  end

  nothing
end

function initializeMP!(support::Kpax3Support,
                       xi::Array{UInt8, 1},
                       xj::Array{UInt8, 1},
                       logγ::Array{Float64, 1},
                       logω::Array{Float64, 1},
                       A::Array{Float64, 2},
                       B::Array{Float64, 2})
  M = zero(Float64)

  # they will be used to compute conditional probabilities
  for j in 1:length(xi)
    support.smMP[1, 1, j] = logγ[1] + logmarglik(xi[j], 1, A[1, j], B[1, j])
    support.smMP[2, 1, j] = logγ[2] + logmarglik(xi[j], 1, A[2, j], B[2, j])
    support.smMP[3, 1, j] = logγ[3] + logω[3] +
                            logmarglik(xi[j], 1, A[3, j], B[3, j])
    support.smMP[4, 1, j] = logγ[4] + logω[4] +
                            logmarglik(xi[j], 1, A[4, j], B[4, j])

    M = max(support.smMP[1, 1, j], support.smMP[2, 1, j],
            support.smMP[3, 1, j], support.smMP[4, 1, j])

    support.smLC[1, j] = M + log(exp(support.smMP[1, 1, j] - M) +
                                 exp(support.smMP[2, 1, j] - M) +
                                 exp(support.smMP[3, 1, j] - M) +
                                 exp(support.smMP[4, 1, j] - M))

    support.smMP[1, 2, j] = logγ[1] + logmarglik(xj[j], 1, A[1, j], B[1, j])
    support.smMP[2, 2, j] = logγ[2] + logmarglik(xj[j], 1, A[2, j], B[2, j])
    support.smMP[3, 2, j] = logγ[3] + logω[3] +
                            logmarglik(xj[j], 1, A[3, j], B[3, j])
    support.smMP[4, 2, j] = logγ[4] + logω[4] +
                            logmarglik(xj[j], 1, A[4, j], B[4, j])

    M = max(support.smMP[1, 2, j], support.smMP[2, 2, j],
            support.smMP[3, 2, j], support.smMP[4, 2, j])

    support.smLC[2, j] = M + log(exp(support.smMP[1, 2, j] - M) +
                                 exp(support.smMP[2, 2, j] - M) +
                                 exp(support.smMP[3, 2, j] - M) +
                                 exp(support.smMP[4, 2, j] - M))
  end

  nothing
end

function initializeV!(support::Kpax3Support,
                      S::Int)
  # they will be the new sample sizes
  support.smV[1] = one(Float64)
  support.smV[2] = one(Float64)
  support.smV[3] = Float64(S + 2.0)
  nothing
end

function split!(mcmcobj::AminoAcidMCMC,
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                data::Array{UInt8, 2},
                ij::Array{Int, 1},
                neighbours::Array{Int, 1},
                S::Int,
                settings::Kpax3Settings,
                support::Kpax3Support)
  # new cluster labels
  hi = mcmcobj.R[ij[1]]
  hj = findfirst(mcmcobj.emptycluster)

  # number of clusters after the split
  k = mcmcobj.k + 1.0

  initializeV!(support, S)
  initializeCO!(support, data[:, ij[1]], data[:, ij[2]])
  initializeMP!(support, data[:, ij[1]], data[:, ij[2]], priorC.logγ,
                [0.0, 0.0, log(k - 1.0) - log(k), -log(k)], priorC.A, priorC.B)

  # sample a new proportion for cluster 'hi'
  ws = StatsBase.rand(settings.distws)

  # logarithm of sequential probabilities
  lq = zero(Float64)

  # temporary / support variables
  g = zero(Int)
  u = zero(Int)
  lcp = zeros(Float64, 2)
  tmp = zeros(Float64, 4)
  M = zero(Float64)
  wv = zeros(Float64, 2)

  # allocate the neighbours of i and j
  for l in 1:S
    u = neighbours[l]
    lcp[1] = lcp[2] = 0.0

    # compute p(x_{u} | x_{hi,1:(u-1)}) and p(x_{u} | x_{hj,1:(u-1)})
    for j in 1:m
      tmp[1] = support.smMP[1, 1, j] - support.smLC[1, j] +
               logcondmarglik(data[j, u], support.smCO[1, j], support.smV[1],
                              priorC.A[1, j], priorC.B[1, j])
      tmp[2] = support.smMP[2, 1, j] - support.smLC[1, j] +
               logcondmarglik(data[j, u], support.smCO[1, j], support.smV[1],
                              priorC.A[2, j], priorC.B[2, j])
      tmp[3] = support.smMP[3, 1, j] - support.smLC[1, j] +
               logcondmarglik(data[j, u], support.smCO[1, j], support.smV[1],
                              priorC.A[3, j], priorC.B[3, j])
      tmp[4] = support.smMP[4, 1, j] - support.smLC[1, j] +
               logcondmarglik(data[j, u], support.smCO[1, j], support.smV[1],
                              priorC.A[4, j], priorC.B[4, j])

      M = max(tmp[1], tmp[2], tmp[3], tmp[4])

      lcp[1] += M + log(exp(tmp[1] - M) + exp(tmp[2] - M) +
                        exp(tmp[3] - M) + exp(tmp[4] - M))

      tmp[1] = support.smMP[1, 2, j] - support.smLC[2, j] +
               logcondmarglik(data[j, u], support.smCO[2, j], support.smV[2],
                              priorC.A[1, j], priorC.B[1, j])
      tmp[2] = support.smMP[2, 2, j] - support.smLC[2, j] +
               logcondmarglik(data[j, u], support.smCO[2, j], support.smV[2],
                              priorC.A[2, j], priorC.B[2, j])
      tmp[3] = support.smMP[3, 2, j] - support.smLC[2, j] +
               logcondmarglik(data[j, u], support.smCO[2, j], support.smV[2],
                              priorC.A[3, j], priorC.B[3, j])
      tmp[4] = support.smMP[4, 2, j] - support.smLC[2, j] +
               logcondmarglik(data[j, u], support.smCO[2, j], support.smV[2],
                              priorC.A[4, j], priorC.B[4, j])

      M = max(tmp[1], tmp[2], tmp[3], tmp[4])

      lcp[2] += M + log(exp(tmp[1] - M) + exp(tmp[2] - M) +
                        exp(tmp[3] - M) + exp(tmp[4] - M))
    end

    wv[1] = 1.0 / (1.0 + exp(log(1.0 - ws) - log(ws) + lcp[2] - lcp[1]))
    wv[2] = 1.0 - wv[1]

    z = StatsBase.WeightVec(wv, 1.0)
    g = StatsBase.sample(1:2, z)

    addunit!(l, g, data[:, u], support)

    lq += log(values(z)[g])
  end

  distwm = Distributions.Beta(settings.parawm + support.smV[1],
                              settings.parawm + support.smV[2])

  ratio = exp(logSplitRatioPriorRow(priorR, support.smV, k) + )

  if ratio >= 1.0 || rand(Bernoulli(ratio)) == 1
    # move units to their new cluster
  end

  nothing
end
