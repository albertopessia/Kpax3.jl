# This file is part of Kpax3. License is MIT.

function split!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                mcmcobj::AminoAcidMCMC)
  hi = mcmcobj.R[ij[1]]

  # number of clusters after the split
  k = length(mcmcobj.cl) + 1

  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  initializeSupport!(support, ij[1], ij[2], S, k, data, priorC.logγ, logω,
                     priorC.A, priorC.B)

  # sample a new proportion for cluster 'hi'
  ws = Distributions.rand(settings.distws)

  # logarithm of the product of sequential probabilities
  lq = 0.0

  # temporary / support variables
  g = 0
  u = 0
  M = 0.0
  lcp = zeros(Float64, 2)
  tmp = zeros(Float64, 4)
  wv = zeros(Float64, 2)

  # allocate the neighbours of i and j
  for l in 1:S
    u = neighbours[l]
    lcp[1] = lcp[2] = 0.0

    # compute p(x_{u} | x_{hi,1:(u-1)}) and p(x_{u} | x_{hj,1:(u-1)})
    for b in 1:size(data, 1)
      support.wi.z[1, b] = support.wi.w[1, b] +
                           logcondmarglik(data[b, u], support.ni[b], support.vi,
                                          priorC.A[1, b], priorC.B[1, b])
      support.wi.z[2, b] = support.wi.w[2, b]+
                           logcondmarglik(data[b, u], support.ni[b], support.vi,
                                          priorC.A[2, b], priorC.B[2, b])
      support.wi.z[3, b] = support.wi.w[3, b]+
                           logcondmarglik(data[b, u], support.ni[b], support.vi,
                                          priorC.A[3, b], priorC.B[3, b])
      support.wi.z[4, b] = support.wi.w[4, b] +
                           logcondmarglik(data[b, u], support.ni[b], support.vi,
                                          priorC.A[4, b], priorC.B[4, b])

      tmp[1] = support.wi.z[1, b] - support.wi.c[b]
      tmp[2] = support.wi.z[2, b] - support.wi.c[b]
      tmp[3] = support.wi.z[3, b] - support.wi.c[b]
      tmp[4] = support.wi.z[4, b] - support.wi.c[b]

      M = max(tmp[1], tmp[2], tmp[3], tmp[4])

      lcp[1] += M + log(exp(tmp[1] - M) + exp(tmp[2] - M) +
                        exp(tmp[3] - M) + exp(tmp[4] - M))

      support.wj.z[1, b] = support.wj.w[1, b] +
                           logcondmarglik(data[b, u], support.nj[b], support.vj,
                                          priorC.A[1, b], priorC.B[1, b])
      support.wj.z[2, b] = support.wj.w[2, b]+
                           logcondmarglik(data[b, u], support.nj[b], support.vj,
                                          priorC.A[2, b], priorC.B[2, b])
      support.wj.z[3, b] = support.wj.w[3, b]+
                           logcondmarglik(data[b, u], support.nj[b], support.vj,
                                          priorC.A[3, b], priorC.B[3, b])
      support.wj.z[4, b] = support.wj.w[4, b] +
                           logcondmarglik(data[b, u], support.nj[b], support.vj,
                                          priorC.A[4, b], priorC.B[4, b])

      tmp[1] = support.wj.z[1, b] - support.wj.c[b]
      tmp[2] = support.wj.z[2, b] - support.wj.c[b]
      tmp[3] = support.wj.z[3, b] - support.wj.c[b]
      tmp[4] = support.wj.z[4, b] - support.wj.c[b]

      M = max(tmp[1], tmp[2], tmp[3], tmp[4])

      lcp[2] += M + log(exp(tmp[1] - M) + exp(tmp[2] - M) +
                        exp(tmp[3] - M) + exp(tmp[4] - M))
    end

    # (w * p1) / (w * p1 + (1 - w) * p2) = 1 / (1 + ((1 - w) / w) * (p2 / p1))
    # => e^(-log(1 + e^(log(1 - w) - log(w) + log(p2) - log(p1))))
    wv[1] = exp(-log1p(exp(log(1 - ws) - log(ws) + lcp[2] - lcp[1])))
    wv[2] = 1 - wv[1]

    z = StatsBase.WeightVec(wv)
    g = StatsBase.sample(1:2, z)

    if g == 1
      updateClusteri!(support, u, data)
    else
      updateClusterj!(support, u, data)
    end

    lq += log(values(z)[g])
  end

  # new index of cluster hi
  idx = 0
  for g in mcmcobj.cl
    idx += 1
    support.filledcluster[idx] = true
    support.v[idx] = g != hi ? mcmcobj.v[g] : support.vi
  end

  support.filledcluster[k] = true
  support.v[k] = support.vj

  for b in 1:size(data, 1)
    support.n1s[k, b] = support.nj[b]

    idx = 0
    for g in mcmcobj.cl
      support.n1s[idx += 1, b] = g != hi ? mcmcobj.n1s[g, b] : support.ni[b]
    end
  end

#  if length(support.unit[gi]) < support.vi
#    support.unit[gi] = zeros(Int, support.vi)
#  end
#  copy!(support.unit[gi], 1, support.ui, 1, support.vi)

  cl = collect(1:k)

  rpostpartitioncols!(support.C, cl, support.v, support.n1s,
                      priorC.logγ, logω, priorC.A, priorC.B)

  distwm = Distributions.Beta(settings.parawm + support.vi,
                              settings.parawm + support.vj)

  logsrR = logSplitRatioPriorRow(k, support.vi, support.vj, priorR)
  logprC = logpriorC(support.C, cl, priorC.logγ, logω)

  lidx = 0
  loglik = 0.0
  for b in 1:size(data, 1)
    for g in cl
      lidx = support.C[g, b] + 4 * (b - 1)
      loglik += logmarglik(support.n1s[g, b], support.v[g], priorC.A[lidx],
                           priorC.B[lidx])
    end
  end

  logpocC = logcondpostC(support.C, cl, support.v, support.n1s, priorC.logγ,
                         logω, priorC.A, priorC.B)

  ratio = exp(logsrR +
              logprC - mcmcobj.logprC +
              loglik - mcmcobj.loglik +
              Distributions.logpdf(distwm, ws) + mcmcobj.logpocC -
              Distributions.logpdf(settings.distws, ws) - lq - logpocC)

  if ratio >= 1.0 || (ratio > 0 &&
                      Distributions.rand(Distributions.Bernoulli(ratio)) == 1)
    hj = findfirst(!mcmcobj.filledcluster)

    if hj > 0
      mcmcobj.C[!mcmcobj.emptycluster, :] = C[1:mcmcobj.k, :]
      mcmcobj.C[hj, :] = C[k, :]

      mcmcobj.v[hi] = support.v[gi]
      mcmcobj.unit[hi] = support.unit[gi]
      mcmcobj.n1s[hi, :] = support.n1s[gi, :]

      mcmcobj.v[hj] = support.v[k]
      mcmcobj.unit[hj] = support.unit[k]
      mcmcobj.n1s[hj, :] = support.n1s[k, :]

      mcmcobj.filledcluster[hj] = true
    else
      hj = k

      # reallocate memory
      len = min(settings.maxclust,
                size(data, 2) - length(mcmcobj.filledcluster)) - 1

      mcmcobj.C = vcat(C, zeros(UInt8, len, m))

      mcmcobj.cluster = vcat(cluster,
                             [KCluster(0, zeros(Int, 1), zeros(Float64, 1))
                              for g in 1:len])

      mcmcobj.emptycluster = vcat(mcmcobj.emptycluster, false, trues(len))
    end

    # move units to their new cluster
    mcmcobj.R[support.gj.unit[1:support.gj.v]] = hj

    mcmcobj.k = k

    mcmcobj.logprR = logdPriorRow(n, k, Int[cluster[g].v for g in 1:k], priorR)
    mcmcobj.logprC = logprC
    mcmcobj.loglik = loglik

    mcmcobj.logpocC = logpocC

    priorC.ω = [1.0, 1.0, (k - 1.0) / k, 1.0 / k]
    priorC.logω = logω
  end

  nothing
end
