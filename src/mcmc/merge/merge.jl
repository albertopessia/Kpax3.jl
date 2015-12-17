# This file is part of Kpax3. License is MIT.

function merge!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                mcmcobj::AminoAcidMCMC)
  hi = mcmcobj.R[ij[1]]
  hj = mcmcobj.R[ij[2]]

  # number of clusters after the merge
  k = length(mcmcobj.cl) - 1

  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  initializeSupport!(support, ij[1], ij[2], S, k, data, priorC.logγ, logω,
                     priorC.A, priorC.B)

  distwm = Distributions.Beta(settings.parawm + mcmcobj.v[hi],
                              settings.parawm + mcmcobj.v[hj])

  # sample a new proportion for cluster 'hi'
  wm = Distributions.rand(distwm)

  # logarithm of the product of sequential probabilities
  lq = 0.0

  # temporary / support variables
  u = 0
  M = 0.0
  lcp = zeros(Float64, 2)
  tmp = zeros(Float64, 4)
  z = 0.0
  p = 0.0

  vi = mcmcobj.v[hi] + mcmcobj.v[hj]
  ni = zeros(Float64, size(data, 1))
  for b in 1:size(data, 1)
    ni[b] = mcmcobj.n1s[hi, b] + mcmcobj.n1s[hj, b]
  end

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
    z = -log1p(exp(log(1.0 - wm) - log(wm) + lcp[2] - lcp[1]))
    p = exp(z)

    if mcmcobj.R[u] == hi
      updateClusteri!(support, u, data)
      lq += z
    else
      updateClusterj!(support, u, data)
      lq += log1p(-p)
    end
  end

  logprC, logpocC = mergesimc!(k, hi, logω, priorC, support, mcmcobj)

  ratio = exp(logsrR +
              logprC - mcmcobj.logprC +
              loglik - mcmcobj.loglik +
              Distributions.logpdf(settings.distws, wm) + lq + mcmcobj.logpocC -
              Distributions.logpdf(distwm, wm) - logpocC)

  if ratio >= 1.0 || (ratio > 0.0 &&
                      Distributions.rand(Distributions.Bernoulli(ratio)) == 1)
    mcmcobj.logpocC = logpocC

    mcmcobj.loglik = loglik
    mcmcobj.logprC = logprC
    mcmcobj.logprR = logdPriorRow(n, k, Int[cluster[g].v for g in 1:k], priorR)

    mcmcobj.k = k
    mcmcobj.emptycluster[hj] = true
    mcmcobj.cluster[hi].v = cluster[hi].v
    mcmcobj.cluster[hi].unit = cluster[hi].unit
    mcmcobj.cluster[hi].n1s = cluster[hi].n1s

    mcmcobj.C[!mcmcobj.emptycluster, :] = C
    mcmcobj.R[support.gj.unit[1:support.gj.v]] = hi

    priorC.ω = [1.0, 1.0, (k - 1.0) / k, 1.0 / k]
    priorC.logω = logω
  end

  nothing
end
