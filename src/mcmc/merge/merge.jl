# This file is part of Kpax3. License is MIT.

function merge!(ij::Array{Int, 1},
                neighbours::Array{Int, 1},
                S::Int,
                data::Array{UInt8, 2},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                mcmcobj::AminoAcidMCMC)
  hi = mcmcobj.R[ij[1]]
  hj = mcmcobj.R[ij[2]]

  m, n = size(data)

  # number of clusters after the merge
  k = mcmcobj.k - 1

  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  initializeSupport!(support, ij[1], data[:, ij[1]], ij[2], data[:, ij[2]], S,
                     priorC.logγ, logω, priorC.A, priorC.B)

  distwm = Distributions.Beta(settings.parawm + mcmcobj.cluster[hi].v,
                              settings.parawm + mcmcobj.cluster[hj].v)

  # sample a new proportion for cluster 'hi'
  wm = Distributions.rand(distwm)

  # logarithm of the product of sequential probabilities
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
    for b in 1:m
      support.wi.z[1, b] = support.wi.w[1, b] +
                           logcondmarglik(data[b, u], support.gi.n1s[b],
                                          support.gi.v, priorC.A[1, b],
                                          priorC.B[1, b])
      support.wi.z[2, b] = support.wi.w[2, b]+
                           logcondmarglik(data[b, u], support.gi.n1s[b],
                                          support.gi.v, priorC.A[2, b],
                                          priorC.B[2, b])
      support.wi.z[3, b] = support.wi.w[3, b]+
                           logcondmarglik(data[b, u], support.gi.n1s[b],
                                          support.gi.v, priorC.A[3, b],
                                          priorC.B[3, b])
      support.wi.z[4, b] = support.wi.w[4, b] +
                           logcondmarglik(data[b, u], support.gi.n1s[b],
                                          support.gi.v, priorC.A[4, b],
                                          priorC.B[4, b])

      tmp[1] = support.wi.z[1, b] - support.wi.c[b]
      tmp[2] = support.wi.z[2, b] - support.wi.c[b]
      tmp[3] = support.wi.z[3, b] - support.wi.c[b]
      tmp[4] = support.wi.z[4, b] - support.wi.c[b]

      M = max(tmp[1], tmp[2], tmp[3], tmp[4])

      lcp[1] += M + log(exp(tmp[1] - M) + exp(tmp[2] - M) +
                        exp(tmp[3] - M) + exp(tmp[4] - M))

      support.wj.z[1, b] = support.wj.w[1, b] +
                           logcondmarglik(data[b, u], support.gj.n1s[b],
                                          support.gj.v, priorC.A[1, b],
                                          priorC.B[1, b])
      support.wj.z[2, b] = support.wj.w[2, b]+
                           logcondmarglik(data[b, u], support.gj.n1s[b],
                                          support.gj.v, priorC.A[2, b],
                                          priorC.B[2, b])
      support.wj.z[3, b] = support.wj.w[3, b]+
                           logcondmarglik(data[b, u], support.gj.n1s[b],
                                          support.gj.v, priorC.A[3, b],
                                          priorC.B[3, b])
      support.wj.z[4, b] = support.wj.w[4, b] +
                           logcondmarglik(data[b, u], support.gj.n1s[b],
                                          support.gj.v, priorC.A[4, b],
                                          priorC.B[4, b])

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
    wv[1] = -log1p(exp(log(1.0 - wm) - log(wm) + lcp[2] - lcp[1]))
    wv[2] =  log1p(-exp(wv[1]))

    if mcmcobj.R[u] == hi
      updateSupport!(support.gi, support.wi, u, data[:, u])
      lq += wv[1]
    else
      updateSupport!(support.gj, support.wj, u, data[:, u])
      lq += wv[2]
    end
  end

  # sample a new candidate for C
  emptycluster = falses(k)

  ec = falses(length(mcmcobj.emptycluster))
  ec[:] = mcmcobj.emptycluster
  ec[hj] = true

  C = zeros(UInt8, k, size(data, 1))

  cluster = [KCluster(mcmcobj.cluster[g].v, mcmcobj.cluster[g].unit,
                      mcmcobj.cluster[g].n1s) for g in find(!ec)]

  # new index of cluster hi
  gi = sum(!ec[1:hi])

  cluster[gi].v = support.gi.v + support.gj.v
  if length(cluster[gi].unit) > cluster[gi].v
    cluster[gi].unit[1:cluster[gi].v] = vcat(ij, neighbours[1:S])
  else
    cluster[gi].unit = vcat(ij, neighbours[1:S])
  end
  cluster[gi].n1s = support.gi.n1s

  rpostpartitioncols!(C, cluster, emptycluster, priorC.logγ, logω, priorC.A,
                      priorC.B)

  logsrR = logMergeRatioPriorRow(k, support.gi.v, support.gj.v, priorR)
  logprC = logpriorC(C, emptycluster, priorC.logγ, logω)

  loglik = zero(Float64)
  for g in 1:k
    linearidx = [(C[g, b] + 4 * (b - 1))::Int for b in 1:m]
    loglik += sum(logmarglik(cluster[g].n1s, cluster[g].v, priorC.A[linearidx],
                             priorC.B[linearidx]))
  end

  logpocC = logcondpostC(C, cluster, emptycluster, priorC.logγ, logω, priorC.A,
                         priorC.B)

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
