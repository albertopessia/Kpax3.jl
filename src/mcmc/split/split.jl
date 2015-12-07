# This file is part of Kpax3. License is MIT.

function split!(ij::Array{Int, 1},
                neighbours::Array{Int, 1},
                S::Int,
                data::Array{UInt8, 2},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport
                mcmcobj::AminoAcidMCMC)
  hi = mcmcobj.R[ij[1]]
  hj = findfirst(mcmcobj.emptycluster)

  m, n = size(data)

  # number of clusters after the split
  k = mcmcobj.k + 1

  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  initializeSupport!(support, ij[1], data[:, ij[1]], ij[2], data[:, ij[2]], S,
                     priorC.logγ, logω, priorC.A, priorC.B)

  # sample a new proportion for cluster 'hi'
  ws = Distributions.rand(settings.distws)

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
    wv[1] = exp(-log1p(exp(log(1.0 - ws) - log(ws) + lcp[2] - lcp[1])))
    wv[2] = 1.0 - wv[1]

    z = StatsBase.WeightVec(wv)
    g = StatsBase.sample(1:2, z)

    if g == 1
      updateSupport!(support.gi, support.wi, u, data[:, u])
    else
      updateSupport!(support.gj, support.wj, u, data[:, u])
    end

    lq += log(values(z)[g])
  end

  # sample a new candidate for C
  emptycluster = falses(k)

  C = zeros(UInt8, k, size(data, 1))

  cluster = vcat([KCluster(mcmcobj.cluster[g].v, mcmcobj.cluster[g].unit,
                           mcmcobj.cluster[g].n1s)
                  for g in find(!mcmcobj.emptycluster)],
                 KCluster(support.gj.v, support.gj.unit[1:support.gj.v],
                          support.gj.n1s))

  # new index of cluster hi
  gi = sum(!mcmcobj.emptycluster[1:hi])

  cluster[gi].v = support.gi.v
  if length(cluster[gi].unit) > support.gi.v
    cluster[gi].unit[1:support.gi.v] = support.gi.unit[1:support.gi.v]
  else
    cluster[gi].unit = support.gi.unit[1:support.gi.v]
  end
  cluster[gi].n1s = support.gi.n1s

  rpostpartitioncols!(C, cluster, emptycluster, priorC.logγ, logω, priorC.A,
                      priorC.B)

  distwm = Distributions.Beta(settings.parawm + support.gi.v,
                              settings.parawm + support.gj.v)

  logsrR = logSplitRatioPriorRow(k, support.gi.v, support.gj.v, priorR)
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
              Distributions.logpdf(distwm, ws) + mcmcobj.logpocC -
              Distributions.logpdf(settings.distws, ws) - lq - logpocC)

  if ratio >= 1.0 || (ratio > 0 &&
                      Distributions.rand(Distributions.Bernoulli(ratio)) == 1)
    # move units to their new cluster
    if hj > 0
      mcmcobj.C[!mcmcobj.emptycluster, :] = C[1:mcmcobj.k, :]
      mcmcobj.C[hj, :] = C[k, :]

      mcmcobj.emptycluster[hj] = false

      mcmcobj.cluster[hi].v = cluster[hi].v
      mcmcobj.cluster[hi].unit = cluster[hi].unit
      mcmcobj.cluster[hi].n1s = cluster[hi].n1s

      mcmcobj.cluster[hj].v = cluster[hj].v
      mcmcobj.cluster[hj].unit = cluster[hj].unit
      mcmcobj.cluster[hj].n1s = cluster[hj].n1s
    else
      hj = k

      # reallocate memory
      len = min(settings.maxclust, n - length(mcmcobj.emptycluster)) - 1

      mcmcobj.C = vcat(C, zeros(UInt8, len, m))
      mcmcobj.emptycluster = vcat(mcmcobj.emptycluster, false, trues(len))
      mcmcobj.cluster = vcat(cluster,
                             [KCluster(0, zeros(Int, 1), zeros(Float64, 1))
                              for g in 1:len])
    end

    mcmcobj.R[support.gj.unit[1:support.gj.v]] = hj

    mcmcobj.k = k
    mcmcobj.logprR = logdPriorRow(n, k, Int[cluster[g].v for g in 1:k],
                                  priorR)
    mcmcobj.logprC = logprC
    mcmcobj.loglik = loglik

    priorC.ω = [1.0, 1.0, (k - 1.0) / k, 1.0 / k]
    priorC.logω = logω
  end

  nothing
end
