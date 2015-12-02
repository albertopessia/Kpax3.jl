# This file is part of Kpax3. License is MIT.

function initializeSupport!(support::KSupport,
                            i::Int,
                            xi::Array{UInt8, 1},
                            j::Int,
                            xj::Array{UInt8, 1},
                            S::Int,
                            logγ::Array{Float64, 1},
                            logω::Array{Float64, 1},
                            A::Array{Float64, 2},
                            B::Array{Float64, 2})
  M = zero(Float64)

  support.gi.v = one(Int)
  if length(support.gi.unit) < S + 1
    support.gi.unit = zeros(Int, S + 1)
  end
  support.gi.unit[1] = i
  support.gi.ll = 0.0

  support.gj.v = one(Int)
  if length(support.gj.unit) < S + 1
    support.gj.unit = zeros(Int, S + 1)
  end
  support.gj.unit[1] = j
  support.gj.ll = 0.0

  for b in 1:length(xi)
    support.gi.n1s[b] = Float64(xi[b])

    support.wi.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(xi[b], 1, A[1, b], B[1, b])
    support.wi.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(xi[b], 1, A[2, b], B[2, b])
    support.wi.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(xi[b], 1, A[3, b], B[3, b])
    support.wi.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(xi[b], 1, A[4, b], B[4, b])

    M = max(support.wi.w[1, b], support.wi.w[2, b],
            support.wi.w[3, b], support.wi.w[4, b])

    support.wi.c[b] = M + log(exp(support.wi.w[1, b] - M) +
                              exp(support.wi.w[2, b] - M) +
                              exp(support.wi.w[3, b] - M) +
                              exp(support.wi.w[4, b] - M))

    support.gj.n1s[b] = Float64(xj[b])

    support.wj.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(xj[b], 1, A[1, b], B[1, b])
    support.wj.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(xj[b], 1, A[2, b], B[2, b])
    support.wj.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(xj[b], 1, A[3, b], B[3, b])
    support.wj.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(xj[b], 1, A[4, b], B[4, b])

    M = max(support.wj.w[1, b], support.wj.w[2, b],
            support.wj.w[3, b], support.wj.w[4, b])

    support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                              exp(support.wj.w[2, b] - M) +
                              exp(support.wj.w[3, b] - M) +
                              exp(support.wj.w[4, b] - M))
  end

  nothing
end

function updateSupport!(h::KCluster,
                        w::KWeight,
                        u::Int,
                        x::Array{UInt8, 1})
  M = zero(Float64)

  h.v += 1
  h.unit[h.v] = u

  for b in 1:length(x)
    h.n1s[b] += Float64(x[b])

    w.w[1, b] = w.z[1, b]
    w.w[2, b] = w.z[2, b]
    w.w[3, b] = w.z[3, b]
    w.w[4, b] = w.z[4, b]

    M = max(w.w[1, b], w.w[2, b], w.w[3, b], w.w[4, b])

    w.c[b] = M + log(exp(w.w[1, b] - M) + exp(w.w[2, b] - M) +
                     exp(w.w[3, b] - M) + exp(w.w[4, b] - M))
  end

  nothing
end

function split!(ij::Array{Int, 1},
                neighbours::Array{Int, 1},
                S::Int,
                data::Array{UInt8, 2},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport
                mcmcobj::AminoAcidMCMC)
  # number of clusters after the split
  k = mcmcobj.k + 1

  logω = [0.0, 0.0, log(k - 1.0) - log(k), -log(k)]

  initializeSupport!(support, ij[1], data[:, ij[1]], ij[2], data[:, ij[2]], S,
                     priorC.logγ, logω, priorC.A, priorC.B)

  # sample a new proportion for cluster 'hi'
  ws = StatsBase.rand(settings.distws)

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
  C = zeros(UInt8, k, size(data, 1))
  cluster = [mcmcobj.cluster[1:mcmcobj.k]; support.gj]
  cluster[mcmcobj.R[ij[1]]] = support.gi
  emptycluster = falses(k)
  rpostpartitioncols!(C, cluster, emptycluster, priorC.logγ, logω, priorC.A,
                      priorC.B)

  distwm = Distributions.Beta(settings.parawm + support.gi.v,
                              settings.parawm + support.gj.v)

  logsrR = logSplitRatioPriorRow(k, support.gi.v, support.gj.v, priorR)
  logprC = logpriorC(C, emptycluster, priorC.logγ, logω)
  loglik = zero(Float64)
  for g in 1:k
    linearidx = [(C[g, b] + 4 * (b - 1))::Int for b in 1:m]
    cluster[g].ll = sum(logmarglik(cluster[g].n1s, cluster[g].v,
                                   priorC.A[linearidx], priorC.B[linearidx]))
    loglik += cluster[g].ll
  end

  ratio = exp(logsrR +
              logprC - mcmcobj.logprC +
              loglik - mcmcobj.loglik +
              Distributions.logpdf(distwm, ws) -
              Distributions.logpdf(settings.distws, ws) -
              lq -
              logcondpostC(C, cluster, emptycluster, priorC.logγ, logω,
                           priorC.A, priorC.B))

  if ratio >= 1.0 || rand(Bernoulli(ratio)) == 1
    # move units to their new cluster
    if size(mcmcobj.C, 1) > k
    else
    end
  end

  nothing
end
